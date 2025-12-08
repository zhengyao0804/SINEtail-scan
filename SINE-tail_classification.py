#!/usr/bin/env python3
"""
改进的SINE尾部结构模式分析工具

该脚本用于分析猪基因组中SINE转座子尾部的复杂结构模式。
根据SINE家族信息选择截断点匹配策略:
- SINEA1-11家族: 取序列倒数70bp进行匹配
- SINEB1-6家族: 取序列倒数50bp进行匹配
- SINEC1-8家族: 取序列倒数50bp进行匹配

对所有家族，统一使用以下截断点序列顺序尝试匹配:
1. GCGGCCC (完全匹配)
2. GCGCGGCC (允许1个错配)
3. GCGGCCT (允许1个错配)
4. GCGGCCC (允许1个错配)
5. GCGGCCA (允许1个错配)
6. GCGGCCG (允许1个错配)
7. GCATGGCCC (完全匹配)
8. GCATGGCCC (允许1个错配)
9. GCGCGGCCC (完全匹配)
10. GCGCGGCCC (允许1个错配)
11. GCGTGGCCC (允许1个错配)
12. TGTGTGGC (完全匹配)
13. TGTGTGGC (允许1个错配)

尾巴结构分类:
- A-rich: 连续A占尾巴序列70%及以上
- (AAAAC)n, (AAAC)n, (AAC)n, (AC)n, AC-composite
- (AAAAT)n, (AAAT)n, (AAT)n, (AT)n, AT-composite
- (AAAAG)n, (AAAG)n, (AAG)n, (AG)n, AG-composite
- 其他类型

用法:
    python sine_tail_pattern.py <fasta文件> --bed_file <BED文件> [--output PREFIX]
"""

import re
import os
import argparse
import numpy as np
from collections import Counter, defaultdict
from itertools import product
from Bio import pairwise2
from Bio.Seq import Seq


def read_family_info(bed_file):
    """
    从BED文件中读取SINE家族信息

    参数:
    bed_file (str): BED文件路径(7列格式)

    返回:
    dict: 序列ID到SINE家族的映射
    """
    family_map = {}
    if not bed_file or not os.path.exists(bed_file):
        return family_map

    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 7:
                # 第4列为Family信息，第7列为序列ID
                seq_id = fields[6]
                family = fields[3]
                family_map[seq_id] = family

    return family_map


def read_fasta(file_path):
    """
    读取FASTA文件中的序列

    参数:
    file_path (str): FASTA文件路径

    返回:
    dict: 序列ID到序列的映射字典
    """
    sequences = {}

    with open(file_path, 'r') as f:
        current_id = ""
        current_seq = ""

        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq and current_id:
                    sequences[current_id] = current_seq

                current_id = line[1:]  # 移除'>'前缀
                current_seq = ""
            else:
                current_seq += line.upper()  # 转换为大写

        if current_seq and current_id:
            sequences[current_id] = current_seq

    return sequences


def fuzzy_search(text, pattern, max_mismatches=1):
    """
    在文本中搜索模式,允许有限数量的错配

    参数:
    text (str): 要搜索的文本
    pattern (str): 要查找的模式
    max_mismatches (int): 允许的最大错配数

    返回:
    list: 匹配位置列表,每个元素为(开始位置,结束位置)
    """
    pattern_len = len(pattern)
    matches = []

    # 在文本中滑动窗口搜索
    for i in range(len(text) - pattern_len + 1):
        window = text[i:i + pattern_len]
        mismatches = sum(1 for a, b in zip(window, pattern) if a != b)

        if mismatches <= max_mismatches:
            matches.append((i, i + pattern_len))

    return matches


def normalize_family_name(family):
    """
    标准化SINE家族名称，处理SscSINEA1/SINEA格式

    参数:
    family (str): 原始家族名称

    返回:
    str: 标准化后的家族名称 (SINEA1, SINEB2等)
    """
    # 处理SscSINEX#/SINEX格式
    if '/' in family:
        # 分割，取左半部分 (SscSINEA1)
        left_part = family.split('/')[0]

        # 如果以Ssc开头，移除前缀
        if left_part.startswith('Ssc'):
            return left_part[3:]  # 返回SINEA1部分
        return left_part

    # 已经是标准格式或其他格式
    return family


def determine_family_group(family):
    """
    确定SINE家族所属的主要组

    参数:
    family (str): 标准化后的家族名称

    返回:
    str: 家族组名 (SINEA, SINEB, SINEC 或 unknown)
    """
    if family.startswith("SINEA"):
        return "SINEA"
    elif family.startswith("SINEB"):
        return "SINEB"
    elif family.startswith("SINEC"):
        return "SINEC"
    else:
        return "unknown"


def extract_tails_by_family_new(sequences, family_map):
    """
    根据SINE家族信息选择截断点策略提取尾部序列 (新逻辑)

    参数:
    sequences (dict): 序列ID到序列的映射字典
    family_map (dict): 序列ID到SINE家族的映射

    返回:
    tuple: (尾部字典, 使用的截断点信息字典)
    """
    tails = {}
    motif_used = {}  # 记录每个序列使用的截断点

    # 定义通用的截断点顺序策略
    cutoff_strategies = [
        # 新增的截断点序列（放在最前面，优先级更高）
        ("AAGAAATAGC", 0, "AAGAAATAGC_exact"),
        ("AAGAAATAGC", 1, "AAGAAATAGC_fuzzy"),
        ("TAGAAAAGGC", 0, "TAGAAAAGGC_exact"),
        ("TAGAAAAGGC", 1, "TAGAAAAGGC_fuzzy"),
        ("TAGAAAAGAC", 0, "TAGAAAAGAC_exact"),
        ("TAGAAAAGAC", 1, "TAGAAAAGAC_fuzzy"),

        # 在TAGAAAAGAC后面、GCGGCCC前面添加新的截断点序列
        ("TAAAAAGAC", 0, "TAAAAAGAC_exact"),
        ("TAAAAAGAC", 1, "TAAAAAGAC_fuzzy"),

        # 原有的截断点序列
        ("GCGGCCC", 0, "GCGGCCC_exact"),
        ("GCGCGGCC", 1, "GCGCGGCC_fuzzy"),
        ("GCGGCCT", 1, "GCGGCCT_fuzzy"),
        ("GCGGCCC", 1, "GCGGCCC_fuzzy"),
        ("GCGGCCA", 1, "GCGGCCA_fuzzy"),
        ("GCGGCCG", 1, "GCGGCCG_fuzzy"),
        ("GCATGGCCC", 0, "GCATGGCCC_exact"),
        ("GCATGGCCC", 1, "GCATGGCCC_fuzzy"),
        ("GCGCGGCCC", 0, "GCGCGGCCC_exact"),
        ("GCGCGGCCC", 1, "GCGCGGCCC_fuzzy"),
        ("GCGTGGCCC", 1, "GCGTGGCCC_fuzzy"),
        ("TGTGTGGC", 0, "TGTGTGGC_exact"),
        ("TGTGTGGC", 1, "TGTGTGGC_fuzzy")
    ]

    # 按家族组分组序列并标准化家族名称
    family_group_sequences = defaultdict(dict)
    for seq_id, seq in sequences.items():
        original_family = family_map.get(seq_id, "unknown")
        normalized_family = normalize_family_name(original_family)
        family_group = determine_family_group(normalized_family)
        family_group_sequences[family_group][seq_id] = (seq, normalized_family, original_family)

    # 为每个家族组应用相应的策略
    for family_group, group_seqs in family_group_sequences.items():
        # 确定搜索区域
        if family_group == "SINEA":
            tail_region_size = 70  # 取序列末尾70bp
        elif family_group in ["SINEB", "SINEC"]:
            tail_region_size = 50  # 取序列末尾50bp
        else:
            tail_region_size = 50  # 默认取末尾50bp

        # 处理每个序列
        for seq_id, (seq, normalized_family, original_family) in group_seqs.items():
            # 提取匹配区域
            seq_len = len(seq)
            if seq_len <= tail_region_size:
                search_region = seq  # 序列太短，使用整个序列
            else:
                search_region = seq[-tail_region_size:]  # 取末尾指定长度的区域

            # 尝试所有截断点策略
            found_match = False
            for pattern, max_mismatches, strategy_name in cutoff_strategies:
                if found_match:
                    break

                matches = fuzzy_search(search_region, pattern, max_mismatches=max_mismatches)
                if matches:
                    # 第一个匹配的位置（相对于search_region）
                    match_position = matches[0][0]
                    match_end = matches[0][1]

                    # 调整为相对于原始序列的位置
                    absolute_match_end = (seq_len - tail_region_size) + match_end

                    # 提取尾部
                    tail = seq[absolute_match_end:]
                    if tail:
                        tails[seq_id] = tail
                        motif_used[seq_id] = f"{family_group}_{strategy_name}"
                        found_match = True

            # 如果没有匹配任何截断点，使用最后30bp
            if not found_match:
                tail = seq[-30:] if len(seq) > 30 else seq
                tails[seq_id] = tail
                motif_used[seq_id] = f"{family_group}_last30bp"

    # 返回尾部字典和使用的截断点信息
    return tails, motif_used

def find_poly_a_with_mismatch(sequence, min_a_count=5, max_mismatches=1):
    """
    查找允许错配的polyA区域
    
    例如：AAAAACAAAAA 会被识别为一个连续的polyA区域（10个A，1个错配C）
    
    参数:
    sequence (str): 输入序列
    min_a_count (int): 最小A的数量（默认5）
    max_mismatches (int): 允许的最大错配数（默认1）
    
    返回:
    list: 包含匹配信息的列表，每个元素模拟re.Match对象
    """
    class PolyAMatch:
        """模拟re.Match对象，用于兼容原有代码"""
        def __init__(self, start, end, group_str):
            self._start = start
            self._end = end
            self._group = group_str
        
        def start(self):
            return self._start
        
        def end(self):
            return self._end
        
        def group(self):
            return self._group
    
    regions = []
    i = 0
    
    while i < len(sequence):
        if sequence[i] == 'A':
            start = i
            a_count = 1
            mismatch_count = 0
            j = i + 1
            
            # 扫描可能的polyA区域
            while j < len(sequence):
                if sequence[j] == 'A':
                    a_count += 1
                    j += 1
                elif mismatch_count < max_mismatches:
                    # 检查错配后面是否还有A（避免在序列末尾无意义地使用错配）
                    if j + 1 < len(sequence):
                        # 向前看，检查后续是否有A
                        has_following_a = False
                        for k in range(j + 1, min(j + 3, len(sequence))):
                            if sequence[k] == 'A':
                                has_following_a = True
                                break
                        
                        if has_following_a:
                            mismatch_count += 1
                            j += 1
                        else:
                            break
                    else:
                        break
                else:
                    break
            
            # 检查A的数量是否满足最小要求
            if a_count >= min_a_count:
                end = j
                match_obj = PolyAMatch(start, end, sequence[start:end])
                regions.append(match_obj)
                i = end
            else:
                i += 1
        else:
            i += 1
    
    return regions


def analyze_new_tail_structure_advanced(sequence):
    """
    根据高级分类标准分析尾部序列结构（优化版：允许polyA错配）

    判断顺序:
    1. 首先进行AC/AT/AG-composite复合结构判断
    2. 然后进行AC,AT,AG系列的单一模式判断（其中(AC)n/(AT)n/(AG)n优先级最低）
    3. 然后进行"A-rich"结构判断（允许1个错配）
    4. 最后归为其他类

    参数:
    sequence (str): 尾部序列

    返回:
    dict: 分析结果
    """
    result = {}

    # 定义基本模式组和它们的正则表达式
    pattern_groups = {
        "AC": [
            ("(AAAAC)n", r'(AAAAC){2,}'),
            ("(AAAC)n", r'(AAAC){2,}'),
            ("(AAC)n", r'(AAC){2,}'),
            ("(AC)n", r'(AC){3,}')
        ],
        "AT": [
            ("(AAAAT)n", r'(AAAAT){2,}'),
            ("(AAAT)n", r'(AAAT){2,}'),
            ("(AAT)n", r'(AAT){2,}'),
            ("(AT)n", r'(AT){3,}')
        ],
        "AG": [
            ("(AAAAG)n", r'(AAAAG){2,}'),
            ("(AAAG)n", r'(AAAG){2,}'),
            ("(AAG)n", r'(AAG){2,}'),
            ("(AG)n", r'(AG){3,}')
        ]
    }

    # 计算A的含量和连续A区域覆盖率（使用允许错配的新方法）
    a_count = sequence.count('A')
    a_percentage = (a_count / len(sequence)) * 100 if len(sequence) > 0 else 0

    # ===== 修改部分：使用允许错配的polyA区域查找 =====
    # 原来的代码：poly_a_regions = list(re.finditer(r'A{5,}', sequence))
    poly_a_regions = find_poly_a_with_mismatch(sequence, min_a_count=5, max_mismatches=1)
    # ================================================
    
    total_poly_a_length = sum(len(match.group()) for match in poly_a_regions)
    poly_a_coverage = (total_poly_a_length / len(sequence)) * 100 if len(sequence) > 0 else 0

    result['a_percentage'] = a_percentage
    result['poly_a_coverage'] = poly_a_coverage
    result['poly_a_regions_count'] = len(poly_a_regions)
    result['longest_poly_a'] = max([len(match.group()) for match in poly_a_regions], default=0) if poly_a_regions else 0

    # 在每个组内搜索模式
    group_matches = {}
    pattern_details = {}  # 存储每种模式的详细信息

    for group_name, patterns in pattern_groups.items():
        group_matches[group_name] = []

        for pattern_name, regex in patterns:
            matches = list(re.finditer(regex, sequence))

            if matches:
                repeat_unit = pattern_name.strip("()n")
                total_repeats = sum(len(match.group()) // len(repeat_unit) for match in matches)

                pattern_details[pattern_name] = {
                    'count': len(matches),
                    'total_repeats': total_repeats,
                    'longest_match': max(len(match.group()) for match in matches),
                    'coverage': sum(len(match.group()) for match in matches) / len(sequence) * 100 if len(
                        sequence) > 0 else 0
                }

            for match in matches:
                group_matches[group_name].append({
                    'pattern': pattern_name,
                    'start': match.start(),
                    'end': match.end(),
                    'length': match.end() - match.start(),
                    'sequence': match.group()
                })

    # 新的分类逻辑顺序

    # 1. 首先查找复合模式
    composite_found = False
    for group_name, matches in group_matches.items():
        if matches:
            unique_patterns = set(match['pattern'] for match in matches)
            if len(unique_patterns) >= 2:
                result['structure_type'] = f"{group_name}-composite"
                result['composite_patterns'] = list(unique_patterns)
                result['pattern_details'] = pattern_details
                composite_found = True
                break

    # 2. 如果没有找到复合模式，查找单一模式 (注意优先级)
    if not composite_found:
        single_pattern_found = False

        # 定义模式的优先级 (数字越小优先级越高)
        priority_map = {
            "(AAAAC)n": 1, "(AAAC)n": 2, "(AAC)n": 3, "(AC)n": 4,
            "(AAAAT)n": 1, "(AAAT)n": 2, "(AAT)n": 3, "(AT)n": 4,
            "(AAAAG)n": 1, "(AAAG)n": 2, "(AAG)n": 3, "(AG)n": 4
        }

        # 收集所有匹配的模式
        all_patterns = []
        for group, matches in group_matches.items():
            for match in matches:
                all_patterns.append((match['pattern'], priority_map.get(match['pattern'], 999), match['start']))

        if all_patterns:
            # 根据优先级和位置排序
            # 1. 按优先级排序（优先级高的在前）
            # 2. 如果优先级相同，按照在序列中的位置排序（靠前的优先）
            all_patterns.sort(key=lambda x: (x[1], x[2]))

            best_pattern = all_patterns[0][0]
            result['structure_type'] = best_pattern
            result['pattern_details'] = pattern_details
            single_pattern_found = True

        # 3. 如果没有找到单一模式，检查是否为A-rich（使用新的poly_a_coverage）
        if not single_pattern_found:
            if poly_a_coverage >= 70:
                result['structure_type'] = "A-rich"
            else:
                # 4. 如果都不符合，归为其他类型
                result['structure_type'] = "Other"

    # 存储所有统计信息
    result['detected_patterns'] = pattern_details

    return result



def new_simplified_classification_advanced(sequence):
    """
    使用高级分类标准对尾部序列进行分类

    参数:
    sequence (str): 尾部序列

    返回:
    str: 分类结果
    """
    # 分析结构特征
    analysis_result = analyze_new_tail_structure_advanced(sequence)

    # 直接返回确定的结构类型
    return analysis_result['structure_type']


def find_repeats(sequence, min_unit_size=1, max_unit_size=5, min_repeats=2):
    """
    在序列中查找重复单元

    参数:
    sequence (str): 输入序列
    min_unit_size (int): 最小重复单元大小
    max_unit_size (int): 最大重复单元大小
    min_repeats (int): 最小重复次数

    返回:
    list: 重复单元及其位置的列表
    """
    repeats = []

    # 对不同大小的重复单元进行搜索
    for unit_size in range(min_unit_size, max_unit_size + 1):
        i = 0
        while i <= len(sequence) - unit_size * min_repeats:
            unit = sequence[i:i + unit_size]

            # 查找连续重复
            j = i + unit_size
            repeat_count = 1

            while j <= len(sequence) - unit_size and sequence[j:j + unit_size] == unit:
                repeat_count += 1
                j += unit_size

            if repeat_count >= min_repeats:
                repeats.append({
                    'unit': unit,
                    'size': unit_size,
                    'count': repeat_count,
                    'start': i,
                    'end': j,
                    'total_length': unit_size * repeat_count
                })
                i = j  # 跳过已找到的重复
            else:
                i += 1

    # 按照总长度排序
    repeats.sort(key=lambda x: x['total_length'], reverse=True)

    return repeats


def find_tandem_repeats(sequence, min_period=2, max_period=6):
    """
    查找序列中的串联重复

    参数:
    sequence (str): 输入序列
    min_period (int): 最小重复周期
    max_period (int): 最大重复周期

    返回:
    list: 串联重复列表
    """
    repeats = []

    for period in range(min_period, max_period + 1):
        i = 0
        while i <= len(sequence) - 2 * period:
            unit = sequence[i:i + period]
            if sequence[i + period:i + 2 * period] == unit:
                # 找到重复，继续查找更多
                j = i + 2 * period
                while j + period <= len(sequence) and sequence[j:j + period] == unit:
                    j += period

                repeat_count = (j - i) // period
                repeats.append({
                    'unit': unit,
                    'period': period,
                    'count': repeat_count,
                    'start': i,
                    'end': j,
                    'sequence': sequence[i:j]
                })

                i = j  # 跳过已找到的重复
            else:
                i += 1

    # 按重复长度排序
    repeats.sort(key=lambda x: x['count'] * x['period'], reverse=True)

    return repeats


def analyze_dinucleotide_frequencies(sequence):
    """
    分析序列中的二核苷酸频率

    参数:
    sequence (str): 输入序列

    返回:
    dict: 二核苷酸频率
    """
    dinucleotides = {}
    for i in range(len(sequence) - 1):
        dinuc = sequence[i:i + 2]
        dinucleotides[dinuc] = dinucleotides.get(dinuc, 0) + 1

    total = sum(dinucleotides.values())
    frequencies = {dinuc: count / total for dinuc, count in dinucleotides.items()} if total > 0 else {}

    # 找出最高频的二核苷酸
    sorted_freqs = sorted(frequencies.items(), key=lambda x: x[1], reverse=True)

    return {
        'counts': dinucleotides,
        'frequencies': frequencies,
        'top': sorted_freqs[:5] if sorted_freqs else []
    }


def analyze_sequence_complexity(sequence):
    """
    分析序列复杂度

    参数:
    sequence (str): 输入序列

    返回:
    dict: 序列复杂度分析结果
    """
    if not sequence:
        return {'complexity': 0, 'entropy': 0}

    # 计算序列熵
    nucleotide_freqs = Counter(sequence)
    total = len(sequence)
    entropy = 0

    for nuc, count in nucleotide_freqs.items():
        p = count / total
        entropy -= p * np.log2(p)

    # 正则化熵(最大熵为2位，对应于4个核苷酸均等出现)
    max_entropy = np.log2(min(4, len(nucleotide_freqs)))
    normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0

    # 计算语言复杂度(使用k-mer种类数除以可能的总数)
    k = 2  # 使用二核苷酸
    observed_kmers = set()

    for i in range(len(sequence) - k + 1):
        observed_kmers.add(sequence[i:i + k])

    # 计算复杂度
    unique_nucleotides = set(sequence)
    possible_kmers = min(len(sequence) - k + 1, 4 ** k)
    complexity = len(observed_kmers) / possible_kmers if possible_kmers > 0 else 0

    return {
        'entropy': entropy,
        'normalized_entropy': normalized_entropy,
        'complexity': complexity,
        'unique_nucleotides': len(unique_nucleotides),
        'unique_dinucleotides': len(observed_kmers)
    }


def analyze_tail(sequence):
    """
    全面分析尾部序列

    参数:
    sequence (str): 尾部序列

    返回:
    dict: 分析结果
    """
    results = {
        'sequence': sequence,
        'length': len(sequence)
    }

    # 使用高级结构分析方法获取详细信息
    structure_analysis = analyze_new_tail_structure_advanced(sequence)

    # 将结构分析的结果合并到主结果中
    for key, value in structure_analysis.items():
        results[key] = value

    # 设置分类结果
    results['classification'] = structure_analysis['structure_type']

    # 查找重复单元
    results['repeats'] = find_repeats(sequence)

    # 查找串联重复
    results['tandem_repeats'] = find_tandem_repeats(sequence)

    # 分析二核苷酸频率
    results['dinucleotide_freq'] = analyze_dinucleotide_frequencies(sequence)

    # 序列复杂度分析
    results['complexity'] = analyze_sequence_complexity(sequence)

    return results


def calculate_sequence_similarity(seq1, seq2):
    """
    计算两个序列之间的相似度

    参数:
    seq1 (str): 第一个序列
    seq2 (str): 第二个序列

    返回:
    float: 相似度得分(0-1)
    """
    # 使用Biopython的全局对齐
    alignments = pairwise2.align.globalxx(seq1, seq2)

    if not alignments:
        return 0

    # 获取最佳对齐
    best_alignment = alignments[0]
    alignment_score = best_alignment.score

    # 归一化得分
    max_length = max(len(seq1), len(seq2))
    similarity = alignment_score / max_length if max_length > 0 else 0

    return similarity


def cluster_sequences(sequences, threshold=0.7):
    """
    根据序列相似性对序列进行聚类

    参数:
    sequences (dict): 序列ID到尾部序列的映射字典
    threshold (float): 聚类相似度阈值

    返回:
    list: 聚类结果
    """
    if not sequences:
        return []

    # 初始化聚类
    clusters = []
    seq_ids = list(sequences.keys())

    # 创建第一个聚类
    first_id = seq_ids[0]
    first_cluster = {
        'representative': first_id,
        'members': [first_id],
        'sequence': sequences[first_id]
    }
    clusters.append(first_cluster)

    # 处理剩余序列
    for seq_id in seq_ids[1:]:
        seq = sequences[seq_id]
        best_similarity = 0
        best_cluster = None

        # 检查与现有聚类的相似性
        for i, cluster in enumerate(clusters):
            similarity = calculate_sequence_similarity(seq, cluster['sequence'])
            if similarity > best_similarity:
                best_similarity = similarity
                best_cluster = i

        # 如果相似性超过阈值，添加到最佳聚类
        if best_similarity >= threshold and best_cluster is not None:
            clusters[best_cluster]['members'].append(seq_id)
        else:
            # 否则创建新的聚类
            new_cluster = {
                'representative': seq_id,
                'members': [seq_id],
                'sequence': seq
            }
            clusters.append(new_cluster)

    # 计算聚类统计
    for cluster in clusters:
        cluster['size'] = len(cluster['members'])
        cluster['percentage'] = (cluster['size'] / len(seq_ids)) * 100

    # 按大小排序
    clusters.sort(key=lambda x: x['size'], reverse=True)

    return clusters


def write_detailed_report(seq_id, analysis_results, output_dir):
    """
    为单个序列写入详细报告

    参数:
    seq_id (str): 序列ID
    analysis_results (dict): 分析结果
    output_dir (str): 输出目录
    """
    os.makedirs(output_dir, exist_ok=True)

    # 处理文件名中的非法字符
    safe_id = seq_id.replace(":", "_").replace(" ", "_").replace("/", "_").replace("\\", "_")
    safe_id = safe_id.replace("*", "_").replace("?", "_").replace("\"", "_").replace("<", "_")
    safe_id = safe_id.replace(">", "_").replace("|", "_")

    with open(os.path.join(output_dir, f"{safe_id}.txt"), 'w') as f:
        f.write(f"SINE尾部结构分析详细报告 - {seq_id}\n")
        f.write("=" * 50 + "\n\n")

        # 基本信息
        f.write(f"序列: {analysis_results['sequence']}\n")
        f.write(f"长度: {analysis_results['length']} bp\n")
        f.write(f"分类: {analysis_results['classification']}\n")

        # 添加家族信息
        if 'family' in analysis_results and analysis_results['family']:
            f.write(f"家族: {analysis_results['family']}\n")

        # 添加使用的截断点信息
        if 'cutoff_used' in analysis_results:
            f.write(f"使用的截断点: {analysis_results['cutoff_used']}\n")
        f.write("\n")

        # 重复分析
        repeats = analysis_results['repeats']
        f.write("重复单元分析:\n")
        if repeats:
            for i, repeat in enumerate(repeats[:3]):  # 只显示前3个
                f.write(f"  重复 #{i + 1}: '{repeat['unit']}' (大小: {repeat['size']} bp)\n")
                f.write(f"    出现次数: {repeat['count']}\n")
                f.write(f"    总长度: {repeat['total_length']} bp\n")
        else:
            f.write("  未检测到显著重复单元\n")
        f.write("\n")

        # 串联重复
        tandem_repeats = analysis_results['tandem_repeats']
        f.write("串联重复分析:\n")
        if tandem_repeats:
            for i, repeat in enumerate(tandem_repeats[:3]):  # 只显示前3个
                f.write(f"  串联重复 #{i + 1}: '{repeat['unit']}' (周期: {repeat['period']} bp)\n")
                f.write(f"    重复次数: {repeat['count']}\n")
        else:
            f.write("  未检测到串联重复\n")
        f.write("\n")

        # 序列复杂度
        complexity = analysis_results['complexity']
        f.write("序列复杂度分析:\n")
        f.write(f"  熵: {complexity['entropy']:.4f}\n")
        f.write(f"  归一化熵: {complexity['normalized_entropy']:.4f}\n")
        f.write(f"  复杂度: {complexity['complexity']:.4f}\n")


def write_summary_report(all_results, clusters, output_file):
    """
    写入汇总报告，包含更多统计信息

    参数:
    all_results (dict): 所有序列的分析结果
    clusters (list): 聚类结果
    output_file (str): 输出文件
    """
    with open(output_file, 'w') as f:
        f.write("SINE尾部结构分析汇总报告\n")
        f.write("======================\n\n")

        # 基本统计
        f.write(f"分析的序列总数: {len(all_results)}\n\n")

        # 统计家族信息
        family_counts = Counter(
            [results.get('family', '') for results in all_results.values() if results.get('family', '')])

        if family_counts:
            f.write("SINE家族分布:\n")
            for family, count in family_counts.most_common():
                percentage = (count / len(all_results)) * 100
                f.write(f"  {family}: {count} ({percentage:.2f}%)\n")
            f.write("\n")

        # 截断点使用情况统计
        cutoff_counts = Counter(
            [results.get('cutoff_used', '') for results in all_results.values() if results.get('cutoff_used', '')])

        f.write("截断点使用情况:\n")
        for cutoff, count in cutoff_counts.most_common():
            percentage = (count / len(all_results)) * 100
            f.write(f"  {cutoff}: {count} ({percentage:.2f}%)\n")
        f.write("\n")

        # 分类统计
        classifications = [results['classification'] for results in all_results.values()]
        class_counts = Counter(classifications)

        f.write("分类分布:\n")
        for cls, count in class_counts.most_common():
            percentage = (count / len(all_results)) * 100
            f.write(f"  {cls}: {count} ({percentage:.2f}%)\n")
        f.write("\n")

        # 按家族分类统计
        if family_counts:
            f.write("按家族划分的分类分布:\n")
            for family, _ in family_counts.most_common():
                f.write(f"  {family}家族:\n")
                family_seq_ids = [seq_id for seq_id, results in all_results.items()
                                  if results.get('family') == family]

                if family_seq_ids:
                    family_classifications = [all_results[seq_id]['classification'] for seq_id in family_seq_ids]
                    family_class_counts = Counter(family_classifications)

                    for cls, count in family_class_counts.most_common():
                        percentage = (count / len(family_seq_ids)) * 100
                        f.write(f"    {cls}: {count} ({percentage:.2f}%)\n")
            f.write("\n")

        # 长度统计
        lengths = [results['length'] for results in all_results.values()]
        avg_length = sum(lengths) / len(lengths) if lengths else 0

        f.write("长度统计:\n")
        f.write(f"  平均长度: {avg_length:.2f} bp\n")
        f.write(f"  中位长度: {sorted(lengths)[len(lengths) // 2] if lengths else 0} bp\n")
        f.write(f"  最小长度: {min(lengths) if lengths else 0} bp\n")
        f.write(f"  最大长度: {max(lengths) if lengths else 0} bp\n\n")

        # polyA相关统计
        f.write("PolyA相关统计:\n")
        a_percentages = []
        poly_a_coverages = []
        longest_poly_as = []
        poly_a_regions_counts = []

        for results in all_results.values():
            if 'a_percentage' in results:
                a_percentages.append(results['a_percentage'])
            if 'poly_a_coverage' in results:
                poly_a_coverages.append(results['poly_a_coverage'])
            if 'longest_poly_a' in results:
                longest_poly_as.append(results['longest_poly_a'])
            if 'poly_a_regions_count' in results:
                poly_a_regions_counts.append(results['poly_a_regions_count'])

        if a_percentages:
            f.write(f"  平均A含量: {sum(a_percentages) / len(a_percentages):.2f}%\n")
            f.write(f"  最高A含量: {max(a_percentages):.2f}%\n")

        if poly_a_coverages:
            f.write(f"  平均polyA覆盖率: {sum(poly_a_coverages) / len(poly_a_coverages):.2f}%\n")
            f.write(f"  最高polyA覆盖率: {max(poly_a_coverages):.2f}%\n")

        if longest_poly_as:
            f.write(f"  平均最长连续A长度: {sum(longest_poly_as) / len(longest_poly_as):.2f} bp\n")
            f.write(f"  最长连续A长度: {max(longest_poly_as)} bp\n")

        if poly_a_regions_counts:
            f.write(f"  平均polyA区域数: {sum(poly_a_regions_counts) / len(poly_a_regions_counts):.2f}\n")
            f.write(f"  最多polyA区域数: {max(poly_a_regions_counts)}\n")

        f.write("\n")

        # 重复模式统计
        f.write("重复模式统计:\n")

        # 单一重复模式统计
        single_patterns = ["(AAAAC)n", "(AAAC)n", "(AAC)n", "(AC)n",
                           "(AAAAT)n", "(AAAT)n", "(AAT)n", "(AT)n",
                           "(AAAAG)n", "(AAAG)n", "(AAG)n", "(AG)n"]

        # 统计各种模式的出现情况
        pattern_stats = {}
        for pattern in single_patterns:
            pattern_stats[pattern] = {
                'sequences': 0,
                'total_repeats': 0,
                'max_repeats': 0,
                'avg_repeats': 0,
                'coverage': []
            }

        for results in all_results.values():
            if 'detected_patterns' in results:
                for pattern, details in results['detected_patterns'].items():
                    if pattern in pattern_stats:
                        pattern_stats[pattern]['sequences'] += 1
                        pattern_stats[pattern]['total_repeats'] += details.get('total_repeats', 0)
                        pattern_stats[pattern]['max_repeats'] = max(
                            pattern_stats[pattern]['max_repeats'],
                            details.get('count', 0)
                        )
                        pattern_stats[pattern]['coverage'].append(details.get('coverage', 0))

        # 计算平均值并输出
        for pattern, stats in pattern_stats.items():
            if stats['sequences'] > 0:
                stats['avg_repeats'] = stats['total_repeats'] / stats['sequences']
                if stats['coverage']:
                    avg_coverage = sum(stats['coverage']) / len(stats['coverage'])
                else:
                    avg_coverage = 0

                f.write(f"  {pattern}:\n")
                f.write(f"    出现在 {stats['sequences']} 个序列中\n")
                f.write(f"    平均重复次数: {stats['avg_repeats']:.2f}\n")
                f.write(f"    最大重复次数: {stats['max_repeats']}\n")
                f.write(f"    平均覆盖率: {avg_coverage:.2f}%\n")

        f.write("\n")

        # 复合模式统计
        composite_counts = {
            "AC-composite": 0,
            "AT-composite": 0,
            "AG-composite": 0
        }

        for results in all_results.values():
            if results['classification'] in composite_counts:
                composite_counts[results['classification']] += 1

        if any(composite_counts.values()):
            f.write("复合模式统计:\n")
            for composite, count in composite_counts.items():
                if count > 0:
                    percentage = (count / len(all_results)) * 100
                    f.write(f"  {composite}: {count} ({percentage:.2f}%)\n")
            f.write("\n")

        # 聚类结果
        if clusters:
            f.write("序列聚类结果:\n")
            f.write(f"  聚类总数: {len(clusters)}\n\n")

            for i, cluster in enumerate(clusters[:5]):  # 只显示前5个聚类
                f.write(f"  聚类 #{i + 1}:\n")
                f.write(f"    大小: {cluster['size']} 序列 ({cluster['percentage']:.2f}%)\n")
                f.write(f"    代表性序列: {cluster['representative']}\n")
                if cluster['size'] <= 5:
                    f.write(f"    成员: {', '.join(cluster['members'])}\n")
                else:
                    f.write(f"    成员: {', '.join(cluster['members'][:5])}... (等 {cluster['size']} 个)\n")
                f.write("\n")

            if len(clusters) > 5:
                f.write(f"  ... 还有 {len(clusters) - 5} 个聚类未显示\n\n")

        # 总体结构特征概述
        f.write("总体结构特征概述:\n")

        # 判断主要分类类型
        main_class = class_counts.most_common(1)[0][0] if class_counts else "Unknown"
        main_class_pct = (class_counts.most_common(1)[0][1] / len(all_results)) * 100 if class_counts else 0

        f.write(f"  主要类型: {main_class} ({main_class_pct:.2f}%)\n")

        # 提供基于数据的总结性描述
        if main_class == "A-rich":
            f.write("  SINE尾部主要呈现富含A特征，多数序列包含高覆盖率的连续A区域。\n")
        elif "AC-composite" in main_class:
            f.write(f"  SINE尾部主要呈现AC系列的复合重复模式特征。\n")
        elif "AT-composite" in main_class:
            f.write(f"  SINE尾部主要呈现AT系列的复合重复模式特征。\n")
        elif "AG-composite" in main_class:
            f.write(f"  SINE尾部主要呈现AG系列的复合重复模式特征。\n")
        elif "AT" in main_class:
            f.write(f"  SINE尾部主要呈现AT重复模式特征。\n")
        elif "AC" in main_class:
            f.write(f"  SINE尾部主要呈现AC重复模式特征。\n")
        elif "AG" in main_class:
            f.write(f"  SINE尾部主要呈现AG重复模式特征。\n")
        else:
            f.write("  SINE尾部呈现多样化的结构特征，没有单一主导模式。\n")


def create_tsv_report(all_results, output_file, family_map=None):
    """
    创建TSV格式的分析结果报告，使用UTF-8编码

    参数:
    all_results (dict): 所有序列的分析结果
    output_file (str): 输出文件
    family_map (dict, optional): 序列ID到家族的映射
    """
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            # 写入表头
            header = [
                "序列ID", "分类", "Family", "使用的截断点", "长度",
                "A含量(%)", "PolyA覆盖率(%)", "最长连续A", "PolyA区域数",
                "重复单元", "重复次数", "重复单元覆盖率(%)",
                "序列复杂度", "序列"
            ]
            f.write("\t".join(header) + "\n")

            # 写入每个序列的数据
            for seq_id, results in all_results.items():
                # 获取重要重复单元
                main_repeat = ""
                if results['repeats']:
                    main_repeat = results['repeats'][0]['unit']

                # 获取家族信息
                family = family_map.get(seq_id, "") if family_map else ""

                # 获取结构分析数据
                classification = results['classification']

                # 提取polyA相关信息
                poly_a_coverage = 0
                poly_a_regions_count = 0
                longest_poly_a = 0
                a_percentage = 0

                # 提取重复模式相关信息
                repeat_count = 0
                pattern_coverage = 0

                # 如果结构分析数据存在
                if 'structure_type' in results:
                    a_percentage = results.get('a_percentage', 0)
                    poly_a_coverage = results.get('poly_a_coverage', 0)
                    poly_a_regions_count = results.get('poly_a_regions_count', 0)
                    longest_poly_a = results.get('longest_poly_a', 0)

                    # 提取重复信息
                    if 'pattern_details' in results:
                        if classification in results.get('pattern_details', {}):
                            pattern_info = results['pattern_details'][classification]
                            repeat_count = pattern_info.get('count', 0)
                            pattern_coverage = pattern_info.get('coverage', 0)
                        elif results.get('repeat_count'):
                            repeat_count = results.get('repeat_count', 0)

                # 准备行数据
                row = [
                    seq_id,
                    classification,
                    family,
                    results.get('cutoff_used', ""),
                    str(results['length']),
                    f"{a_percentage:.2f}",
                    f"{poly_a_coverage:.2f}",
                    str(longest_poly_a),
                    str(poly_a_regions_count),
                    main_repeat,
                    str(repeat_count),
                    f"{pattern_coverage:.2f}",
                    f"{results['complexity']['normalized_entropy']:.3f}",
                    results['sequence']
                ]

                f.write("\t".join(row) + "\n")
    except PermissionError:
        # 使用备用文件名，处理文件权限问题
        import os
        dir_name, file_name = os.path.split(output_file)
        backup_file = os.path.join(dir_name, f"backup_{file_name}")
        print(f"警告: 无法写入文件 {output_file}，尝试使用备用文件名 {backup_file}")

        with open(backup_file, 'w', encoding='utf-8') as f:
            # 写入表头
            header = [
                "序列ID", "分类", "Family", "使用的截断点", "长度",
                "A含量(%)", "PolyA覆盖率(%)", "最长连续A", "PolyA区域数",
                "重复单元", "重复次数", "重复单元覆盖率(%)",
                "序列复杂度", "序列"
            ]
            f.write("\t".join(header) + "\n")

            # 写入每个序列的数据 (与上面相同的逻辑)
            for seq_id, results in all_results.items():
                # 获取重要重复单元
                main_repeat = ""
                if results['repeats']:
                    main_repeat = results['repeats'][0]['unit']

                # 获取家族信息
                family = family_map.get(seq_id, "") if family_map else ""

                # 获取结构分析数据
                classification = results['classification']

                # 提取polyA相关信息
                poly_a_coverage = 0
                poly_a_regions_count = 0
                longest_poly_a = 0
                a_percentage = 0

                # 提取重复模式相关信息
                repeat_count = 0
                pattern_coverage = 0

                # 如果结构分析数据存在
                if 'structure_type' in results:
                    a_percentage = results.get('a_percentage', 0)
                    poly_a_coverage = results.get('poly_a_coverage', 0)
                    poly_a_regions_count = results.get('poly_a_regions_count', 0)
                    longest_poly_a = results.get('longest_poly_a', 0)

                    # 提取重复信息
                    if 'pattern_details' in results:
                        if classification in results.get('pattern_details', {}):
                            pattern_info = results['pattern_details'][classification]
                            repeat_count = pattern_info.get('count', 0)
                            pattern_coverage = pattern_info.get('coverage', 0)
                        elif results.get('repeat_count'):
                            repeat_count = results.get('repeat_count', 0)

                # 准备行数据
                row = [
                    seq_id,
                    classification,
                    family,
                    results.get('cutoff_used', ""),
                    str(results['length']),
                    f"{a_percentage:.2f}",
                    f"{poly_a_coverage:.2f}",
                    str(longest_poly_a),
                    str(poly_a_regions_count),
                    main_repeat,
                    str(repeat_count),
                    f"{pattern_coverage:.2f}",
                    f"{results['complexity']['normalized_entropy']:.3f}",
                    results['sequence']
                ]

                f.write("\t".join(row) + "\n")


def main():
    """
    主函数
    """
    # 设置命令行参数解析器
    parser = argparse.ArgumentParser(
        description="分析猪基因组中SINE转座子尾部的结构模式，根据SINE家族选择截断点策略"
    )
    parser.add_argument(
        "fasta_file",
        help="包含SINE序列的FASTA文件"
    )
    parser.add_argument(
        "--bed_file",
        required=True,
        help="提供SINE家族信息的BED文件(7列格式，第4列为家族，第7列为序列ID)"
    )
    parser.add_argument(
        "--output",
        default="sine_tail_analysis",
        help="输出目录和文件前缀(默认: sine_tail_analysis)"
    )
    parser.add_argument(
        "--cluster",
        action="store_true",
        help="执行序列聚类分析"
    )

    # 解析参数
    args = parser.parse_args()

    # 创建输出目录
    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    details_dir = os.path.join(output_dir, "details")

    # 从FASTA文件读取序列
    print(f"从{args.fasta_file}读取序列...")
    sequences = read_fasta(args.fasta_file)
    print(f"读取了{len(sequences)}个序列。")

    # 读取SINE家族信息
    print(f"从{args.bed_file}读取SINE家族信息...")
    family_map = read_family_info(args.bed_file)
    print(f"读取了{len(family_map)}个SINE家族信息。")

    # 输出一些家族信息样例，帮助调试
    print("家族名称示例:")
    sample_families = list(set(family_map.values()))[:5] if family_map else []
    for i, family in enumerate(sample_families):
        normalized = normalize_family_name(family)
        print(f"  原始名称: {family} -> 标准化名称: {normalized}")

    # 根据家族提取尾部 (使用新逻辑)
    print("根据SINE家族信息提取尾部...")
    tails, motif_used = extract_tails_by_family_new(sequences, family_map)

    print(f"找到{len(tails)}个尾部({len(tails) / len(sequences) * 100:.2f}%的序列)。")

    if not tails:
        print("未找到尾部。请检查您的序列和家族信息。")
        return

    # 分析每个尾部
    print("分析尾部结构模式...")
    all_results = {}

    for seq_id, tail in tails.items():
        results = analyze_tail(tail)
        # 记录使用的截断点
        results['cutoff_used'] = motif_used[seq_id]
        # 添加家族信息
        results['family'] = family_map.get(seq_id, "")
        all_results[seq_id] = results

        # 写入详细报告
        write_detailed_report(seq_id, results, details_dir)

    # 序列聚类(如果启用)
    clusters = []
    if args.cluster:
        print("执行序列聚类...")
        clusters = cluster_sequences(tails)

    # 写入汇总报告
    summary_file = os.path.join(output_dir, "summary_report.txt")
    print(f"生成汇总报告: {summary_file}")
    write_summary_report(all_results, clusters, summary_file)

    # 创建TSV格式报告
    tsv_file = os.path.join(output_dir, "tail_analysis.tsv")
    print(f"生成TSV分析表: {tsv_file}")
    create_tsv_report(all_results, tsv_file, family_map)

    print("\n分析完成!")
    print(f"详细结果保存在 {output_dir} 目录中")


if __name__ == "__main__":
    main()
