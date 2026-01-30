#!/usr/bin/env python3
import re
import os
import argparse
import numpy as np
from collections import Counter, defaultdict
from itertools import product
from Bio import pairwise2
from Bio.Seq import Seq

def read_family_info(bed_file):
    family_map = {}
    if not bed_file or not os.path.exists(bed_file):
        return family_map

    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 7:
                seq_id = fields[6]
                family = fields[3]
                family_map[seq_id] = family

    return family_map

def read_fasta(file_path):
    sequences = {}

    with open(file_path, 'r') as f:
        current_id = ""
        current_seq = ""

        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq and current_id:
                    sequences[current_id] = current_seq

                current_id = line[1:]
                current_seq = ""
            else:
                current_seq += line.upper()

        if current_seq and current_id:
            sequences[current_id] = current_seq

    return sequences

def fuzzy_search(text, pattern, max_mismatches=1):
    pattern_len = len(pattern)
    matches = []

    for i in range(len(text) - pattern_len + 1):
        window = text[i:i + pattern_len]
        mismatches = sum(1 for a, b in zip(window, pattern) if a != b)

        if mismatches <= max_mismatches:
            matches.append((i, i + pattern_len))

    return matches

def normalize_family_name(family):
    if '/' in family:
        left_part = family.split('/')[0]

        if left_part.startswith('Ssc'):
            return left_part[3:]
        return left_part

    return family

def determine_family_group(family):
    if family.startswith("SINEA"):
        return "SINEA"
    elif family.startswith("SINEB"):
        return "SINEB"
    elif family.startswith("SINEC"):
        return "SINEC"
    else:
        return "unknown"

def extract_tails_by_family_new(sequences, family_map):
    tails = {}
    motif_used = {}

    cutoff_strategies = [
        ("AAGAAATAGC", 0, "AAGAAATAGC_exact"),
        ("AAGAAATAGC", 1, "AAGAAATAGC_fuzzy"),
        ("TAGAAAAGGC", 0, "TAGAAAAGGC_exact"),
        ("TAGAAAAGGC", 1, "TAGAAAAGGC_fuzzy"),
        ("TAGAAAAGAC", 0, "TAGAAAAGAC_exact"),
        ("TAGAAAAGAC", 1, "TAGAAAAGAC_fuzzy"),

        ("TAAAAAGAC", 0, "TAAAAAGAC_exact"),
        ("TAAAAAGAC", 1, "TAAAAAGAC_fuzzy"),

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

    family_group_sequences = defaultdict(dict)
    for seq_id, seq in sequences.items():
        original_family = family_map.get(seq_id, "unknown")
        normalized_family = normalize_family_name(original_family)
        family_group = determine_family_group(normalized_family)
        family_group_sequences[family_group][seq_id] = (seq, normalized_family, original_family)

    for family_group, group_seqs in family_group_sequences.items():
        if family_group == "SINEA":
            tail_region_size = 70
        elif family_group in ["SINEB", "SINEC"]:
            tail_region_size = 50
        else:
            tail_region_size = 50

        for seq_id, (seq, normalized_family, original_family) in group_seqs.items():
            seq_len = len(seq)
            if seq_len <= tail_region_size:
                search_region = seq
            else:
                search_region = seq[-tail_region_size:]

            found_match = False
            for pattern, max_mismatches, strategy_name in cutoff_strategies:
                if found_match:
                    break

                matches = fuzzy_search(search_region, pattern, max_mismatches=max_mismatches)
                if matches:
                    match_position = matches[0][0]
                    match_end = matches[0][1]

                    absolute_match_end = (seq_len - tail_region_size) + match_end

                    tail = seq[absolute_match_end:]
                    if tail:
                        tails[seq_id] = tail
                        motif_used[seq_id] = f"{family_group}_{strategy_name}"
                        found_match = True

            if not found_match:
                tail = seq[-30:] if len(seq) > 30 else seq
                tails[seq_id] = tail
                motif_used[seq_id] = f"{family_group}_last30bp"

    return tails, motif_used

def find_poly_a_with_mismatch(sequence, min_a_count=5, max_mismatches=1):
    class PolyAMatch:
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
            
            while j < len(sequence):
                if sequence[j] == 'A':
                    a_count += 1
                    j += 1
                elif mismatch_count < max_mismatches:
                    if j + 1 < len(sequence):
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
    result = {}

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

    a_count = sequence.count('A')
    a_percentage = (a_count / len(sequence)) * 100 if len(sequence) > 0 else 0

    poly_a_regions = find_poly_a_with_mismatch(sequence, min_a_count=5, max_mismatches=1)
    total_poly_a_length = sum(len(match.group()) for match in poly_a_regions)
    poly_a_coverage = (total_poly_a_length / len(sequence)) * 100 if len(sequence) > 0 else 0

    result['a_percentage'] = a_percentage
    result['poly_a_coverage'] = poly_a_coverage
    result['poly_a_regions_count'] = len(poly_a_regions)
    result['longest_poly_a'] = max([len(match.group()) for match in poly_a_regions], default=0) if poly_a_regions else 0

    group_matches = {}
    pattern_details = {}

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

    if not composite_found:
        single_pattern_found = False

        priority_map = {
            "(AAAAC)n": 1, "(AAAC)n": 2, "(AAC)n": 3, "(AC)n": 4,
            "(AAAAT)n": 1, "(AAAT)n": 2, "(AAT)n": 3, "(AT)n": 4,
            "(AAAAG)n": 1, "(AAAG)n": 2, "(AAG)n": 3, "(AG)n": 4
        }

        all_patterns = []
        for group, matches in group_matches.items():
            for match in matches:
                all_patterns.append((match['pattern'], priority_map.get(match['pattern'], 999), match['start']))

        if all_patterns:
            all_patterns.sort(key=lambda x: (x[1], x[2]))

            best_pattern = all_patterns[0][0]
            result['structure_type'] = best_pattern
            result['pattern_details'] = pattern_details
            single_pattern_found = True

        if not single_pattern_found:
            if poly_a_coverage >= 70:
                result['structure_type'] = "A-rich"
            else:
                result['structure_type'] = "Other"

    result['detected_patterns'] = pattern_details

    return result

def new_simplified_classification_advanced(sequence):
    analysis_result = analyze_new_tail_structure_advanced(sequence)

    return analysis_result['structure_type']

def find_repeats(sequence, min_unit_size=1, max_unit_size=5, min_repeats=2):
    repeats = []

    for unit_size in range(min_unit_size, max_unit_size + 1):
        i = 0
        while i <= len(sequence) - unit_size * min_repeats:
            unit = sequence[i:i + unit_size]

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
                i = j
            else:
                i += 1

    repeats.sort(key=lambda x: x['total_length'], reverse=True)

    return repeats

def find_tandem_repeats(sequence, min_period=2, max_period=6):
    repeats = []

    for period in range(min_period, max_period + 1):
        i = 0
        while i <= len(sequence) - 2 * period:
            unit = sequence[i:i + period]
            if sequence[i + period:i + 2 * period] == unit:
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

                i = j
            else:
                i += 1

    repeats.sort(key=lambda x: x['count'] * x['period'], reverse=True)

    return repeats

def analyze_dinucleotide_frequencies(sequence):
    dinucleotides = {}
    for i in range(len(sequence) - 1):
        dinuc = sequence[i:i + 2]
        dinucleotides[dinuc] = dinucleotides.get(dinuc, 0) + 1

    total = sum(dinucleotides.values())
    frequencies = {dinuc: count / total for dinuc, count in dinucleotides.items()} if total > 0 else {}

    sorted_freqs = sorted(frequencies.items(), key=lambda x: x[1], reverse=True)

    return {
        'counts': dinucleotides,
        'frequencies': frequencies,
        'top': sorted_freqs[:5] if sorted_freqs else []
    }

def analyze_sequence_complexity(sequence):
    if not sequence:
        return {'complexity': 0, 'entropy': 0}

    nucleotide_freqs = Counter(sequence)
    total = len(sequence)
    entropy = 0

    for nuc, count in nucleotide_freqs.items():
        p = count / total
        entropy -= p * np.log2(p)

    max_entropy = np.log2(min(4, len(nucleotide_freqs)))
    normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0

    k = 2
    observed_kmers = set()

    for i in range(len(sequence) - k + 1):
        observed_kmers.add(sequence[i:i + k])

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
    results = {
        'sequence': sequence,
        'length': len(sequence)
    }

    structure_analysis = analyze_new_tail_structure_advanced(sequence)

    for key, value in structure_analysis.items():
        results[key] = value

    results['classification'] = structure_analysis['structure_type']

    results['repeats'] = find_repeats(sequence)

    results['tandem_repeats'] = find_tandem_repeats(sequence)

    results['dinucleotide_freq'] = analyze_dinucleotide_frequencies(sequence)

    results['complexity'] = analyze_sequence_complexity(sequence)

    return results

def calculate_sequence_similarity(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)

    if not alignments:
        return 0

    best_alignment = alignments[0]
    alignment_score = best_alignment.score

    max_length = max(len(seq1), len(seq2))
    similarity = alignment_score / max_length if max_length > 0 else 0

    return similarity

def cluster_sequences(sequences, threshold=0.7):
    if not sequences:
        return []

    clusters = []
    seq_ids = list(sequences.keys())

    first_id = seq_ids[0]
    first_cluster = {
        'representative': first_id,
        'members': [first_id],
        'sequence': sequences[first_id]
    }
    clusters.append(first_cluster)

    for seq_id in seq_ids[1:]:
        seq = sequences[seq_id]
        best_similarity = 0
        best_cluster = None

        for i, cluster in enumerate(clusters):
            similarity = calculate_sequence_similarity(seq, cluster['sequence'])
            if similarity > best_similarity:
                best_similarity = similarity
                best_cluster = i

        if best_similarity >= threshold and best_cluster is not None:
            clusters[best_cluster]['members'].append(seq_id)
        else:
            new_cluster = {
                'representative': seq_id,
                'members': [seq_id],
                'sequence': seq
            }
            clusters.append(new_cluster)

    for cluster in clusters:
        cluster['size'] = len(cluster['members'])
        cluster['percentage'] = (cluster['size'] / len(seq_ids)) * 100

    clusters.sort(key=lambda x: x['size'], reverse=True)

    return clusters

def write_detailed_report(seq_id, analysis_results, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    safe_id = seq_id.replace(":", "_").replace(" ", "_").replace("/", "_").replace("\\", "_")
    safe_id = safe_id.replace("*", "_").replace("?", "_").replace("\"", "_").replace("<", "_")
    safe_id = safe_id.replace(">", "_").replace("|", "_")

    with open(os.path.join(output_dir, f"{safe_id}.txt"), 'w') as f:
        f.write(f"SINE尾部结构分析详细报告 - {seq_id}\n")
        f.write("=" * 50 + "\n\n")

        f.write(f"序列: {analysis_results['sequence']}\n")
        f.write(f"长度: {analysis_results['length']} bp\n")
        f.write(f"分类: {analysis_results['classification']}\n")

        if 'family' in analysis_results and analysis_results['family']:
            f.write(f"家族: {analysis_results['family']}\n")

        if 'cutoff_used' in analysis_results:
            f.write(f"使用的截断点: {analysis_results['cutoff_used']}\n")
        f.write("\n")

        repeats = analysis_results['repeats']
        f.write("重复单元分析:\n")
        if repeats:
            for i, repeat in enumerate(repeats[:3]):
                f.write(f"  重复 #{i + 1}: '{repeat['unit']}' (大小: {repeat['size']} bp)\n")
                f.write(f"    出现次数: {repeat['count']}\n")
                f.write(f"    总长度: {repeat['total_length']} bp\n")
        else:
            f.write("  未检测到显著重复单元\n")
        f.write("\n")

        tandem_repeats = analysis_results['tandem_repeats']
        f.write("串联重复分析:\n")
        if tandem_repeats:
            for i, repeat in enumerate(tandem_repeats[:3]):
                f.write(f"  串联重复 #{i + 1}: '{repeat['unit']}' (周期: {repeat['period']} bp)\n")
                f.write(f"    重复次数: {repeat['count']}\n")
        else:
            f.write("  未检测到串联重复\n")
        f.write("\n")

        complexity = analysis_results['complexity']
        f.write("序列复杂度分析:\n")
        f.write(f"  熵: {complexity['entropy']:.4f}\n")
        f.write(f"  归一化熵: {complexity['normalized_entropy']:.4f}\n")
        f.write(f"  复杂度: {complexity['complexity']:.4f}\n")

def write_summary_report(all_results, clusters, output_file):
    with open(output_file, 'w') as f:
        f.write("SINE尾部结构分析汇总报告\n")
        f.write("======================\n\n")

        f.write(f"分析的序列总数: {len(all_results)}\n\n")

        family_counts = Counter(
            [results.get('family', '') for results in all_results.values() if results.get('family', '')])

        if family_counts:
            f.write("SINE家族分布:\n")
            for family, count in family_counts.most_common():
                percentage = (count / len(all_results)) * 100
                f.write(f"  {family}: {count} ({percentage:.2f}%)\n")
            f.write("\n")

        cutoff_counts = Counter(
            [results.get('cutoff_used', '') for results in all_results.values() if results.get('cutoff_used', '')])

        f.write("截断点使用情况:\n")
        for cutoff, count in cutoff_counts.most_common():
            percentage = (count / len(all_results)) * 100
            f.write(f"  {cutoff}: {count} ({percentage:.2f}%)\n")
        f.write("\n")

        classifications = [results['classification'] for results in all_results.values()]
        class_counts = Counter(classifications)

        f.write("分类分布:\n")
        for cls, count in class_counts.most_common():
            percentage = (count / len(all_results)) * 100
            f.write(f"  {cls}: {count} ({percentage:.2f}%)\n")
        f.write("\n")

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

        lengths = [results['length'] for results in all_results.values()]
        avg_length = sum(lengths) / len(lengths) if lengths else 0

        f.write("长度统计:\n")
        f.write(f"  平均长度: {avg_length:.2f} bp\n")
        f.write(f"  中位长度: {sorted(lengths)[len(lengths) // 2] if lengths else 0} bp\n")
        f.write(f"  最小长度: {min(lengths) if lengths else 0} bp\n")
        f.write(f"  最大长度: {max(lengths) if lengths else 0} bp\n\n")

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

        f.write("重复模式统计:\n")

        single_patterns = ["(AAAAC)n", "(AAAC)n", "(AAC)n", "(AC)n",
                           "(AAAAT)n", "(AAAT)n", "(AAT)n", "(AT)n",
                           "(AAAAG)n", "(AAAG)n", "(AAG)n", "(AG)n"]

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

        if clusters:
            f.write("序列聚类结果:\n")
            f.write(f"  聚类总数: {len(clusters)}\n\n")

            for i, cluster in enumerate(clusters[:5]):
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

        f.write("总体结构特征概述:\n")

        main_class = class_counts.most_common(1)[0][0] if class_counts else "Unknown"
        main_class_pct = (class_counts.most_common(1)[0][1] / len(all_results)) * 100 if class_counts else 0

        f.write(f"  主要类型: {main_class} ({main_class_pct:.2f}%)\n")

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
    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            header = [
                "序列ID", "分类", "Family", "使用的截断点", "长度",
                "A含量(%)", "PolyA覆盖率(%)", "最长连续A", "PolyA区域数",
                "重复单元", "重复次数", "重复单元覆盖率(%)",
                "序列复杂度", "序列"
            ]
            f.write("\t".join(header) + "\n")

            for seq_id, results in all_results.items():
                main_repeat = ""
                if results['repeats']:
                    main_repeat = results['repeats'][0]['unit']

                family = family_map.get(seq_id, "") if family_map else ""

                classification = results['classification']

                poly_a_coverage = 0
                poly_a_regions_count = 0
                longest_poly_a = 0
                a_percentage = 0

                repeat_count = 0
                pattern_coverage = 0

                if 'structure_type' in results:
                    a_percentage = results.get('a_percentage', 0)
                    poly_a_coverage = results.get('poly_a_coverage', 0)
                    poly_a_regions_count = results.get('poly_a_regions_count', 0)
                    longest_poly_a = results.get('longest_poly_a', 0)

                    if 'pattern_details' in results:
                        if classification in results.get('pattern_details', {}):
                            pattern_info = results['pattern_details'][classification]
                            repeat_count = pattern_info.get('count', 0)
                            pattern_coverage = pattern_info.get('coverage', 0)
                        elif results.get('repeat_count'):
                            repeat_count = results.get('repeat_count', 0)

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
        import os
        dir_name, file_name = os.path.split(output_file)
        backup_file = os.path.join(dir_name, f"backup_{file_name}")
        print(f"警告: 无法写入文件 {output_file}，尝试使用备用文件名 {backup_file}")

        with open(backup_file, 'w', encoding='utf-8') as f:
            header = [
                "序列ID", "分类", "Family", "使用的截断点", "长度",
                "A含量(%)", "PolyA覆盖率(%)", "最长连续A", "PolyA区域数",
                "重复单元", "重复次数", "重复单元覆盖率(%)",
                "序列复杂度", "序列"
            ]
            f.write("\t".join(header) + "\n")

            for seq_id, results in all_results.items():
                main_repeat = ""
                if results['repeats']:
                    main_repeat = results['repeats'][0]['unit']

                family = family_map.get(seq_id, "") if family_map else ""

                classification = results['classification']

                poly_a_coverage = 0
                poly_a_regions_count = 0
                longest_poly_a = 0
                a_percentage = 0

                repeat_count = 0
                pattern_coverage = 0

                if 'structure_type' in results:
                    a_percentage = results.get('a_percentage', 0)
                    poly_a_coverage = results.get('poly_a_coverage', 0)
                    poly_a_regions_count = results.get('poly_a_regions_count', 0)
                    longest_poly_a = results.get('longest_poly_a', 0)

                    if 'pattern_details' in results:
                        if classification in results.get('pattern_details', {}):
                            pattern_info = results['pattern_details'][classification]
                            repeat_count = pattern_info.get('count', 0)
                            pattern_coverage = pattern_info.get('coverage', 0)
                        elif results.get('repeat_count'):
                            repeat_count = results.get('repeat_count', 0)

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

    args = parser.parse_args()

    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)
    details_dir = os.path.join(output_dir, "details")

    print(f"从{args.fasta_file}读取序列...")
    sequences = read_fasta(args.fasta_file)
    print(f"读取了{len(sequences)}个序列。")

    print(f"从{args.bed_file}读取SINE家族信息...")
    family_map = read_family_info(args.bed_file)
    print(f"读取了{len(family_map)}个SINE家族信息。")

    print("家族名称示例:")
    sample_families = list(set(family_map.values()))[:5] if family_map else []
    for i, family in enumerate(sample_families):
        normalized = normalize_family_name(family)
        print(f"  原始名称: {family} -> 标准化名称: {normalized}")

    print("根据SINE家族信息提取尾部...")
    tails, motif_used = extract_tails_by_family_new(sequences, family_map)

    print(f"找到{len(tails)}个尾部({len(tails) / len(sequences) * 100:.2f}%的序列)。")

    if not tails:
        print("未找到尾部。请检查您的序列和家族信息。")
        return

    print("分析尾部结构模式...")
    all_results = {}

    for seq_id, tail in tails.items():
        results = analyze_tail(tail)
        results['cutoff_used'] = motif_used[seq_id]
        results['family'] = family_map.get(seq_id, "")
        all_results[seq_id] = results

        write_detailed_report(seq_id, results, details_dir)

    clusters = []
    if args.cluster:
        print("执行序列聚类...")
        clusters = cluster_sequences(tails)

    summary_file = os.path.join(output_dir, "summary_report.txt")
    print(f"生成汇总报告: {summary_file}")
    write_summary_report(all_results, clusters, summary_file)

    tsv_file = os.path.join(output_dir, "tail_analysis.tsv")
    print(f"生成TSV分析表: {tsv_file}")
    create_tsv_report(all_results, tsv_file, family_map)

    print("\n分析完成!")
    print(f"详细结果保存在 {output_dir} 目录中")

if __name__ == "__main__":
    main()
