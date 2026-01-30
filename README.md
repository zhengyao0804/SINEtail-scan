# SINE-tail Classification Tool

A Python tool for analyzing and classifying tail structures of SINE (Short Interspersed Nuclear Elements) transposons in pig genomes.

## Description

This tool analyzes the complex structural patterns of SINE transposon tails in pig genomes. It applies family-specific truncation point matching strategies based on SINE family information:

- **SINEA1-11 families**
- **SINEB1-6 families**
- **SINEC1-8 families**

### Tail Structure Classification

The tool classifies tail structures into the following categories:

| Category                         | Description                                 |
| -------------------------------- | ------------------------------------------- |
| A-rich                           | Consecutive A content ≥70% of tail sequence |
| (AAAAC)n, (AAAC)n, (AAC)n, (AC)n | AC-format repeats                           |
| AC-composite                     | Composite patterns of AC-series             |
| (AAAAT)n, (AAAT)n, (AAT)n, (AT)n | AT-format repeats                           |
| AT-composite                     | Composite patterns of AT-series             |
| (AAAAG)n, (AAAG)n, (AAG)n, (AG)n | AG-format repeats                           |
| AG-composite                     | Composite patterns of AG-series             |
| Other                            | Other structural patterns                   |

## Requirements

### Software

- **Python** ≥ 3.6

### Python Dependencies

```
numpy
biopython
```

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/zhengyao0804/SINEtail-scan
cd SINE-tail-classification
```

### 2. Install dependencies

Using pip:

```bash
pip install numpy biopython
```

Or using conda:

```bash
conda install numpy biopython
```

Or install from requirements file:

```bash
pip install -r requirements.txt
```

### 3. Verify installation

```bash
python SINE-tail_classification.py --help
```

## Usage

### Basic Command

```bash
python SINE-tail_classification.py <fasta_file> --bed_file <bed_file> [options]
```

### Arguments

| Argument     | Required | Description                                                  |
| ------------ | -------- | ------------------------------------------------------------ |
| `fasta_file` | Yes      | FASTA file containing SINE sequences                         |
| `--bed_file` | Yes      | BED file providing SINE family information (7-column format: column 4 = family, column 7 = sequence ID) |
| `--output`   | No       | Output directory and file prefix (default: `sine_tail_analysis`) |

### Examples

**Basic analysis:**

```shell
python SINE-tail_classification.py sequences.fasta --bed_file annotations.bed
```

**Specify output directory:**

```shell
python SINE-tail_classification.py sequences.fasta --bed_file annotations.bed --output my_results
```

## Input File Formats

### FASTA File

Standard FASTA format containing SINE sequences:

```
>SscSINEA1#SINE/SINEA
GGAGTTCCCGTCGTGGCGCAGTGGTTAACGAATCCGACTAGGAACCATGAGGTTGCGGGTTCGGTCCCTGCCCTTGCTCAGTGGGTTAACGATCCGGCGTTGCCGTGAGCTGTGGTGTAGGTTGCAGACGCGGCTCGGATCCCGCGTTGCTGTGGCTCTGGCGTAGGCCGGTGGCTACAGCTCCGATTCGACCCCTAGCCTGGGAACCTCCATATGCCGCGGGAGCGGCCCAAGAAATAGCAACAACAACAACAACAAAAAAGACAAAAAGACCAAAAAAAAAAAAAAAAAAAA
>SscSINEA2#SINE/SINEA
GGAGTTCCCGTCGTGGCGCAGTGGTTAACGAATCCGACTAGGAACCATGAGGTTGCGGGTTCGGTCCCTGCCCTTGCTCAGTGGGTTAACGATCCGGCGTTGCCGTGAGCTGTGGTGTAGGTTGCAGACGCGGCTCGGATCCCGCGTTGCTGTGGCTCTGGCGTAGGCCGGCGGCTACAGCTCCGATTCGACCCCTAGCCTGGGAACCTCCATATGCCGCGGGAGCGGCCCAAGAAATAGCAAAAAGACAAAAAAAAAAAAAA
```

### BED File (8-column format)

Tab-separated file with the following columns:

| Column | Description                                 |
| ------ | ------------------------------------------- |
| 1      | Chromosome                                  |
| 2      | Start position                              |
| 3      | End position                                |
| 4      | SINE family (e.g., SscSINEA1/SINEA, SINEB2) |
| 5      | Length                                      |
| 6      | Strand                                      |
| 7      | -1.0                                        |
| 8      | Serial number                               |

Example:

```
chr1    469     714     SscSINEA9       245     -       SINE    SINEA   -1.0    1
chr1    1076    1323    SscSINEA7       247     -       SINE    SINEA   -1.0    2
chr1    6359    6581    SscSINEB4       222     +       SINE    SINEB   -1.0    3
chr1    6942    7205    SscSINEA4       263     +       SINE    SINEA   -1.0    4
```

## Output Files

The tool generates the following output files in the specified output directory:

```
output_directory/
├── summary_report.txt      # Summary statistics and analysis
├── tail_analysis.tsv       # Detailed results in TSV format
└── details/                # Individual sequence analysis reports
    ├── sequence_id_1.txt
    ├── sequence_id_2.txt
    └── ...
```

### Output File Descriptions

#### `summary_report.txt`

Contains overall statistics including:

- Total number of sequences analyzed
- Distribution of tail structure classifications
- Family-specific statistics

#### `tail_analysis.tsv`

Tab-separated file with columns:

| Column            | Description                              |
| ----------------- | ---------------------------------------- |
| 序列ID            | Sequence identifier                      |
| 分类              | Tail structure classification            |
| Family            | SINE family                              |
| 使用的截断点      | Truncation point used                    |
| 长度              | Tail length                              |
| A含量(%)          | Adenine content percentage               |
| PolyA覆盖率(%)    | Poly-A coverage percentage               |
| 最长连续A         | Longest consecutive A stretch            |
| PolyA区域数       | Number of poly-A regions                 |
| 重复单元          | Main repeat unit                         |
| 重复次数          | Repeat count                             |
| 重复单元覆盖率(%) | Repeat unit coverage                     |
| 序列复杂度        | Sequence complexity (normalized entropy) |
| 序列              | Tail sequence                            |

## Example

An `example/` directory is provided to test the tool:

```bash
# Run example analysis
python SINE-tail_classification.py example/input/example_sequences.fasta \
    --bed_file example/input/example_annotations.bed \
    --output example/output/test_results
```

## Citation

If you use this tool in your research, please cite:

> Not yet

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

- **Author**: Yao Zheng
- **Email**: MZ120180996@yzu.edu.cn

## Version History

- **v1.0.0** (2024-XX-XX): Initial release accompanying *Genes* manuscript submission
