# Example Files for SINE-tail Classification

This directory contains example files to test the SINE-tail classification tool.

## Quick Start

Run the following command from the repository root directory:

```bash
python SINE-tail_classification.py example/input/example_sequences.fasta \
    --bed_file example/input/example_annotations.bed \
    --output example/output/my_results
```

## Input Files

### `input/example_sequences.fasta`

Contains 5 example SINE sequences with different tail types:

| Sequence ID | Expected Tail Type |
|-------------|-------------------|
| SINE_demo_001 | A-rich |
| SINE_demo_002 | (AAAAC)n |
| SINE_demo_003 | (AAAAT)n |
| SINE_demo_004 | (AAAAG)n |
| SINE_demo_005 | (AC)n |

### `input/example_annotations.bed`

BED file (8-column format) containing family annotations for each sequence:

- Column 1: Chromosome
- Column 2: Start position
- Column 3: End position
- Column 4: SINE family
- Column 5: Length
- Column 6: Strand
- Column 7: -0.1
- Column 8: Sequence ID

## Expected Output

After running the tool, you should see the following files in your output directory:

```
my_results/
├── summary_report.txt      # Overall statistics
├── tail_analysis.tsv       # Detailed classification results
└── details/                # Individual sequence reports
    ├── SINE_demo_001.txt
    ├── SINE_demo_002.txt
    ├── SINE_demo_003.txt
    ├── SINE_demo_004.txt
    └── SINE_demo_005.txt
```

### Expected Classification Results

The `tail_analysis.tsv` file should contain classifications similar to:

| Sequence ID | Classification | Family | Tail Length |
|-------------|---------------|--------|-------------|
| SINE_demo_001 | A-rich | SscSINEA1/SINEA | 21 bp |
| SINE_demo_002 | (AAAAC)n | SscSINEA2/SINEA | 23 bp |
| SINE_demo_003 | A-rich | SscSINEB1/SINEB | 18 bp |
| SINE_demo_004 | (AAAAG)n | SscSINEB2/SINEB | 20 bp |
| SINE_demo_005 | (AC)n | SscSINEC1/SINEC | 20 bp |

## Verify Your Installation

If the tool runs successfully and produces output files with the expected classifications, your installation is working correctly.
