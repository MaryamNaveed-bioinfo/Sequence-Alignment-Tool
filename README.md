# Sequence Alignment Tool

![Python](https://img.shields.io/badge/Python-3.11-blue?logo=python&logoColor=white)
![Tkinter](https://img.shields.io/badge/GUI-Tkinter-orange)
![License](https://img.shields.io/badge/License-MIT-green)

## Table of Contents
1. [Project Overview](#project-overview)
2. [Features](#features)
3. [Folder Structure](#folder-structure)
4. [Installation](#installation)
5. [Usage](#usage)
6. [Algorithms](#algorithms)
7. [Libraries Used](#libraries-used)
8. [License](#license)

## Project Overview
Sequence Alignment Tool is a GUI-based bioinformatics desktop application for performing and visualizing biological sequence alignments. Built with Python and Tkinter, it supports Global Alignment, Local Alignment, and Multiple Sequence Alignment (MSA) with customizable scoring parameters and color-coded visual outputs.

This project demonstrates:
- Python coding skills
- Implementation of bioinformatics algorithms
- Building interactive GUI applications with Tkinter

It is suitable for students and researchers in bioinformatics who need an easy-to-use tool for DNA/protein sequence analysis — no programming knowledge required.

## Features
- Manual sequence input or FASTA file upload
- Global Alignment (Needleman-Wunsch algorithm)
- Local Alignment (Smith-Waterman algorithm)
- Multiple Sequence Alignment (MSA)
- Customizable scoring (match, mismatch, gap penalties)
- Color-coded alignment visualization
- Statistical analysis (identity, similarity, gaps)
- Export results as .txt or .aln files

## Folder Structure
```
Sequence-Alignment-Tool/
├── msa_gui.py        # GUI module — input handling and visualization
├── msa_tool.py       # Backend — alignment algorithms and logic
└── README.md         # This file
```

## Installation

1. Clone the repository:
```bash
git clone https://github.com/MaryamNaveed-bioinfo/Sequence-Alignment-Tool.git
cd Sequence-Alignment-Tool
```

2. Install dependencies:
```bash
pip install biopython
```

## Usage

1. Run the application:
```bash
python msa_gui.py
```

2. Input your sequences manually or upload a FASTA file
3. Set your scoring parameters (match, mismatch, gap penalties)
4. Select alignment type (Global, Local, or MSA)
5. Click **Run Alignment** to see results
6. Export results as needed

## Algorithms
- **Global Alignment** — Needleman-Wunsch algorithm with affine gap penalties
- **Local Alignment** — Smith-Waterman algorithm for optimal local regions
- **Multiple Sequence Alignment** — Progressive alignment using pairwise similarity

## Libraries Used
- [Tkinter](https://docs.python.org/3/library/tkinter.html) — For the graphical user interface
- [Biopython](https://biopython.org/) — For biological sequence handling

## License
This project is open-source and available under the MIT License.
