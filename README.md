---
title: FishAMR-Link
emoji: 🐟
colorFrom: blue
colorTo: indigo
sdk: streamlit
app_file: FishAMR_Link_v15.py
pinned: false
license: mit
---

# FishAMR-Link (v15.0)

An advanced, single-module computational genomics platform engineered for the multi-dimensional tracking, functional annotation, and machine learning-driven risk stratification of epizootic plasmid vectors conveying Antimicrobial Resistance (AMR) across aquaculture ecosystems.



## Architecture Overview

FishAMR-Link coordinates data processing pipelines across three distinct analytical domains:
1. **Dynamic Functional Annotation:** Orchestrates local multi-threaded subprocess layers executing `prodigal` and `ncbi-blast+` to perform direct structural feature predictions on raw genomic assemblies.
2. **Network Topology Graphing:** Maps spatial and functional associations using proximity analysis to calculate genetic distances and clusters.
3. **Machine Learning Risk Prediction:** Drives an optimized Random Forest classifier utilizing continuous probability arrays to compute precise risk stratification indices.

## Local Deployment Instructions

### Prerequisites
* Miniconda or Anaconda Distribution installed locally.
* Visual Studio Code containing the Python extension pack.

### Environment Setup & Installation
Activate your existing environment inside your terminal and configure the package indexes:

```bash
# Activate your workspace
conda activate bio-app

# Append channels for genomic data handling
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict

# Install command line utilities 
conda install -y prodigal blast-plus

# Install primary runtime application packages
pip install -r requirements.txt
