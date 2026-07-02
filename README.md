---
title: FishAMR-Link v17.0
emoji: 🧬
colorFrom: green
colorTo: indigo
sdk: streamlit
sdk_version: 1.35.0
app_file: FishAMR-Link_v17.0.py
pinned: false
license: mit
---

# FishAMR-Link v17.0

**Automated Inference of Plasmid-Mediated Antimicrobial Resistance Transmission in Aquatic Pathogens**

FishAMR-Link v17.0 is an enterprise-grade, multidisciplinary structural bioinformatics pipeline designed to detect, track, and model plasmid-mediated antimicrobial resistance (AMR), mobile genetic elements (MGEs), and virulence factors across epizootic aquaculture reservoirs. 

The application integrates deterministic biophysical sequence modeling with an explainable ensemble machine learning classifier, backed by an advanced Heuristic Reservoir Similarity Scoring (HRSS) topology map to infer horizontal gene transfer (HGT) pathways among aquatic pathogens.

---

## ✨ Key Architectural Features

### 1. Multi-Tiered Structural Genomics Engine
* **Native Production Track:** Dynamically intercepts the runtime environment using `shutil.which` to execute local installations of NCBI's `AMRFinderPlus` for core AMR alignment against curated databases.
* **Deterministic Fallback Array:** If the environment lacks local binaries (e.g., lightweight Hugging Face Space basic instances), an internal translation engine activates. It processes six-frame reading frames under standard **Bacterial Translation Table 11** mechanics and applies strict peptide motif regex parsing to isolate key targets (e.g., `tetA`, `sul1`, `floR`, `blaCTX-M`, `mcr-1`).

### 2. K-Mer Sketching & Graph Transmission Topologies
* Utilizes **21-mer bottom-400 MinHash sketching** with cryptographic `Blake2b` hashes to construct cross-isolate containment coefficients.
* Combines **Jaccard containment distances** (40%), **GC profile variances** (40%), and **spatial/temporal metadata tracking vectors** (20%) into a unified **Heuristic Reservoir Similarity Score (HRSS)**, visualizing transmission linkages through interactive directional `NetworkX` and `Plotly` graphs.

### 3. Explainable Machine Learning Engine
* Employs a custom `Random Forest` architecture trained over an internal relational SQLite data warehouse containing historic epizootic reference baselines.
* Provides rigorous, honest validation diagnostics via **Leave-One-Out Cross-Validation (LOOCV)** to guard against model overfitting on small synthetic sample matrices and logs real-time feature weight contributions.

### 4. Biosecurity Provenance Auditing
* Implements a secure `ReportLab` design block to export comprehensive, production-grade PDF compliance records detailing geographic localization vectors, compositional skews, identified elements, and calculated threat probability scores.

---

## 🚀 Quick Start & Deployment

### Run Locally (Standard Python Environment)
1. Clone the repository and install the dependencies:
   ```bash
   git clone [https://github.com/your-username/FishAMR-Link.git](https://github.com/your-username/FishAMR-Link.git)
   cd FishAMR-Link
   pip install -r requirements.txt
