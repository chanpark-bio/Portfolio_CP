# CP In Silico Hub: Automated De Novo Protein Design and Developability Pipeline

![Python](https://img.shields.io/badge/Python-3.9+-blue.svg)
![Bash](https://img.shields.io/badge/Bash-Scripting-black.svg)
![AI](https://img.shields.io/badge/AI_Engine-RFdiffusion%20%7C%20ProteinMPNN-orange)
![Physics](https://img.shields.io/badge/Physics_Engine-FoldX-green)
![QC](https://img.shields.io/badge/Cascade_QC-NetMHCpan%20%7C%20A3D-red)

## 1. Overview
CP In Silico Hub is an automated, end-to-end computational pipeline for *de novo* protein binder design. Integrating generative AI models with rigorous physics-based and multi-parametric developability filters (Cascade QC), this architecture evaluates sequence stability, aggregation propensity, and immunogenicity prior to wet-lab synthesis. This repository demonstrates the pipeline's workflow using the SARS-CoV-2 Spike RBD (PDB: 6M0J) as a proof-of-concept target.

## 2. Architecture & Workflow
The workflow is orchestrated via a centralized `config.yaml` and bash scripts, enabling execution across six modular phases:

* **Target Profiling:** Automated parsing and structural pre-processing of the input target PDB.
* **Scaffold Generation:** Employs RFdiffusion to generate target-aware backbones based on user-defined contigs, followed by an initial triage module to filter scaffolds via pseudo-RMSD and pLDDT metrics.
* **Sequence Design:** Executes high-throughput sequence generation using ProteinMPNN on the filtered elite scaffolds.
* **Developability Filtering (Fast Cascade QC):** Discards highly unstable sequences, identifies hydrophobic surface patches predicting DSP solubility issues via Aggrescan3D, and applies a strict 0.0 epitope cutoff using NetMHCpan to minimize Anti-Drug Antibody (ADA) risks.
* **Physics Engine:** Calculates the complex binding energy ($\Delta \Delta G$) using FoldX to estimate binding affinity.
* **Dossier Generation:** Aggregates all structural and physicochemical metrics into a consolidated CSV/HTML report for candidate selection.

## 3. Technical Highlights
* **Integration of DSP Criteria:** Binders with high affinity often fail during downstream processing or clinical stages due to aggregation or immunogenicity. This pipeline mitigates these risks by embedding developability checks directly into the early discovery loop.
* **Reproducibility & Resilience:** A SQLite-based tracking system logs the status of each binder at every computational phase. This architecture allows the system to safely resume interrupted tasks or overwrite specific modules without redundant processing.

## 4. Usage (Demo)
The current repository is configured for a scaled-down run targeting `6M0J`.

```bash
# Initialize project workspace and define contig (e.g., E333-526/0 15-50)
./init.sh

# Execute the pipeline
./run.sh

Output data, including scatter plots and CSV reports, will be generated in the 03_Workspace/Portfolio_CoV2/04_Fast_QC/02_qc_reports directory.
5. Future Directions: Modality Engineering

Subsequent phases currently under development aim to bridge computational discovery with practical bioprocessing operations:

    DSP Condition Prediction: Modeling optimal AC, TFF, IEX, and SEC parameters based on surface charge (pI) and hydrophobicity profiles.

    Expression Optimization: Codon optimization tailored for CHO cell systems.

    PK/PD Considerations: Exploring Fc/HSA fusion modalities for in vivo half-life extension.

Note: This is a sanitized version of the pipeline intended for portfolio demonstration. Proprietary institutional targets, internal data, and specialized hyper-parameters have been removed or set to standard public defaults.
