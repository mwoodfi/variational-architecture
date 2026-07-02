# Research Programme Materials: Variational Foundations of Agency and Social Order

**Author:** Marc D. Woodfield  
**Programme Version:** Working Paper / Pre-print  
**License:** Creative Commons Attribution–NonCommercial 4.0 International (CC-BY-NC 4.0)

## Overview
This repository contains code, computational materials, figures, diagnostic artefacts, implementation notes, and related technical materials for the research programme “Variational Foundations of Agency and Social Order: Decision Architectures, Physics of Choice, and Constitutional Agency.”

The research programme develops a multi-level variational architecture linking:
  - Part A: Individual choice and policy selection
  - Part B: Constitutional agency, identity, preferences, and precision
  - Part C: Intragroup coordination and collective order
  - Part D: Intergroup dynamics, strategic interaction, and institutional stabilisation

Citable manuscripts and formal working-paper versions are archived on Zenodo under their respective DOIs. This GitHub repository is maintained as the computational and technical companion to the research programme.

## Scope and Purpose
This repository serves three distinct purposes:

1.	Programme-level transparency and organisation: Public coordination of the computational and technical materials associated with the research programme, including:
    - Code and simulation scripts
    - Replication materials
    - Figures and diagrams
    - Diagnostic artefacts
    - Implementation notes
    - Supplementary technical materials

2.	Component-level implementation and replication: Where available, provision of computational materials linked to the individual Parts of the programme:
    - Part A: Individual choice and Expected Free Energy policy selection
    - Part B: Constitutional agency and slow parameter formation
    - Part C: Intragroup dynamics and collective inference
    - Part D: Intergroup dynamics and strategic interaction

3.	Transparency and auditability: Repository materials are provided to support inspection, replication, extension, and independent re-analysis of the computational components of the programme.

This repository does not claim empirical validity, behavioural realism, or descriptive adequacy for real-world data. Computational materials should be interpreted only in relation to the assumptions, scope conditions, and interpretation limits stated in the corresponding manuscripts.

## Reproducibility and Determinism
Where simulation scripts or diagnostic pipelines are provided, reproducibility conditions are documented in the relevant component-level README files or replication notes.

Unless otherwise stated, deterministic reproducibility is conditional on:
- Identical random seeds
- Identical random number generation settings
- Identical software versions
- Identical numerical backend
- The repository state, commit, release, or archive associated with the cited manuscript version

Numerical equality should be assessed using the tolerances documented in the corresponding manuscript or component-level replication materials.

Because the live GitHub repository may change over time, exact replication should rely on the code state associated with the cited manuscript version.

## Repository Structure
The repository is organised as a programme-level coordination repository.

The exact structure may evolve as the research programme develops. Component-specific folders and README files should be consulted for exact replication instructions, software requirements, file structure, and execution order.

Citable manuscripts and formal working-paper versions are archived on Zenodo. GitHub is used for code, computational materials, figures, diagnostic artefacts, implementation notes, and related technical materials.

## Software Requirements
Software requirements are component-specific.

Where code is provided, the relevant component-level README or replication note specifies:
  - Required programming language
  - Required software version
  - Package dependencies, if any
  - Execution instructions
  - Output locations
  - Reproducibility conventions

No programme-wide software environment is assumed unless explicitly stated.

## How to Use This Repository

### Step 1: Identify the relevant component
Determine which part of the research programme the material belongs to:

  - Part A: Individual choice
  - Part B: Constitutional agency
  - Part C: Intragroup dynamics
  - Part D: Intergroup dynamics

### Step 2: Read the corresponding manuscript
Repository materials should be interpreted together with the relevant manuscript assumptions, scope limits, and formal definitions.

The manuscripts remain the authoritative source for:
  - Formal claims
  - Assumption sets
  - Theorem statements
  - Proof structure
  - Interpretation limits
  - Citation requirements

### Step 3: Consult the component-level README
Where computational materials are available, the component-level README or replication note provides the relevant execution instructions.

### Step 4: Run, inspect, or extend the materials
Code and artefacts are provided to support:
  - Identification of errors, boundary cases, or implementation failures
  - Audit
  - Replication
  - Independent re-analysis
  - Extension under modified assumptions

### Output Artefacts
Depending on the component and script executed, outputs may include:
	- Simulation trajectories (.rds or equivalent)
	- Parameter recovery tables (.csv)
	- Confusion matrices and diagnostic audits (.csv)
	- Subject-level or run-level diagnostics (.csv)
	- Figures and diagrams (.png, .pdf, or equivalent)
	- Reproducibility metadata
	- Implementation notes

All artefacts should be read in relation to the corresponding manuscript. Simulation outputs do not, by themselves, constitute external empirical validation.

### Relationship to the Manuscripts
	* Research Programme document: Defines the overarching architecture linking Parts A–D.
	* Part A: Individual choice and micro-level policy selection.
	* Part B: Constitutional agency and slow formation of identity, preferences, and precision.
	* Part C: Intragroup dynamics and collective inference.
	* Part D: Intergroup dynamics and strategic interaction.

Manuscripts are archived on Zenodo. Repository materials provide computational and technical support where available.

### Citation
If you use or adapt these materials, please cite the relevant manuscript and, where applicable, the repository or archived code state.

Research programme:
Woodfield, M. D. (2026). Research Program - Variational Foundations of Agency and Social Order: Decision Architectures, Physics of Choice, and Constitutional Agency. Working Paper. Geneva. DOI: https://doi.org/10.5281/zenodo.18513521

Repository:
Woodfield, M. D. (2026). Variational Foundations of Agency and Social Order: Code, computational materials, figures, diagnostic artefacts, and implementation notes [Repository]. Available at: https://github.com/mwoodfi/variational-architecture

Component papers and component-specific software materials should be cited separately where applicable.

### Disclaimer
This repository contains research materials provided for academic and informational purposes only. The materials are provisional, may contain errors or omissions, and are subject to revision. Code is provided “as is”, without warranty of any kind. Use is entirely at your own risk. For full legal terms, see the LICENSE file.

Nothing in this repository constitutes professional advice, policy guidance, empirical validation, or normative recommendation.

### License
All code and materials are released under the Creative Commons Attribution–NonCommercial 4.0 International (CC-BY-NC 4.0), unless otherwise stated. You are free to use, modify, and redistribute the materials with appropriate attribution, within the limitations of the granted license. Commercial use is not permitted.

### Final Note
This repository is intended to support transparency, auditability, replication, and critical scrutiny of the research programme. Any failure to reproduce computational results under the stated conditions constitutes evidence against the relevant implementation, numerical environment, or maintained assumptions, not a hidden tuning choice. No undocumented steps are intended to be required.
