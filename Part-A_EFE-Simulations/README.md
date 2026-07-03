# Replication Materials: Expected Free Energy as a Structural Architecture of Individual Choice

**Author:** Marc D. Woodfield  
**Paper Version:** Working Paper / Pre-print  
**License:** Creative Commons Attribution–NonCommercial 4.0 International (CC-BY-NC 4.0)

## Overview
This repository contains the reference implementation and full replication materials for the manuscript “Part A - Expected Free Energy as a Structural Architecture of Individual Choice.”

The codebase implements the generative environments, Expected Free Energy (EFE) decision functional, and diagnostic pipelines used to produce all numerical illustrations, identification checks, and discriminability diagnostics reported in:
- Section 4 (Numerical Illustrations)
- Appendix D (Replication Standards)
- Appendix E (Identification and Cross-Fitting Diagnostics W1–W4)

## Scope and Purpose

This repository serves three distinct purposes:

1. **Numerical illustration of structural mechanisms:** Replication of the simulation environments used to illustrate:
   - Identity hysteresis
   - Uncertainty-dependent choice distortion
   - Catastrophic recovery dynamics

2. **Implementation-level identification diagnostics (W1–W4):** Verification that the implemented EFE architecture:
   - Is internally coherent (W1)
   - Supports likelihood-based estimation under correct specification (W2)
   - Exhibits finite-sample recoverability under controlled designs (W3)
   - Is structurally discriminable and falsifiable out-of-sample (W4)

3. **Transparency and auditability:** All outputs are generated from fully specified generative models with fixed seeds, explicit task manipulations, and machine-readable artefacts suitable for hostile-referee inspection.

This repository does not claim empirical validity, behavioural realism, or descriptive adequacy for real-world data. All simulations are stylised and used solely to test internal properties of the decision architecture.

## Reproducibility and Determinism

All scripts enforce explicit random number generation (RNG) control using:

```r
RNGkind(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)
```

and explicit `set.seed()` calls at the script or subject level.

Deterministic reproducibility is intended conditional on:
- Identical random seeds
- Identical RNG settings
- Identical R version and numerical backend

Numerical equality should be assessed using tolerances documented in Appendix D of the manuscript.

## Repository Structure

```text
.
├── src/
│   ├── Part-A_EFE Simulations.R        # Core generative models and EFE decision architecture
│   ├── 01_cr_constant_omega.R          # W1: Constant-Ω structural dynamics
│   ├── 02_fit_cr_params.R              # W2: Single-subject estimation coherence
│   ├── 03_run_W3_recovery.R            # W3: Finite-sample parameter recovery
│   ├── 04_run_W4_crossfit.R            # W4: Cross-fitting discriminability & falsification
│   └── outputs/                        # Auto-generated; initially empty
│       ├── data/
│       ├── figures/
│       └── tables/
└── README.md
```

The `src/outputs/` directory is created automatically at runtime. All reported figures, tables, confusion matrices, and diagnostics are written there verbatim.

## Software Requirements
- R (version 4.0.0 or higher recommended)
- Base R only (no external packages)
This design choice ensures long-term archival stability and removes dependency risk.

## How to Run the Replications

### Option 1: Interactive (RStudio)
1. Open the repository root as an RStudio project.
2. Ensure the working directory is set to the `src/` folder. The scripts assume execution from `src/` and reference paths relative to it.
3. Source the core implementation:
   ```r
   source("Part-A_EFE Simulations.R")
   ```
4. Run individual diagnostics or pipelines as needed, for example:
   ```r
   source("01_cr_constant_omega.R")   # W1
   source("02_fit_cr_params.R")       # W2
   source("03_run_W3_recovery.R")     # W3
   source("04_run_W4_crossfit.R")     # W4
   ```
Each script is self-contained and writes its outputs automatically.

### Option 2: Command Line
From the repository root:
```bash
cd src
Rscript 01_cr_constant_omega.R
Rscript 02_fit_cr_params.R
Rscript 03_run_W3_recovery.R
Rscript 04_run_W4_crossfit.R
```

Scripts may be executed independently or sequentially. Later stages assume outputs only in the conceptual sense, not as file dependencies.

### Output Artefacts
Depending on the script executed, outputs include:
- Deterministic simulation trajectories (`.rds`)
- Parameter recovery tables (`.csv`)
- Confusion matrices and NA audits (`.csv`)
- Subject-level diagnostics (`.csv`)
- Publication-ready figures (`.png`)
- Reproducibility metadata and replication anchors (`_meta.csv`)

Data, tables, diagnostics, and metadata are machine-readable and suitable for independent inspection without re-running simulations. Figures are exported as publication-ready visual artefacts.

### Relationship to the Manuscript
- Section 4: Numerical illustrations are generated directly by the simulation code.
- Appendix D: Reproducibility conventions implemented here.
- Appendix E: W1–W4 diagnostics correspond exactly to scripts 01–04.

Figure numbers and table references in the manuscript correspond to the exported artefacts by filename.

### Citation
If you use or adapt these materials, please cite:

Woodfield, M. D. (2026). *Part A - Expected Free Energy as a Structural Architecture of Individual Choice.* Working Paper. Geneva. DOI: https://doi.org/10.5281/zenodo.18009668

Woodfield, M. D. (2026). *Part A - Expected Free Energy as a Structural Architecture of Individual Choice: Reference implementation and identification diagnostics* [Software]. Available at: https://github.com/mwoodfi/variational-architecture/Part-A_EFE-Simulations

### Disclaimer
This repository contains research code provided for academic and informational purposes only. The software is provided “as is”, without warranty of any kind. Use is entirely at your own risk. For full legal terms, see the `LICENSE` file.

### License
All code and materials are released under the Creative Commons Attribution–NonCommercial 4.0 International (CC-BY-NC 4.0). You are free to use, modify, and redistribute the materials with appropriate attribution, within the limitations of the granted license. Commercial use is not permitted.

### Final Note
This repository is intended to meet elite-journal replication and transparency standards. Any failure to reproduce results under the stated conditions constitutes evidence against the implementation or assumptions, not a hidden tuning choice. No undocumented steps are required.
