# HREA: Hierarchy Ranking Method for Multimodal Multi-objective Problems with Local Pareto Fronts

[![Paper](https://img.shields.io/badge/DOI-10.1109%2FTEVC.2022.3155757-blue)](https://doi.org/10.1109/TEVC.2022.3155757)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![MATLAB](https://img.shields.io/badge/MATLAB-PlatEMO%20v4-orange)

MATLAB code for the paper:

> W. Li, X. Yao, T. Zhang, R. Wang, and L. Wang, "Hierarchy Ranking Method for Multimodal Multi-objective Optimization With Local Pareto Fronts," *IEEE Transactions on Evolutionary Computation*, vol. 27, no. 1, pp. 98–110, 2023.
> DOI: [10.1109/TEVC.2022.3155757](https://doi.org/10.1109/TEVC.2022.3155757)

The paper PDFs are included in this repository (`202302-TEVC-...pdf` and `sm-tevc-hrea.pdf`).

## Overview

Multimodal multi-objective optimization requires finding **all equivalent Pareto-optimal sets**, not just one. HREA introduces a hierarchy ranking method that treats global and local Pareto fronts separately, so that locally non-dominated but globally dominated solutions are preserved as useful diversity instead of being discarded — a key capability on problems with local Pareto fronts (e.g., the IDMP test suite).

## Repository Structure

```
HREA/
├── HREA/                  % the HREA algorithm (drop into PlatEMO)
├── IDMP_e/                % IDMP_e test problems with local Pareto fronts
├── 202302-TEVC-*.pdf      % published paper
└── sm-tevc-hrea.pdf       % supplementary material
```

## Usage

This code is developed based on **PlatEMO v4**:

1. Download [PlatEMO](https://github.com/BIMK/PlatEMO) v4.x.
2. Copy the `HREA` folder into `PlatEMO/Algorithms/Multi-objective optimization/`.
3. Copy the `IDMP_e` folder into `PlatEMO/Problems/`.
4. Run PlatEMO, select algorithm **HREA** and problem **IDMP_e**, and start.

## Citation

If HREA helps your research, please cite:

```bibtex
@article{Li2023HREA,
  author  = {Li, Wenhua and Yao, Xingyi and Zhang, Tao and Wang, Rui and Wang, Ling},
  title   = {Hierarchy Ranking Method for Multimodal Multi-objective Optimization With Local Pareto Fronts},
  journal = {IEEE Transactions on Evolutionary Computation},
  volume  = {27},
  number  = {1},
  pages   = {98--110},
  year    = {2023},
  doi     = {10.1109/TEVC.2022.3155757}
}
```

## Contact

If you have any question, please do not hesitate to contact me: `liwenhua1030@aliyun.com`
