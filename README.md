# Multiverse-IC

This repository is for my study, where I document importnat information and processes I do not want to forget. It contains Jupyter notebooks and scripts for generating initial conditions (IC) for using CAMB and MUSIC for cosmological simulations, **Multiverse**. The structure is as follows:

## Directory Structure:
```
.
├── notebooks
│   ├── 00.box-cell-setup.ipynb
│   ├── 01.generate_music_input_from_camb.ipynb
│   ├── 02.check-power-sepctrum.ipynb
│   └── 03.about_normalization.ipynb
└── scripts
    ├── camb_transfer_*.txt
    ├── check-power-spectrum.ipynb
    └── generate_music_input_from_camb.ipynb
```

## Contents

1. `notebooks/`: Contains Jupyter notebooks that guide through various steps in setting up and checking initial conditions.
   - `00.box-cell-setup.ipynb`: Initial setup for determining the box size and resolution for simulation.
   - `01.generate_music_input_from_camb.ipynb`: Generates a MUSIC input file based on CAMB transfer functions.
   - `02.check-power-sepctrum.ipynb`: Compares power spectra from CAMB and MUSIC.
   - `03.about_normalization.ipynb`: Discusses the normalization of power spectrum in CAMB and MUSIC.

2. `scripts/`: Contains the condense scripts.
   - `camb_transfer_*.txt`: CAMB transfer function files (aka. MUSIC input files).
   - `check-power-spectrum.ipynb`: Script version for comparing power spectra from CAMB and MUSIC.
   - `generate_music_input_from_camb.ipynb`: Script version for generating a MUSIC input file based on CAMB transfer functions.

## Dependencies

- Python (tested with Python 3.9)
- Jupyter Notebook (tested Jupyte Notebook 6.4.3)
- CAMB (1.5.8)
- MUSIC (1.53)

## Contact

For any questions or issues, please contact Gain Lee (gainlee.cosmos@gmail.com).
