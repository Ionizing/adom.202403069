# Reproduction data for the paper "Spin Valley Dynamics Entangled with Optical Fields, Phonons, and Spin-Orbit Coupling in Monolayer MoSe2"

Paper link: https://onlinelibrary.wiley.com/doi/10.1002/adom.202403069

---

## Tool requirements
- NAMD-LMI: https://github.com/Ionizing/NAMD-LMI (0.2.0 or upper)
- rsgrad: https://github.com/Ionizing/rsgrad (0.5.0 or upper)

## Electronic structure at 0 K

- Location: `0K/`

This folder contains the VASP input files to calculate the electronic structure
of monolayer MoSe2 in primitive cell and supercell.

- Reproduction steps:

Just run VASP at each folder containing `INCAR`, remember to copy CHGCAR from corresponding
scf folders when calculating band structures.

## Photo-excitation process at 0 K

- Location: `0K/namd`

This folder contains the NAMD-LMI input files to reproduce the NAMD result at 0 Kelvin.

- Reproduction steps:

1. Go to `0K/2x3/scf` and run `rsgrad model-nac --brange 213 220` to obtain "NAC-0K.h5".
   `rsgrad` can be downloaded from https://github.com/Ionizing/rsgrad/releases ,
   **0.5.5** or higher versions are needed.

2. Go to `0K/namd/single-electron` and run `namd_lmi hamil -c 02_hamil_config.toml` to
   generate `HAMIL-0K-single-electron-efield.h5`.

3. Go to `0K/namd/single-electron/excitation` and run `namd_lmi surfhop -c 03_surfhop_config.toml`
   to run the actual NAMD process.

4. Head to `0K/namd/single-electron/excitation/result` and run `../surfhop_plot.py` to obtain figures.

There are also simulations for systems with multiple electrons at `0K/namd/multi-electron`, and the
reproduction process is similar to the steps shown in the above.

## AIMD at 50 K and 100 K

For 50 K:

- Location: `50K/nve` and `50K/static_ncl`

These two folders contains the NVE trajectory at 50 Kelvin.

- Reproduction steps:

1. Run VASP at `50K/nve` to obtain XDATCAR.

2. Run `genstatic.py` to generate the `run/` folder and move it to `50K/static_ncl`

3. Submit a slurm job using `sub_vasp_namd` at `50K/static_ncl`
   (you may need to modify it to adopt you job scheduling system)

4. Run `mixWave.py` to perform linear unitary transform on degenerated states in WAVECARs.

The process for 100 K is same as 50 K.

## NAMD at 50 K and 100 K

For 50 K:

1. Go to `50K/namd` and run `namd_lmi nac -c 01_nac_config.toml` to calculate non-adiabatic coupling files,
   the couplings will be written to `NAC-50K.h5`

2. Run NAMD at `50K/namd/single-electron` and `50K/namd/multi-electron` with similar processes shown
   in step 2 to step 4 at 0 K.

The process for 100 K is same as 50 K.

You may need the following reference result files, which will be uploaded
to https://github.com/Ionizing/adom.202403069/releases due to GitHub file size limit.

- `50K/namd/NAC-50K.h5`
- `50K/namd/single-electron/HAMIL-50K-single-electron-efield.h5`
- `50K/namd/single-electron/HAMIL-50K-single-electron-noefield.h5`
- `50K/namd/multi-electron/HAMIL-50K-multi-electron-efield.h5`
- `50K/namd/multi-electron/HAMIL-50K-multi-electron-noefield.h5`

- `100K/namd/NAC-100K.h5`
- `100K/namd/single-electron/HAMIL-100K-single-electron-efield.h5`
- `100K/namd/single-electron/HAMIL-100K-single-electron-noefield.h5`
- `100K/namd/multi-electron/HAMIL-100K-multi-electron-efield.h5`
- `100K/namd/multi-electron/HAMIL-100K-multi-electron-noefield.h5`

There are also `hamil_plot.py` to plot figures for `HAMIL-xxx.h5`, use it with `./hamil_plot.py HAMIL-xxx.h5`.
