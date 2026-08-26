# SPR Kinetics and Affinity Analysis Scripts

Python scripts for processing Surface Plasmon Resonance (SPR) data: kinetics preprocessing/plotting, and steady-state affinity fitting. These scripts were used to generate the SPR figures and fitted binding parameters reported in [Pharmacological inhibition of UBASH3B induces mitotic arrest and selectively reduces tumor growth / DOI — to fill in once available].

The two scripts form a pipeline:

```
raw *_Kinetics.txt  →  spr_kinetics_analysis.py   →  *_SteadyState.txt
                                                            │
                                                            ▼
                                              spr_affinity_fit.py  →  KD, Rmax, Rmin (+ fit plot)
```

## Scripts

### 1. `spr_kinetics_analysis.py`

Reads raw SPR kinetics traces, subtracts a reference signal, normalizes by an analyte/ligand molecular-weight coefficient, extracts the steady-state signal at each concentration, and generates diagnostic/summary plots.

## What it does

For each pair of SPR data files (sample + reference):
1. Reads raw kinetics traces (`*_Kinetics.txt`) and optionally a reference trace for double-referencing.
2. Subtracts the reference signal, scaled by a user-defined coefficient (`CoeffMBP`).
3. Normalizes the subtracted signal by a coefficient computed from ligand immobilization level and analyte/ligand molecular weights (`CoeffNorm_A_L = Immob × MW_analyte / MW_ligand`).
4. Extracts the steady-state signal at each concentration (averaged around the first time point past a defined threshold) and writes it to a `*_SteadyState_NormAL.txt` file.
5. Generates:
   - a multi-panel PDF with every individual sensorgram (data / reference / subtracted) per concentration,
   - a single overlay figure (EPS) with all normalized subtracted curves for the dataset.
     
**Input:** a '*<basename>_Kinetics.txt' file containing Time and Response data from any kind of kinetic data. The file might contain data recorded for several analyte concentrations, and should include the reference.


**Key parameters (top of script)**
| Parameter | Description |
|---|---|
| `DOWNSAMPLING_STEP` | Downsampling factor applied to raw traces |
| `POINTS_TO_REMOVE_AT_END` | Number of trailing points discarded (instrument artifacts) |
| `STEADY_STATE_THRESHOLD` | Time (s) after which the steady-state signal is measured |
| `STEADY_STATE_WINDOW` | Half-width (points) of the averaging window around the steady-state point |
| `LIGANDS` / `ANALYTES` | Molecular weights (kDa) used for normalization |

**Output used downstream:**
- `<file1>_SteadyState_NormAL.txt` — concentration vs. averaged steady-state signal
- `<file1>_Kinetics_Individual_MBP.pdf` — multi-panel diagnostic plots (if `plot_individual_pdf=True`)
- `<file1>_Kinetics_MBPsubtracted.eps` — overlay figure of normalized curves

See inline docstrings in the script for full parameter details.

```python
MAIN_PROCESS(
    file1="<sample_file_basename>",
    file2="<reference_file_basename>",  # "" if no reference
    CoeffMBP=1.0,
    Immob=1000.,
    plot_individual_pdf=True
)
```

### 2. `spr_affinity_fit.py`

Fits a 1:1 binding isotherm (`R = Rmin + (Rmax - Rmin) * C / (C + KD)`) to steady-state SPR data (concentration vs. response), with optional Monte Carlo resampling to estimate parameter uncertainty from replicate/duplicate concentration points.

**Input:** a `*_SteadyState.txt` file (as produced by `spr_kinetics_analysis.py`, or any tab/whitespace-separated file with the same column structure: concentration, response, [optional immobilization level for normalization]).

**Command-line usage:**

```bash
python spr_affinity_fit.py -file=<SteadyStateFile>.txt -mc=1000
```

**Arguments:**

| Argument | Description | Default |
|---|---|---|
| `-file` | Input steady-state file (required) | — |
| `-Rmax` | Fix Rmax in the fit (leave unset to fit it) | not fixed |
| `-Rmin` | Fix Rmin in the fit (leave unset to fit it) | not fixed |
| `-norm` | 1 to normalize by immobilization level (requires a 3rd data column) | 0 |
| `-mc` | Number of Monte Carlo resampling cycles for uncertainty estimation | 0 |

**Output:** fitted parameters printed to console (KD, Rmax, Rmin, RMSD), plus a `.eps` plot of the data and fitted curve.

> **Note:** the script opens an interactive plot window (`plt.show()`) at the end — run it in an environment with display support, or comment out that line for headless/batch execution.

## Requirements

```bash
pip install -r requirements.txt
```

Tested with Python ≥ 3.9.

## Input data format

- `spr_kinetics_analysis.py`: tab-separated `.txt` files exported from the SPR instrument software, named `<basename>_Kinetics.txt` (decimal comma, ISO-8859-1 encoding).
- `spr_affinity_fit.py`: whitespace-separated `.txt` file with one header line, followed by rows of `concentration<TAB>response[<TAB>immobilization_level]`.

> Raw SPR data files used to generate the article's figures are not included in this repository.

## Citation

If you use these scripts, please cite the associated article: [citation — to fill in once the article is published]. See `CITATION.cff` for machine-readable citation metadata.

## License

MIT

## Authors

Yves Nominé — CBI, IGBMC
