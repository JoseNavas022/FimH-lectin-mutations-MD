# FimH lectin-domain mutations: supplementary data and code

Supporting data, input files and analysis code for the study of three substitutions (R60P, Y137A, Q133K) in the FimH lectin domain (PDB 4XO8) versus wild type, using replicate molecular dynamics (Amber) and stability predictions (FoldX, Rosetta, ERIS-2).

## Contents

| Folder | Description |
|---|---|
| `estructuras/` | Starting crystal structure (`4XO8.pdb`). |
| `construccion_sistemas/` | System preparation: `00_prep/` (PDBFixer, mutagenesis and Pro60 ring-reconstruction scripts) and `mutants/<WT,R60P,Y137A,Q133K>/` (prepared/solvated PDBs, `*.parm7` topologies, `*.rst7` coordinates, `ions.txt`; `build.sh`). |
| `entradas_MD/` | Amber input files: nine-stage restrained equilibration (`s01`-`s09`), production (`prod_first.in`, `prod_next.in`) and replicate-specific inputs (`rep_especificos/`). Random seeds are set via `ig=`. |
| `datos_crudos_dat/` | Raw per-replicate cpptraj time series (`Figure2_cpptraj/`, `Figure3_series_bolsillo/`). |
| `tablas_fuente/` | Source tables for the figures: PCA/FEL (`Figure3_PCA_FEL/`) and DCCM, DSSP and interactions (`Figure5_DCCM_DSSP_interacciones/`). |
| `resultados_crudos_energias/` | Stability predictions: `FoldX/`, `Rosetta/` (relaxed models and scores), `ERIS2_clasico_4XO0/`, and `comparacion_FoldX_Rosetta.tsv`. |
| `scripts/` | Figure-generation scripts (`make_figure2.py`, `make_figure3.py`, `make_figure5.py`, `make_panel_fig1.py`) and figure summary tables. |
| `Figuras_finales/`, `Figura 1/`, `Figure 2/`, `Figure 3/`, `Figure 5/` | Final figures and panels (SVG/PDF/PNG/TIFF) with their generating scripts. |
| `metodos/` | Full methodology, protocol log, traceability notes and completeness check. |
| `Suplementarios_FimH.docx` | Supplementary Material of the article (Figure S1, notes, Tables S1-S3). |

## Notes

- Scripts contain absolute paths from the authors' machine; adjust the root paths before re-running.
- Raw trajectories (`.nc`) are not included because of size.
- tLEaP input and log files used for system building are not included.
