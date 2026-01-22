# Antiproliferative Activity Analysis

This repository contains the full data analysis workflow used to evaluate **cell viability and proliferation (Ki-67)** under different experimental conditions.
The project is organized as a **reproducible, notebook-driven pipeline**, from raw data inspection to exploratory data analysis (EDA) and figure generation for publication.

All analyses were implemented in **Python 3.12** using standard scientific libraries, with an emphasis on transparency, reproducibility, and traceability.

---

## Project Structure

```text
Antiproliferative_activity_analysis/
│
├── environment.yml              # Conda environment for full reproducibility
├── LICENSE                      # GNU GPL v3 (non-commercial use)
├── README.md                    # Project documentation
│
├── raw_data/                    # Original experimental data (immutable)
│   ├── Cell Viability.xlsx
│   ├── KI67.xlsx
│   └── viabilidad vs proliferacion_completo.xlsx
│
├── notebooks/                   # Analysis notebooks (executed sequentially)
│   ├── 1_sanity_checks.ipynb
│   ├── 2_eda_viability.ipynb
│   ├── 3_eda_proliferation.ipynb
│   ├── 4_eda_viability_proliferation.ipynb
│   │
│   └── src/                     # Shared utilities used by notebooks
│       ├── __init__.py
│       └── auxiliar_functions.py
│
├── figures_paper/               # Final, publication-ready figures
│   ├── figure_viability_panels.png/.pdf
│   ├── figure_proliferation_panels.png/.pdf
│   └── figure_relation_viability_proliferation.png/.pdf
│
└── .gitignore
```

---

## Reproducibility Principles

* **Raw data are never modified** (`raw_data/` is read-only by convention)
* All processing steps are explicitly documented in notebooks
* Shared logic is centralized in `notebooks/src/`
* A pinned Conda environment ensures consistent execution
* Figures are generated programmatically and saved to disk

---

## Environment Setup

The recommended way to run this project is via **Conda**.

### 1. Create the environment

```bash
conda env create -f environment.yml
```

### 2. Activate the environment

```bash
conda activate antiproliferative_activity_analysis
```

> The environment includes Python 3.12 and all required scientific libraries
> (NumPy, Pandas, Matplotlib, Seaborn, SciPy, etc.).

---

## ▶How to Run the Analysis

All notebooks are designed to be executed **sequentially**.

### Recommended execution order

| Step | Notebook                              | Purpose                                                   |
| ---: | ------------------------------------- | --------------------------------------------------------- |
|    1 | `1_sanity_checks.ipynb`               | Data integrity checks, column validation, basic summaries |
|    2 | `2_eda_viability.ipynb`               | Exploratory analysis of cell viability data               |
|    3 | `3_eda_proliferation.ipynb`           | Exploratory analysis of Ki-67 proliferation data          |
|    4 | `4_eda_viability_proliferation.ipynb` | Integrated analysis of viability vs proliferation         |

### Launch Jupyter

```bash
jupyter lab
```

or

```bash
jupyter notebook
```

Then open the notebooks from the `notebooks/` directory.

---

## Figures and Outputs

* All **publication-ready figures** are exported automatically to:

  ```
  figures_paper/
  ```

* Figures are saved in both **PNG** (raster) and **PDF** (vector) formats.

* No manual figure editing is required.

---

## Shared Utilities (`notebooks/src/`)

The file `auxiliar_functions.py` contains:

* Data loading helpers
* Standardized preprocessing routines
* Plotting utilities
* Reusable statistical descriptors

This avoids code duplication and ensures consistency across notebooks.

---

## Notes on Statistical Analysis

* When biological or technical replicates are not available (e.g. Ki-67),
  results are treated **descriptively**, without inferential statistics.
* All assumptions and limitations are explicitly stated in the notebooks.
* Statistical tests (when applicable) follow standard best practices
  (normality checks, variance homogeneity, multiple-testing awareness).

---

## License

This project is distributed under the **GNU General Public License v3 (GPL-3.0)**
and is intended for **non-commercial academic use**.

See the `LICENSE` file for details.

---

## Citation & Usage

If you use or adapt this workflow for academic purposes, please cite the associated manuscript or acknowledge the repository.

listo papi — abajo te dejo **la sección de Contact** para **agregar tal cual al README**, escrita en tono académico, clara y estándar para repos científicos. Dejé **campos explícitos** para que completes correos/afiliaciones sin ambigüedades.

---

## 📬 Contact

For questions related to the **experimental data**, **biological interpretation**, or **manuscript content**, please contact:

### **First Author**

* **Name:** *Giannina Maya Hun*
* **Email:** *gianny.maya@gmail.com*

---

### **Corresponding Author**

* **Name:** *Barbara Andrews*
* **Affiliation:** *Departamento de Ingeniería Química, Biotecnología y Materiales, Universidad de Chile.*
* **Email:** *bandrews@ing.uchile.cl*

---

For questions related to the **data analysis pipeline**, **Python scripts**, **notebooks**, or **computational reproducibility**, please contact:

### **Computational Analysis & Code Development**

* **Name:** **David Medina**
* **Affiliation:** *Departamento de Ingeniería en Computación, Universidad de Magallanes*
* **Email:** *david.medina@umag.cl*
