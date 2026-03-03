# TME-types

**Transcriptomic density calculation from MERFISH slices to a 3D brain atlas**

TME-types is a scientific analysis project for deriving and visualizing spatial transcriptomic density maps from MERFISH imaging data and aligning them with a three-dimensional brain coordinate framework.

This repository contains scripts and Jupyter notebooks to:

- preprocess MERFISH slice data  
- compute transcriptomic density maps  
- align densities with a 3D brain atlas  
- visualize spatial patterns across brain regions  

---

## Overview

Spatial transcriptomics technologies — such as MERFISH — produce high-resolution, molecule-level maps of gene expression across biological tissue sections. TME-types bridges 2D tissue data and a reference 3D brain atlas by calculating *transcriptomic densities* across space and visualizing them in anatomical context.

This enables:

- comparison of spatial patterns between samples  
- integration of multi-slice data into three dimensions  
- deeper insight into structure–expression relationships  

---

## Repository Structure

TME-types/
│
├── notebooks/ # Jupyter notebooks for analysis & visualization
├── scripts/ # Python scripts for data processing and density computation
├── README.md # This file


---

## Getting Started

### Prerequisites

You will need:

- Python (>=3.9 recommended)
- Common scientific packages (NumPy, Pandas, Matplotlib, etc.)
- Jupyter Notebook

---

## Installation

Clone the repository:

```bash
git clone https://github.com/cveraszto/TME-types.git
cd TME-types
```

(Optional) Create a virtual environment:

```bash
python -m venv venv
source venv/bin/activate  # macOS/Linux
venv\Scripts\activate     # Windows
```

Install dependencies (if using a requirements file):

```bash
pip install -r requirements.txt
```

## Usage

Open notebooks in the notebooks/ folder:

These notebooks walk through:

- importing MERFISH slice data
- calculating transcriptomic densities
- visualizing results in 2D and 3D

Run Scripts

The scripts/ folder includes utilities for command-line processing:

```bash
python scripts/your_script.py --input data/ --output results/
```

Replace `your_script.py` with the relevant script.

## Contributing

Contributions are welcome:

Fork the repository

Create a feature branch

Submit a pull request

Please follow standard coding and documentation practices.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

## Citation

If you use this code in published research, please cite the repository appropriately.

## Acknowledgments

Special thanks to all contributors and collaborators who helped bring this project to life!