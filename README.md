# TME-types

**Transcriptomic density calculation from MERFISH slices to a 3D brain atlas**

TME-types is a scientific analysis project for deriving and visualizing spatial transcriptomic density maps from MERFISH imaging data (from the Allen Institute) and aligning them with a three-dimensional brain coordinate framework.

This repository contains scripts and Jupyter notebooks to:

- preprocess MERFISH slice data  
- compute transcriptomic density maps  
- align densities with a 3D brain atlas  
- visualize spatial patterns across brain regions
- validate results against literature

---

## Overview

This it the first part of the pipeline for the paper: *A multimodal spatial atlas of transcriptomic, morphological, and electrophysiological cell type densities in the mouse brain*.

- Notebooks lead you through each step. To improve performance we use multiprocessing, given the current Python version. < soon legacy >
- For larger computational steps we rely on Python scripts. Shell scripts are provided to see resource requirements. 

At the end of this part: we will process MERFISH data from the Allen Institute with density values for all 5k transcriptomic cell types in each brain region.
End products are `.csv` and binary `.pickle` files. 
The next part of the pipeline can be found here: https://github.com/YannRoussel/probabilistic_mapping_extention
For the common coordinate framework (brain atlas) we relied on this work: https://github.com/BlueBrain/ccfv3a-extended-atlas


---

## Repository Structure
```bash
TME-types/
│
├── notebooks/ # Jupyter notebooks for analysis & visualization
├── scripts/ # Python scripts for data processing and density computation
├── README.md # This file
```

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
Paper under review: A multimodal spatial atlas of transcriptomic, morphological, and electrophysiological cell type densities in the mouse brain

## Acknowledgments

The development of this software was supported by funding to the Blue Brain Project, a research center of the École polytechnique fédérale de Lausanne (EPFL), from the Swiss government's ETH Board of the Swiss Federal Institutes of Technology. Copyright (c) 2024 Blue Brain Project/EPFL

Special thanks to all contributors and collaborators who helped!