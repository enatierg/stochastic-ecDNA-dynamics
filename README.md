# Modelling ecDNA Dynamics in Dividing Cell Populations

This project was developed as part of the **Stochastic Modelling Coursework (2024)** during my MMath at the University of Bath.  

It models the dynamics of **extra-chromosomal circular DNA (ecDNA)** in a population of dividing cells. The goal is to understand how ecDNA is inherited, diluted, or lost through stochastic partitioning during cell division, and how this shapes the long-term behaviour of the cell population.

---

## Model Overview

We let **Nk** denote the number of cells that contain *k* copies of ecDNA, where *k* is a non-negative integer. The system evolves according to two coupled differential equations:

1. **Cells with zero ecDNA copies (N₀):**  
   - Increase due to the natural division of cells already lacking ecDNA (assumed to occur at a constant rate of one).  
   - Increase when cells with ecDNA lose all copies during division. This is captured by a probability factor that accounts for all ecDNA copies going to a single daughter cell.  

2. **Cells with k > 0 ecDNA copies (Nk):**  
   - Decrease as these cells divide and are replaced.  
   - Increase due to daughter cells produced during division of other parent cells.  
     - This involves summing over all possible parent cells with sufficiently many ecDNA copies (at least *k/2*), and weighting by both the combinatorial number of ways *k* copies can be partitioned and the probability of that distribution occurring.  

Together, these equations capture the stochastic inheritance of ecDNA, balancing loss and propagation across generations of dividing cells.

---

## Key Features
- **Mathematical Modelling:** Differential equations derived to represent stochastic partitioning of ecDNA.  
- **Biological Relevance:** Provides a framework for understanding how ecDNA dynamics evolve in cancer and other biological contexts.  
- **Coursework Focus:** This project emphasised the use of stochastic models to bridge probability, biology, and applied mathematics.  

---

## Files
- `model_equations.py` – implementation of the differential equations.  
- `simulation.ipynb` – Jupyter notebook with derivations, simulation experiments, and visualisations.  
- `report.pdf` – coursework write-up with full details.  

---

## Requirements
- Python 3.8+  
- NumPy  
- SciPy  
- Matplotlib  

Install all dependencies with:
```bash
pip install -r requirements.txt
