# Estimating Velocities in Spatio-Temporal Point Processes

**Authors**: Fernando Rodríguez Avellaneda, Jorge Mateu and Paula Moraga.  
**Paper**: *Estimating velocities of infectious disease spread through spatio-temporal log-Gaussian Cox point processes* ([arXiv:2409.05036](https://arxiv.org/abs/2409.05036))

This repository provides code supporting the methods and results presented in the paper above, with a focus on estimating velocity fields from spatio-temporal point process data.

---

## Repository Structure

```
estimating-velocities/
├── compute_velocities.R   # Main script for estimating velocity fields
├── create_dataset.R       # Script to simulate or load data
├── modeling.R             # LGCP Spatio temporal modeling using R-INLA
├── params.R               # Parameter definitions for reuse across scripts
├── utilis.R               # Utility functions
├── data/                  # Input data
│   ├── dp.csv             # Population or spatial grid data
│   ├── obs.csv            # Spatio-temporal event data
│   └── lambda1.Rdata      # Precomputed intensity surface
├── images/                # Output plots
│   ├── simulation.pdf
│   ├── intensity.pdf
│   └── velocity.pdf
```

---

## How to Run

1. **Simulate spatio-temporal point pattern**:
   ```r
   source("create_dataset.R")
   ```

2. **Estimate intensity function**:
   ```r
   source("modeling.R")
   ```

3. **Compute velocities using finite differences**:
   ```r
   source("compute_velocities.R")
   ```

---


## Outputs

After running the main scripts, output visualizations can be found in the `images/` folder:

- `simulation.pdf`: simulated data points  
- `intensity.pdf`: fitted intensity surface  
- `velocity.pdf`: estimated velocity 

---

## Citation & BibTeX

If you use this work, please cite:

Rodríguez Avellaneda, F., Mateu, J., & Moraga, P. (2024). Estimating velocities of infectious disease spread through spatio-temporal log-Gaussian Cox point processes. arXiv preprint arXiv:2409.05036. ([https://arxiv.org/abs/2409.05036](https://arxiv.org/abs/2409.05036))

```bibtex
@misc{avellaneda2024estimatingvelocitiesinfectiousdisease,
      title={Estimating velocities of infectious disease spread through spatio-temporal log-Gaussian Cox point processes}, 
      author={Fernando Rodriguez Avellaneda and Jorge Mateu and Paula Moraga},
      year={2024},
      eprint={2409.05036},
      archivePrefix={arXiv},
      primaryClass={stat.AP},
      url={https://arxiv.org/abs/2409.05036}, 
}
```

---

## Contact

**Fernando Rodríguez Avellaneda**  
📧 [fernando.rodriguezavellaneda@kaust.edu.sa](mailto:fernando.rodriguezavellaneda@kaust.edu.sa)

**Jorge Mateu**  
📧 [mateu@uji.es](mailto:mateu@uji.es)

**Paula Moraga**  
📧 [paula.moraga@kaust.edu.sa](mailto:paula.moraga@kaust.edu.sa)
