# salmon-ibm-model

The SalmonIRA (Salmon Inflation Reduction Act) project simulates the migration of Chinook salmon smolts along the U.S. West Coast using Lagrangian particle tracking methods with particles released into MOM6 oceanographic flow fields.

---

##  Project Structure
This repo contains only the model and files required to run simulations locally or on HPC. For analysis notebooks and figures, see the companion repository: SalmonIRA

---

## Requirements

- Python 3.9+
- [Parcels](https://oceanparcels.org/) (particle tracking)
- xarray
- numpy
- matplotlib

---

##  How to Run

Navigate to the project directory and run the model script with the required parameters:

```bash
python MOM6salmonira.py <seed> <nfish> <kf> <kb> <ki> <pavg>
```

### Example Run

```bash
python MOM6salmonira.py 321 10 0 1 0.01 0 0.1
```

*Make sure your Conda environment is activated and all required flow field data files (surface u velocity, surface v velocity, (later: surface temperature), geometry) are available.*

###  Parameter Description

| Parameter | Description | Example Value |
|:---------:|:------------|:-------------:|
| `<seed>` | random seed | `321` |
| `<nfish>` | number of fish (particle-type) | `10` |
| `<kf>` | kappaF: relative contribution of other fish to fish swimming direction | `0` |
| `<kb>` | kappaB: relative contribution of bathymetry to fish swimming direction | `0.5` |
| `<ki>` | kappaI: relative contribution of inertia to fish swimming direction (constant) | `0.01` |
| `<pavg>` | mean prey index per grid cell | `0.1` |

Original model (OP_tuna/tunaFADpreyF.py & OP_tuna/scratchtunaFADpreykernels.py) from https://doi.org/10.1016/j.ecolmodel.2022.110188

---

##  Outputs

Running the model will generate a `.zarr` output file containing particle trajectories (later: and temperature data).

---

##  Acknowledgments

This work is supported by:
- NOAA Northwest Fisheries Science Center
- Oregon State University's Lagrangian Ocean Laboratory
