# Spatiotemporal Modeling & SMC-ABC Calibration Pipeline for Transboundary Crop Pathogens

A data-driven epidemiological modeling framework designed to simulate, calibrate, and forecast the 30-year spatiotemporal invasion trajectory of Cassava Mosaic Disease (CMD/UG variant) across mosaic agroecosystems in West Africa (Côte d'Ivoire, Sierra Leone, and Guinea).This repository integrates high-resolution remote sensing layers (Cassava production density), transportation networks, and field-sampled diagnostics into a gravity-weighted network model. The spatial transmission kernel parameters ($\beta$, $\gamma$) are empirically calibrated using a Sequential Monte Carlo Approximate Bayesian Computation (SMC-ABC) loop running on an adaptive multi-core engine.

## 📂 Repository Architecture

```bash
├── data/                       # Foundational Raw Input Data
│   ├── epidata.xlsx            # Field diagnostics, whitefly counts, & coordinates
│   ├── CassavaMap_Prod_v1.tif  # High-resolution host crop production raster
│   ├── countriesMaps/          # National administrative boundary shapefiles
│   └── roadMaps/               # Regional OpenStreetMap transit shapefiles
├── R/
│   └── library.R               # Unified environment package dependencies
├── scripts/                    # Modulated Execution Pipelines
│   ├── data_cleaning.R         # Coordinate filtering and diagnostic re-coding
│   ├── powerLaw.R              # Vector dispersal kernel mechanics
│   ├── cassMap.R               # Spatial cropping, masking, and terra integration
│   ├── roads.R                 # Transit accessibility and focal density scaling
│   ├── riskMap.R               # Gravity-weighted structural risk surface mapping
│   ├── network_nodes.R         # Scale-aware spatial cell aggregation (10x)
│   ├── ibm_smc_abc.R           # Parallelized Lenormand SMC-ABC Model Calibration
│   └── simulate_progression.R  # 30-Year stochastic projection & high-res animation
├── outputs/                    # Runtime Binaries and Publication Assets
│   ├── nodes.rds               # Aggregated landscape network sites
│   ├── distance_matrix.rds     # Pairwise geodesic Euclidean distances (km)
│   ├── movement_matrix.rds     # Gravity gravity-connectivity weights
│   ├── optimal_abc_parameters.rds # Calibrated posteriors point-estimates
│   ├── UG_HighContrast_RiskMap.png # Map of foundational structural risk
│   └── grounded_epidemic_progression_highres.mp4 # 30-Year video simulation
└── logs/                       # Automated execution tracking and profiling text files

```

## ⚙️ Core Pipeline Execution Sequence

The analytical pipeline must be run sequentially, as each step structurally feeds downstream models via optimized RDS binaries saved to `outputs/`:

### 1. Data Ingestion & Geometry Harmonization
 
* `scripts/data_cleaning.R`: Ingests `data/epidata.xlsx`, strips missing coordinates, re-codes raw diagnostic values into binary disease classifications (e.g., CMD Severity > 1, positive target viruses), converts structures to Simple Features (`sf`), and projects coordinates to *UTM Zone 30N (EPSG:32630)* to lock measurements into raw metric bounds
* `scripts/cassMap.R`: Standardizes multi-country vector shapefiles, resolves boundary self-intersections (`st_make_valid`), and crops/masks the continuous global `CassavaMap_Prod_v1.tif` using highly optimized `terra::rast` workflows.

### 2. Connectivity Network Infrastructure 

* `scripts/roads.R`: Filters regional transit shapefiles down to epidemiologically active vectors (motorways through residential routes), computes focal density over a 5x5 structural grid canvas using a Gaussian blur representation, and extracts normalized road risk layers.
* `scripts/riskMap.R`: Constructs a non-linear multiplicative gravity index balancing localized host presence, market connectivity, and regional long-distance dispersal kernels ($80\text{ km}$ scale factor). Outbreaks print high-resolution visual reference surfaces straight to `outputs/UG_HighContrast_RiskMap.png`.
* `scripts/network_nodes.R`: Mitigates $O(N^2)$ memory limits during distance calculation by aggregating the high-res risk surface by a spatial factor of 10. Computes and serializes the explicit inter-farm Euclidean distance matrix and gravity-based trade interaction models.

### 3. Simulation Calibration & Forecasting

* `scripts/ibm_smc_abc.R`: Hosts an individual-based Susceptible-Infected (SI) mechanistic simulation engine. The model is calibrated against observed spatial field tracking data utilizing an *adaptive Sequential Monte Carlo Approximate Bayesian Computation (SMC-ABC Lenormand)* loop:

$$P(\text{infection}_{i,t}) = 1 - \prod_{j \in I_t} \left(1 - \frac{\beta}{1 + (d_{i,j}/1000)^\gamma}\right)$$

This step leverages multi-core parallel computing to approximate joint posterior distribution vectors for transmission rate ($\beta$) and distance decay ($\gamma$).
* `scripts/simulate_progression.R`: Consumes the optimized point-estimates derived from the ABC step, initialises historical outbreaks from field configurations, and rolls out a forward 30-year stochastic epidemic simulation. The trajectory renders via `gganimate` to compile a high-fidelity video tracking container stored at `outputs/grounded_epidemic_progression_highres.mp4`.

