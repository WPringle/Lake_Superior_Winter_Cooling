[![DOI](https://zenodo.org/badge/265254045.svg)](https://doi.org/10.5281/zenodo.23429888)

# Pringle-etal_2026_GRL

![superior_map](figures/lake_superior_map.png)

**Drivers of Deep-Water Cooling and Model Bias in Lake Superior**

William J. Pringle<sup>1\*</sup>, [COMPASS-GLM Team 1](https://compass.pnnl.gov/GLM/Team)<sup>1,2,3</sup>

<sup>1 </sup>Environmental Science Division, Argonne National Laboratory, Lemont, IL, USA.
<sup>2 </sup>Department of Civil, Environmental and Geospatial Engineering, Michigan Technological University, Houghton, MI, USA.
<sup>3 </sup>Pacific Northwest National Laboratory, Richland, WA, USA.

\* corresponding author:  wpringle@anl.gov

## Abstract
The transition from fall overturn to winter stratification governs ice formation and spring warming timing in deep dimictic lakes. Here, we investigate the drivers of deep-water warm biases reported in prior lake modeling studies during this transition in Lake Superior using subsurface mooring observations and high-quality surface flux measurements. Large-eddy simulations, a k−ε column model, and a one-dimensional analytical model were evaluated over 60-day periods. At the 212-m Eastern Mooring, incorporating thermobaricity and realistic surface flux forcing largely resolves the warm bias, though deep-water mixing differences emerge between models and with wind direction. In contrast, at the 384-m Southern Mooring, warm biases persist across all models; applied surface cooling rates prove insufficient to match observed water-column cooling. We hypothesize that horizontal processes along steep bathymetric troughs, absent in one-dimensional frameworks, drive this discrepancy. Results highlight the sensitivity of winter stratification to cooling rates, wind forcing, thermobaricity and, potentially, three-dimensional dynamics.

## Journal reference
- Pringle, W. J. et al. Drivers of Deep-Water Winter Cooling and Model Bias in Lake Superior. TBD

## Code reference
- Pringle, W. J. (2026). Metarepository for "Drivers of Deep-Water Winter Cooling and Model Bias in Lake Superior" [Metarepository]. Zenodo. https://doi.org/10.5281/zenodo.23429888

## Data reference

### Input data
- Austin, J., & Elmer, C (2023). Lake Superior moored temperature and currents, Sep 2005-May 2015 [Data set]. University Digital Conservancy. https://doi.org/10.13020/zqw9-mk81
- Spence, C., Lenters, J., & Nicholls, E. (2024). US-GL1: Stannard Rock, 2008 - 2022 [Data set]. Ameriflux. https://ameriflux.lbl.gov/sites/siteinfo/US-GL1

### Preprocessing
Input data has been preprocessed and added to `data` directory.
It is also possible to run the following scripts to preprocess the raw subsurface mooring and surface flux [input data](#input-data), e.g., if modification to preprocessing is desired.

| Script Name | Description | How to Run |
| --- | --- | --- |
| `scripts/fill_fluxes.py` | Gap-fills surface flux data | `python3 fill_fluxes.py [input.nc] [output.nc]` |
| `scripts/process_lake_superior_southern_eastern_moorings.jl` | Processes mooring temperature profiles | `julia process_lake_superior_southern_eastern_moorings.jl` |


### Output data
- Pringle, W. J. (2026). Model Outputs for "Drivers of Deep-Water Cooling and Model Bias in Lake Superior"  [Data set]. Zenodo. https://doi.org/10.5281/zenodo.23429567

## Contributing modeling software
| Model | Version | Code | Documentation |
|-------|---------|-----------------|-----|
| Oceananigans.jl | [0.95.3](https://github.com/CliMA/Oceananigans.jl/releases/tag/v0.95.3) | [Github](https://github.com/CliMA/Oceananigans.jl) | [HomePage](https://clima.github.io/OceananigansDocumentation/) |

## Reproduce my experiment
1. Install the software components required to conduct the experiment from [contributing modeling software](#contributing-modeling-software)
2. Optional: Download and preprocess raw [input data](#input-data) required to conduct the experiment (already in `data` directory).
3. Run the Julia `.jl` scripts located in the `scripts` directory (`TURB_experiment` (vertical mixing column model) and `LES_experiment` subdirectories) directly or adapt the PBS job scripts provided to run through all the simulated experiments from the paper. 
4. Download and unzip the [output data](#output-data) from my experiment to compare outputs to those from the publication

## Reproduce my figures
Use the Julia `.jl` scripts found in the `figures` directory to reproduce the figures used in this publication. Filenames of the scripts are numbered to match the figure numbers from the paper.
