[![DOI](https://zenodo.org/badge/265254045.svg)](https://zenodo.org/doi/10.5281/zenodo.10442485)

# Pringle-etal_2026_GRL

**Drivers of Deep-Water Cooling and Model Bias in Lake Superior**

William J. Pringle<sup>1\*</sup>, [GLEAM Team 1](https://compass.pnnl.gov/GLM/Team)<sup>1,2,3</sup>

<sup>1 </sup>Environmental Science Division, Argonne National Laboratory, Lemont, IL, USA.
<sup>2 </sup>Department of Civil, Environmental and Geospatial Engineering, Michigan Technological University, Houghton, MI, USA.
<sup>3 </sup>Pacific Northwest National Laboratory, Richland, WA, USA.

\* corresponding author:  wpringle@anl.gov

## Abstract
The transition from fall overturn to winter stratification influences subsequent ice formation and spring warming timing in deep dimictic lakes. We investigate the drivers of modeled deep-water temperature bias during this transition in Lake Superior using observations from two subsurface moorings and high-quality surface flux measurements. Large-eddy simulations (LES), a $k-\varepsilon$ turbulence-closure column model, and an analytical bottom-temperature model are evaluated over 60-day winter periods. At the 212-m Eastern Mooring, modeled temperature profiles and bottom-temperature evolution generally agree well with observations, although the hydrodynamic models, particularly LES, underestimate mixing. Including thermobaricity and accurate wind direction improves the simulated thermal structure and reduces deep-water warm bias. At the 384-m Southern Mooring, substantial warm biases persist, indicating that surface forcing and vertical mixing alone are insufficient. Horizontal processes associated with nearby steep bathymetry may be important. Results highlight the sensitivity of winter stratification to wind forcing, thermobaricity, and three-dimensional lake dynamics.

## Journal reference
_your journal reference_

## Code reference
- Pringle, W. J. (2026). Lake_Superior_Winter_Cooling [Software]. Zenodo. https://doi.org/some-doi-number

## Data reference

### Input data
- Austin, J., & Elmer, C (2023). Lake Superior moored temperature and currents, Sep 2005-May 2015 [Data set]. University Digital Conservancy. https://doi.org/10.13020/zqw9-mk81
- Spence, C., Lenters, J., & Nicholls, E. (2024). US-GL1: Stannard Rock, 2008 - 2022 [Data set]. Ameriflux. https://ameriflux.lbl.gov/sites/siteinfo/US-GL1

### Output data
- Pringle, W. J. (2026). My output dataset name [Data set]. Zenodo. https://doi.org/some-doi-number

## Contributing modeling software
| Model | Version | Github | Documentation |
|-------|---------|-----------------|-----|
| Oceananigans.jl | [0.95.3](https://github.com/CliMA/Oceananigans.jl/releases/tag/v0.95.3) | [Github](https://github.com/CliMA/Oceananigans.jl) | [HomePage](https://clima.github.io/OceananigansDocumentation/) |

## Reproduce my experiment
Fill in detailed info here or link to other documentation to thoroughly walkthrough how to use the contents of this repository to reproduce your experiment. Below is an example.


1. Install the software components required to conduct the experiment from [contributing modeling software](#contributing-modeling-software)
2. Download and install the supporting [input data](#input-data) required to conduct the experiment
3. Run the following scripts in the `workflow` directory to re-create this experiment:

| Script Name | Description | How to Run |
| --- | --- | --- |
| `step_one.py` | Script to run the first part of my experiment | `python3 step_one.py -f /path/to/inputdata/file_one.csv` |
| `step_two.py` | Script to run the second part of my experiment | `python3 step_two.py -o /path/to/my/outputdir` |

4. Download and unzip the [output data](#output-data) from my experiment 
5. Run the following scripts in the `workflow` directory to compare my outputs to those from the publication

| Script Name | Description | How to Run |
| --- | --- | --- |
| `compare.py` | Script to compare my outputs to the original | `python3 compare.py --orig /path/to/original/data.csv --new /path/to/new/data.csv` |

## Reproduce my figures
Use the scripts found in the `figures` directory to reproduce the figures used in this publication.

| Figure Number(s) | Script Name | Description | How to Run |
| --- | --- | --- | --- |
| 1, 2 | `generate_plot.py` | Description of figure, ie. "Plots the difference between our two scenarios" | `python3 generate_plot.py -input /path/to/inputs -output /path/to/outuptdir` |
| 3 | `generate_figure.py` | Description of figure, ie. "Shows how the mean and peak differences are calculated" | `python3 generate_figure.py -input /path/to/inputs -output /path/to/outuptdir` |

