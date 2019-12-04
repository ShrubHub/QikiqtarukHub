# Qikiqtaruk Ecological Monitoring

### Description
_This repository contains code and data necessary to replicate data analysis, figures, and tables in __"Eighteen years of ecological monitoring reveals multiple lines of evidence for tundra vegetation change"__._

### Paper reference
Myers‐Smith, I.H., Grabowski, M.M., Thomas, H.J., Angers‐Blondin, S., Daskalova, G.N., Bjorkman, A.D., Cunliffe, A.M., Assmann, J.J., Boyle, J.S., McLeod, E. and McLeod, S., 2019. Eighteen years of ecological monitoring reveals multiple lines of evidence for tundra vegetation change. Ecological Monographs, 89(2), p.e01351.
https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecm.1351

### Dataset DOI
https://doi.org/10.5281/zenodo.2397996

### Authors
Isla H. Myers-Smith, Meagan M. Grabowski, Haydn J.D. Thomas, Sandra Angers-Blondin, Gergana N. Daskalova, Anne D. Bjorkman, Andrew M. Cunliffe, Jakob J. Assmann, Joseph Boyle, Edward McLeod, Samuel McLeod, Ricky Joe, Paden Lennie, Deon Arey, Richard Gordon, Cameron Eckert

### Contact: Isla Myers-Smith isla.myers-smith@ed.ac.uk

### Data use guidelines
Data are publicly available using a Creative Commons Attribution 4.0 International copyright (CC BY 4.0). Data are fully public but should be appropriately referenced by citing the paper. Although not mandatory, we additionally suggest that data users contact and collaborate with data contributors if this dataset will contribute a substantial proportion of observations used in a particular paper or analysis.

### Data availability & access
This dataset will be maintained at this GitHub repository (https://github.com/ShrubHub/QikiqtarukHub).

### Citation
Myers-Smith IH, MM Grabowski, HJD Thomas, S Angers-Blondin, GN Daskalova, AD Bjorkman, AM Cunliffe, JJ Assmann, J Boyle, E McLeod, S McLeod, R Joe, P Lennie, D Arey, R Gordon, C Eckert. Eighteen years of ecological monitoring reveals multiple lines of evidence for tundra vegetation change. Ecological Monographs. In press.

### Acknowledgements from the manuscript
We thank the Herschel Island-Qikiqtaruk Territorial Park management, Catherine Kennedy, Dorothy Cooley, and Dr. Jill F. Johnstone for establishing and maintaining the phenology and composition data from Qikiqtaruk. We thank previous rangers including LeeJohn Meyook, Jordan McLeod, Pierre Foisy, Colin Gordon, Jeremy Hansen, Albert Rufus and field assistants including Santeri Lehtonen, William Palmer, Louise Beveridge, Clara Flintrop, John Godlee, Eleanor Walker, Catherine Henry and Anika Trimble. We thank Sigrid S. Nielsen and Prof. Christopher Burn for providing feedback on the manuscript. We thank the Alfred Wegener Institute, Helmholtz Centre for Polar and Marine Research for providing logistical support for this research and in particular Prof. Hugues Lantuit. Funding was provided by Herschel Island-Qikiqtaruk Territorial Park and the UK Natural Environment Research Council ShrubTundra Grant (NE/M016323/1), and we thank the NERC GEF for loan of GNSS equipment (NERC GEF:1063 and GEF:1064). Haydn Thomas and Jakob Assmann were funded by a NERC doctoral training partnership grant (NE/L002558/1). Sandra Angers-Blondin was funded by NSERC and the Canadian Centennial Scholarship Fund. Meagan Grabowski was funded by NSERC and Yukon Parks. We thank the Inuvialuit People for the opportunity to conduct research on their traditional lands.

### Please use the following acknowledgement statement if you use these data
We thank Isla H. Myers-Smith, Meagan M. Grabowski, Haydn J.D. Thomas, Sandra Angers-Blondin, Gergana N. Daskalova, Anne D. Bjorkman, Andrew M. Cunliffe, Jakob J. Assmann, Joseph Boyle, Edward McLeod, Samuel McLeod, Ricky Joe, Paden Lennie, Deon Arey, Richard Gordon, Cameron Eckert for collecting and compiling these data. We thank the Herschel Island-Qikiqtaruk Territorial Park management, Catherine Kennedy, Dorothy Cooley, and Dr. Jill F. Johnstone for establishing and maintaining the phenology and composition data from Qikiqtaruk. We thank previous rangers including LeeJohn Meyook, Jordan McLeod, Pierre Foisy, Colin Gordon, Jeremy Hansen, Albert Rufus and field assistants including Santeri Lehtonen, William Palmer, Louise Beveridge, Clara Flintrop, John Godlee, Eleanor Walker, Catherine Henry and Anika Trimble. We thank the Inuvialuit People for the opportunity to conduct research on their traditional lands.

### Other links
The phenology data are also included as a part of the 'Long-term phenology data for 47 tundra plant species at 19 high-latitude sites, 1992-2017' dataset:
https://www.polardata.ca/pdcsearch/?doi_id=12722

For more information on the Team Shrub research group see:
https://teamshrub.com/

For more information on this dataset and the perspectives of Qikiqtaruk Park Rangers please see:
https://teamshrub.com/2017/09/28/qikiqtaruk-perspectives-by-ranger-edward-mcleod/

# Data

All data for our analyses can be found in the `data` folder. The folder contains the following files:

```
+ qhi_SAC_accumulated.csv	
# Species accumulation data for the Herschel and Komakuk vegetation types

+ qhi_act_layer_all_data.csv
# Long-term active layer depth data

+ qhi_active_layer_2017_all.csv
# Active layer depth in the Herschel and Komakuk vegetation types from 2017

+ qhi_cover_ITEX_1999_2017.csv	
# Percentage cover data for different species within the long-term monitoring plots

+ qhi_cwa_1968_2017.csv	
# Sea ice concentration data

+ qhi_dendro_2017.csv	
# Dendroecological data

+ qhi_frs_2017.csv	
# Frost frequency data

+ qhi_herbivores.csv	
# Herbivore data

+ qhi_phen_with_before_2017.csv	
# Plant phenology data

+ qhi_point_fraiming_ITEX_1999-2017.csv	
# Species richness and community composition data

+ qhi_temp_2017.csv	
# Temperature data

+ qhi_temp_multimeter_2017.csv
Soil temperature data from boreholes

```

# Scripts

```
+ QHI_Monitoring_Code_Final.R
# Performs statistical analyses and plots figures
```

# Figures

The figures generated in `R` are stored in the `figures` folder.

# JAGS models

The `JAGS` models for the phenology models can be found in the `models` folder. The folder contains the follwing files:

```
+ growing_season_pheno.jags
+ int_cens_P2.jags
+ int_cens_P3.jags
+ int_cens_P5.jags
```

# Model outputs

Full model outputs for all statistical analyses are stored in the `model_outputs` folder:

```
+ 2017_mcmc_outputs.csv
# Model outputs for all models ran using the `MCMCglmm` package

+ 2017_standardised_effects.csv
# All standardised effect sizes

+ HE_plot_all_linear_outputs.csv
+ KO_plot_all_linear_outputs.csv
# Coefficients for species-specific changes in abundance over time in the Herschel and Komakuk vegetation types

+ Phenology_JAGS_coefficients_all_models.csv
# Coefficients from all JAGS phenology models

+ Phenology_JAGS_EffectSizes_all_models.csv
# Effect sizes from all JAGS phenology models
```

# Requirements

### Software
R version 3.3.3 or greater

### Packages
`MCMCglmm, ggplot2, plyr, gridExtra, dplyr, tidyr, rjags, R2jags, stargazer`

_Note that for full functionality of the code, the `plyr` package has to be loaded before the `dplyr` package._
