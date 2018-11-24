# Qikiqtaruk Ecological Monitoring

_This repository contains code and data necessary to replicate data analysis, figures, and tables in __"Eighteen years of ecological monitoring reveals multiple lines of evidence for tundra vegetation change"__._

Authors: Isla H. Myers-Smith, Meagan M. Grabowski, Haydn J.D. Thomas, Sandra Angers-Blondin, Gergana N. Daskalova, Anne D. Bjorkman, Andrew M. Cunliffe, Jakob J. Assmann, Joseph Boyle, Edward McLeod, Samuel McLeod, Ricky Joe, Paden Lennie, Deon Arey, Richard Gordon, Cameron Eckert

### Contact: Isla Myers-Smith isla.myers-smith@ed.ac.uk

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
