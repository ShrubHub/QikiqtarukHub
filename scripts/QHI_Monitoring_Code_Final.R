# Qikiqtaruk Ecological Monitoring manuscript script
# Code for all modelling and data visualisation within the manuscript

# Written by Isla Myers-Smith, Anne Bjorkman, Haydn Thomas, Sandra Angers-Blondin and Gergana Daskalova
# E-mail: isla.myers-smith@ed.ac.uk
# 2018-01-30 

# Packages ----
library(MCMCglmm)
library(ggplot2)
library(plyr) # Load plyr before dplyr
library(gridExtra)
library(dplyr)
library(tidyr)
library(rjags)
library(R2jags)
library(stargazer)

# Defining functions used within the script ----

# Function to extract MCMCglmm model summary outputs ----
clean.MCMC <- function(x) {
  sols <- summary(x)$solutions  # pull out relevant info from model summary
  Gcovs <- summary(x)$Gcovariances
  Rcovs <- summary(x)$Rcovariances
  
  fixed <- data.frame(row.names(sols), sols, row.names = NULL)  # convert to dataframes with the row.names as the first col
  random <- data.frame(row.names(Gcovs), Gcovs, row.names = NULL)
  residual <- data.frame(row.names(Rcovs), Rcovs, row.names = NULL)
  
  names(fixed)[names(fixed) == "row.names.sols."] <- "variable"  # change the columns names to variable, so they all match
  names(random)[names(random) == "row.names.Gcovs."] <- "variable"
  names(residual)[names(residual) == "row.names.Rcovs."] <- "variable"
  
  fixed$effect <- "fixed"  # add ID column for type of effect (fixed, random, residual)
  random$effect <- "random"
  residual$effect <- "residual"
  
  modelTerms <- as.data.frame(bind_rows(fixed, random, residual))  # merge it all together
}

# Function for when models have no random effects
clean.MCMC.2 <- function(x) {
  sols <- summary(x)$solutions  # pull out relevant info from model summary
  Gcovs <- summary(x)$Gcovariances
  Rcovs <- summary(x)$Rcovariances
  
  fixed <- data.frame(row.names(sols), sols, row.names = NULL)  # convert to dataframes with the row.names as the first col
  residual <- data.frame(row.names(Rcovs), Rcovs, row.names = NULL)
  
  names(fixed)[names(fixed) == "row.names.sols."] <- "variable"  # change the columns names to variable, so they all match
  names(residual)[names(residual) == "row.names.Rcovs."] <- "variable"
  
  fixed$effect <- "fixed"  # add ID column for type of effect (fixed, random, residual)
  residual$effect <- "residual"
  
  modelTerms <- as.data.frame(bind_rows(fixed, residual))  # merge it all together
}

getName.MCMC <- function(x) deparse(substitute(x))  # adding the model name

# Customised ggplot2 theme function ----
theme_QHI <- function(){
  theme_bw() +
    theme(axis.text = element_text(size = 16), 
          axis.title = element_text(size = 20),
          axis.text.x = element_text(angle = -45, hjust = -0.05),
          axis.line.x = element_line(color = "black"), 
          axis.line.y = element_line(color = "black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 20, vjust = 1, hjust = 0),
          legend.text = element_text(size = 16, face = "italic"),          
          legend.title = element_blank(),                              
          legend.position = c(0.9, 0.9), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 4, linetype = "blank"))
}

# Loading data ----

# Temperature - CRU DATA
qhi.tmp <- read.csv("Qikiqtaruk_manuscript/data/qhi_temp_2017.csv")

# Frost frequency days 
qhi.frs <- read.csv("Qikiqtaruk_manuscript/data/qhi_frs_2017.csv")
cutoff <- 1979
qhi.frs <- filter(qhi.frs, year > cutoff)
# Summing frost frequency days
qhi.frs.sum <- ddply(qhi.frs, "year", summarize, sum.ice = sum(value))

# Sea Ice Concentration
cwa <- read.csv("Qikiqtaruk_manuscript/data/qhi_cwa_1968_2017.csv")

cwa$Year <- as.numeric(cwa$Year)
cwa.sub <- cwa[cwa$Year > 1979,]

# Temperature at depth 12m and 15m
boreholes <- read.csv("Qikiqtaruk_manuscript/data/qhi_temp_multimeter_2017.csv")
soil.temp.qhi <- boreholes %>% filter(Depth == "12" | Depth == "15" | Depth == "16")  
date_split <- strsplit(as.character(soil.temp.qhi$date),'-')
date_split <- do.call(rbind, date_split)
soil.temp.qhi$year <- as.numeric(date_split[,3]) + 2000
# remove outlier point for first borehole 15m in 2012
soil.temp.qhi <- filter(soil.temp.qhi, temperature > -9.9 )

# Phenology
qiphen <- read.csv("Qikiqtaruk_manuscript/data/qhi_phen_with_before_2017.csv", stringsAsFactors = F)

# Snowmelt
snow <- dplyr::select(qiphen[!is.na(qiphen$P1),], 1:4)

# Radial growth
dendro <- read.csv("Qikiqtaruk_manuscript/data/qhi_dendro_2017.csv", stringsAsFactors = F)

# Vegetation cover and canopy height ----

HEcover <- read.csv("Qikiqtaruk_manuscript/data/qhi_cover_ITEX_1999_2017.csv")
levels(HEcover$name) # Checking all spelling is correct

# Zero fill data
HEcover$cover <- as.numeric(as.character(HEcover$cover))

HEcover <- HEcover %>% 
  complete(plot, sub_name, year, name, fill = list(cover = 0)) %>% 
  group_by(name, year, sub_name, plot) %>%
  summarise(cover = mean(cover))

# Total biomass - Cover
species_IDs <- as.data.frame(unique(HEcover$name))
HEcoverveg <- subset(HEcover, name !="XXXlitter" & name !="XXXlitter " & 
                       name !="XXXbareground" & name !="XXXbareground " & 
                       name !="XXXrock" & name !="XXXrock " & name !="XXXfeces" & 
                       name !="XXXfeces " & name !="XXXstandingwater" & 
                       name !="XXXstandingwater " & name !="XXXspider")

biomass_cover <- ddply(HEcoverveg, .(year, plot, sub_name), summarise,
                       Biomass = sum(cover))

HEcoverHE <- subset(HEcover, sub_name == "QHI:HE")
HEcoverKO <- subset(HEcover, sub_name == "QHI:KO")

HEcoverHE$name <- as.character(HEcoverHE$name)
HEcoverHE$name[HEcoverHE$name == "XXXcarex:QHI "] <- "Carex sp."

# Total biomass - Hits
abundance <- read.csv("Qikiqtaruk_manuscript/data/qhi_point_fraiming_ITEX_1999-2017.csv")
abundance$unique_coords <- paste(abundance$X, abundance$Y, sep = "")

# Check no. points per plot as sometimes 90
n_hits <- ddply(abundance, .(YEAR, PLOT, SUBSITE), summarise,
                n_hits = length(unique(unique_coords)))

# Remove non-veg
abundance_veg <- subset(abundance, SPP !="XXXlitter" & SPP !="XXXlitter " & 
                          SPP !="XXXbareground" & SPP !="XXXbareground " & 
                          SPP !="XXXrock" & SPP !="XXXrock " & SPP !="XXXfeces" & 
                          SPP !="XXXfeces " & SPP !="XXXstandingwater" & 
                          SPP !="XXXstandingwater " & SPP !="XXXspider")

# Calculate biomass per plot
biomass_hits <- ddply(abundance_veg, .(YEAR, PLOT, SUBSITE), summarise,
                      Biomass = sum(Abundance))

# Combine with number of points per plot
biomass_hits <- merge(biomass_hits, n_hits)
biomass_hits$Biomass <- biomass_hits$Biomass/biomass_hits$n_hits

# ** Bareground cover data ----

abundance$uniqueID <- paste(abundance$YEAR, abundance$SUBSITE, 
                            abundance$PLOT, abundance$X, abundance$Y, sep="")  # Assign unique ID to every point
bareground_IDs <- abundance[abundance$SPP=="XXXlitter" | abundance$SPP=="XXXlitter " | 
                              abundance$SPP=="XXXbareground" | abundance$SPP=="XXXbareground " | 
                              abundance$SPP=="XXXrock" | abundance$SPP =="XXXrock " | 
                              abundance$SPP =="XXXfeces" | abundance$SPP =="XXXfeces " | 
                              abundance$SPP =="XXXstandingwater" | abundance$SPP =="XXXstandingwater " | 
                              abundance$SPP =="XXXspider",17]  # Identify points with non-vegetation indicators
abundance_BGs <- abundance[abundance$uniqueID %in% bareground_IDs,]  # Extract only points that have non-veg indicators

Out=NULL  # Set up loop
for(i in unique(abundance_BGs$uniqueID)){  # For each point
  a <- subset(abundance_BGs, uniqueID==i)  # create dataframe of all entries for that point
  b <- a[a$SPP == "XXXlitter" | a$SPP == "XXXlitter " | a$SPP == "XXXbareground" | 
           a$SPP == "XXXbareground " | a$SPP == "XXXrock" | a$SPP == "XXXrock " | 
           a$SPP == "XXXfeces" | a$SPP == "XXXfeces " | a$SPP == "XXXstandingwater" | 
           a$SPP == "XXXstandingwater " | a$SPP == "XXXspider",]  # Identify how many entries are not vegetation
  c <- nrow(a) - nrow(b)  # Finnd out if any entries are vegetation (i.e. total rows - non-veg rows)
  Out <- rbind(Out, c(i, c))  # Extract point name and number vegetation entries
}

BGs <- as.data.frame(Out)  # Convert into data frame
BGs <- subset(BGs,V2=="0")  # Extract points for which there are only non-veg data (i.e.i.e. total rows - non-veg rows = 0)
BGs <- BGs[,1]  # Extract only first column (unique points)

bareground <- abundance_BGs <- abundance[abundance$uniqueID %in% BGs,]  # Create dataframe of bare ground points 
bareground <- ddply(bareground,.(YEAR, PLOT, SUBSITE), summarise,
                    Bareground = sum(Abundance))  # Count number of bare ground points per plot

# Create dummy dataframe with all plots because if plots have no bare ground they wont be included
bareground_full <- ddply(abundance,.(YEAR, SUBSITE, PLOT), summarise,
                         Bareground = 0)  

# Replace dummy bareground with real baregound
bareground_full$Bareground <- bareground$Bareground[match(paste(bareground_full$YEAR,
                                                                bareground_full$SUBSITE,
                                                                bareground_full$PLOT),
                                                          paste(bareground$YEAR,
                                                                bareground$SUBSITE,
                                                                bareground$PLOT))] 
bareground <- bareground_full  # Rename to original
bareground[is.na(bareground$Bareground),]$Bareground <- 0  # Replace NAs from match with zeros

# ** Canopy height ----

abundance$Height..cm. <- as.numeric(as.character(abundance$Height..cm.))
heights <- subset(abundance,!is.na(Height..cm.))

canopy <- ddply(heights,.(YEAR, SUBSITE, PLOT, unique_coords), summarise,
                MaxHeight = max(Height..cm.)) #Take max height at each unique point

# Average by plot
avg_heights <- ddply(canopy,.(YEAR, SUBSITE, PLOT), summarise,
                     Mean.C.H = mean(MaxHeight))

# ** Salix pulchra Height ----

salpuls <- subset(abundance, SPP=="SALPUL"|SPP == "Salix pulchra") #Take only heights for Salix pulchra

salpuls <- subset(salpuls,!is.na(Height..cm.))

meansp <- ddply(salpuls,.(YEAR), summarise,
                mean.height = mean(Height..cm.),
                sd = sd(Height..cm.))

avg_salpuls <- ddply(salpuls,.(YEAR, SUBSITE, PLOT), summarise,
                     Mean.C.H = mean(Height..cm.))

# ** Community composition and diversity measures ----

# Remove rows with no cover
diversity <- subset(HEcover, cover>0)
# Remove non veg
diversity <- subset(diversity, name !="XXXlitter" & name !="XXXlitter " & 
                 name !="XXXbareground" & name !="XXXbareground " & 
                 name !="XXXrock" & name !="XXXrock " & 
                 name !="XXXfeces" & name !="XXXfeces " & 
                 name !="XXXstandingwater" & name !="XXXstandingwater " & 
                 name !="XXXspider"& name !="Xxxspider"& name !="XXXspider ")
# Add unique plots
diversity$plot_unique <- paste(diversity$sub_name,diversity$plot,diversity$year,sep="")

# Convert to relative cover
plot_cover <- ddply(diversity,.(plot_unique), summarise,
                    total_cover = sum(cover))
diversity$total_cover <- plot_cover$total_cover[match(diversity$plot_unique, plot_cover$plot_unique)]
diversity$rel_cover <- diversity$cover/diversity$total_cover*100

# ** Species pool data ----

SAC.acc <- read.csv("Qikiqtaruk_manuscript/data/qhi_SAC_accumulated.csv")

# ** Active layer depth data ----

act.layer <- read.csv("Qikiqtaruk_manuscript/data/qhi_act_layer_all_data.csv")
act.layer.2017 <- read.csv("Qikiqtaruk_manuscript/data/qhi_active_layer_2017_all.csv")

# Modelling and data visualisation ----

# Defining parameter-expanded priors for MCMCglmm models

# For models with one random effect
prior1 <- list(R = list(V = 1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 1, 
                                  alpha.mu = 0, alpha.v = 10000)))

# For models with two random effects
prior2 <- list(R = list(V = 1, nu = 0.002), 
               G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000), 
                        G2 = list(V = 1, nu = 1, alpha.mu = 0, alpha.v = 10000)))

# Figure 3. Climate, sea ice concentration and soil temperatures ----

# Temperature

# Data manipulation
# Use season instead of year
qhi.tmp$season <- ifelse(qhi.tmp$month %in% c(10:12), qhi.tmp$year + 1, qhi.tmp$year)

# Calculate mean per season
qhi.tmp.spring <- ddply(qhi.tmp[qhi.tmp$month %in% c(4, 5),], "season", summarize, meanT = mean(temp, na.rm = T))
head(qhi.tmp.spring)
qhi.tmp.spring$Time <- "Spring"

qhi.tmp.summer <- ddply(qhi.tmp[qhi.tmp$month %in% c(6, 7),], "season", summarize, meanT = mean(temp, na.rm = T))
head(qhi.tmp.summer)
qhi.tmp.summer$Time <- "Summer"

qhi.tmp.fall <- ddply(qhi.tmp[qhi.tmp$month %in% c(8, 9),], "season", summarize, meanT = mean(temp, na.rm = T))
head(qhi.tmp.fall)
qhi.tmp.fall$Time <- "Fall"

qhi.tmp.winter <- ddply(qhi.tmp[qhi.tmp$month %in% c(10 ,11, 12, 1, 2, 3),], "season", summarize, meanT = mean(temp, na.rm = T))
head(qhi.tmp.winter)
qhi.tmp.winter$Time <- "Winter"

qhi.tmp.all <- rbind(qhi.tmp.summer, qhi.tmp.winter, qhi.tmp.fall, qhi.tmp.spring)

# Filtering data to include records from 1979 onwards
qhi.tmp.spring.2 <- filter(qhi.tmp.spring, season > cutoff)
qhi.tmp.summer.2 <- filter(qhi.tmp.summer, season > cutoff)
qhi.tmp.fall.2 <- filter(qhi.tmp.fall, season > cutoff)
qhi.tmp.winter.2 <- filter(qhi.tmp.winter, season > cutoff)

# Modelling
# Including I(season - 1979) to transform year records from e.g. 1980, 1981 to 1, 2, 3, so that the model
# doesn't calculate predictions for all years leading up the 1980s

# Spring
spring_m <- MCMCglmm(meanT ~ I(season - 2000), data = qhi.tmp.spring.2, 
                     family = "gaussian", pr=TRUE, nitt = 100000, burnin = 20000)
summary(spring_m)

# Assessing model convergence
# plot(spring_m$VCV)
# plot(spring_m$Sol)
autocorr(spring_m$VCV)

# Summer
summer_m <- MCMCglmm(meanT ~ I(season - 2000), data = qhi.tmp.summer.2, 
                     family = "gaussian", pr = TRUE, nitt = 100000, burnin = 20000)
summary(summer_m)
# plot(summer_m$VCV)
# plot(summer_m$Sol)
autocorr(summer_m$VCV)

# Fall
fall_m <- MCMCglmm(meanT ~ I(season - 2000), data = qhi.tmp.fall.2, family = "gaussian", 
                   pr = TRUE, nitt = 100000, burnin = 20000)
summary(fall_m)
# plot(fall_m$VCV)
# plot(fall_m$Sol)
autocorr(fall_m$VCV)

# Winter
winter_m <- MCMCglmm(meanT ~ I(season - 2000), data = qhi.tmp.winter.2, family = "gaussian",
                     pr = TRUE, nitt = 100000, burnin = 20000)
summary(winter_m)
# plot(winter_m$VCV)
# plot(winter_m$Sol)
autocorr(winter_m$VCV)

# Temperature plots

# Spring
# Calculating model predictions
nyears <- 17
niter <- length(spring_m$Sol[,"(Intercept)"])
spring_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    spring_preds[i,j] <- spring_m$Sol[i,"(Intercept)"] + spring_m$Sol[i,"I(season - 2000)"]*j
  }
}

spring_preds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  spring_preds_df[i,] <- quantile(spring_preds[,i], c(0.025, 0.5, 0.975))
}

spring_preds_df <- cbind.data.frame(lower = spring_preds_df[,1], 
                                    mean = spring_preds_df[,2], upper = spring_preds_df[,3], year = seq(1:17))

# Spring graph
(spring <- ggplot() +
    geom_point(data = qhi.tmp.spring.2, aes(x = season, y = meanT), alpha = 0.8, colour = "#b70000", size = 4) +
    scale_y_continuous(limits = c(-14, -3)) + 
    scale_x_continuous(limits = c(2001, 2017), breaks = c(2001, 2006, 2011, 2017)) +
    geom_ribbon(data = spring_preds_df, aes(x = year + 2000, ymin = lower, ymax = upper), 
                fill = "#b70000", alpha = 0.2) +
    geom_line(data = spring_preds_df, aes(x = year + 2000, y = mean), colour = "#b70000") +
  theme_QHI() +
    labs(x = "", y = "Temperature (°C)\n", title = "(a) Spring\n") + 
  annotate("text", x = 2002, y = -3.2, label = "April - May", size = 5, hjust = 0))

# Summer
# Calculating model predictions
nyears <- 17
niter <- length(summer_m$Sol[,"(Intercept)"])
summer_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    summer_preds[i,j] <- summer_m$Sol[i,"(Intercept)"] + summer_m$Sol[i,"I(season - 2000)"]*j
  }
}

summer_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  summer_preds_df[i,] <- quantile(summer_preds[,i], c(0.025, 0.5, 0.975))
}

summer_preds_df <- cbind.data.frame(lower = summer_preds_df[,1], 
                                    mean = summer_preds_df[,2], upper = summer_preds_df[,3], year = seq(1:17))

# Summer graph
(summer <- ggplot() +
    geom_point(data = qhi.tmp.summer.2, aes(x = season, y = meanT), alpha = 0.8, colour = "#b70000", size = 4) +
    scale_y_continuous(limits = c(0, 14)) + 
    scale_x_continuous(limits = c(2001, 2017), breaks = c(2001, 2006, 2011, 2017)) +
    geom_ribbon(data = summer_preds_df, aes(x = year + 2000, ymin = lower, ymax = upper), 
                fill = "#b70000", alpha = 0.2) +
    geom_line(data = summer_preds_df, aes(x = year + 2000, y = mean), colour = "#b70000") +
    theme_QHI() +
    labs(x = "", y = "Temperature (°C)\n", title = "(b) Summer\n") +
  annotate("text", x = 2002, y = 13.7, label = "June - July", size = 5, hjust = 0))

# Fall
# Calculating model predictions
nyears <- 17
niter <- length(fall_m$Sol[,"(Intercept)"])
fall_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    fall_preds[i,j] <- fall_m$Sol[i,"(Intercept)"] + fall_m$Sol[i,"I(season - 2000)"]*j
  }
}

fall_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  fall_preds_df[i,] <- quantile(fall_preds[,i], c(0.025, 0.5, 0.975))
}

fall_preds_df <- cbind.data.frame(lower = fall_preds_df[,1], 
                                  mean = fall_preds_df[,2], upper = fall_preds_df[,3], year = seq(1:17))

# Fall graph
(fall <- ggplot() +
    geom_point(data = qhi.tmp.fall.2, aes(x = season, y = meanT), alpha = 0.8, colour = "#b70000", size = 4) +
    scale_y_continuous(limits = c(2, 10)) + 
    scale_x_continuous(limits = c(2001, 2017), breaks = c(2001, 2006, 2011, 2017)) +
    geom_ribbon(data = fall_preds_df, aes(x = year + 2000, ymin = lower, ymax = upper), 
                fill = "#b70000", alpha = 0.2) +
    geom_line(data = fall_preds_df, aes(x = year + 2000, y = mean), colour = "#b70000") +
    theme_QHI() +
    labs(x = "", y = "Temperature (°C)\n", title = "(c) Fall\n") +
    annotate("text", x = 2002, y = 9.8, label = "August - September", size = 5, hjust = 0))

# Winter
# Calculating model predictions
nyears <- 17
niter <- length(winter_m$Sol[,"(Intercept)"])

winter_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    winter_preds[i,j] <- winter_m$Sol[i,"(Intercept)"] + winter_m$Sol[i,"I(season - 2000)"]*j
  }
}

winter_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  winter_preds_df[i,] <- quantile(winter_preds[,i], c(0.025, 0.5, 0.975))
}

winter_preds_df <- cbind.data.frame(lower = winter_preds_df[,1], 
                                    mean = winter_preds_df[,2], upper = winter_preds_df[,3], year = seq(1:17))

# Winter graph
(winter <- ggplot() +
    geom_point(data = qhi.tmp.winter.2, aes(x = season, y = meanT), alpha = 0.8, colour = "#b70000", size = 4) +
    scale_y_continuous(limits = c(-25, -10)) + 
    scale_x_continuous(limits = c(2001, 2017), breaks = c(2001, 2006, 2011, 2017)) +
    geom_ribbon(data = winter_preds_df, aes(x = year + 2000, ymin = lower, ymax = upper), 
                fill = "#b70000", alpha = 0.2) +
    geom_line(data = winter_preds_df, aes(x = year + 2000, y = mean), colour = "#b70000") +
    theme_QHI() +
    labs(x = "", y = "Temperature (°C)\n", title = "(d) Winter\n") +
    annotate("text", x = 2002, y = -10, label = "October - March", size = 5, hjust = 0, vjust=0))

# Frost frequency MCMCglmm model
frs_m <- MCMCglmm(sum.ice ~ I(year - 1979), data = qhi.frs.sum, family = "gaussian", 
                  pr=TRUE, nitt = 100000, burnin = 20000)
summary(frs_m)
# plot(frs_m$VCV)
# plot(frs_m$Sol)
autocorr(frs_m$VCV)

# Calculating model predictions
nyears <- 36
niter <- length(frs_m$Sol[,"(Intercept)"])
frs_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    frs_preds[i,j] <- frs_m$Sol[i,"(Intercept)"] + frs_m$Sol[i,"I(year - 1979)"]*j
  }
}

frs_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  frs_preds_df[i,] <- quantile(frs_preds[,i], c(0.025, 0.5, 0.975))
}

frs_preds_df <- cbind.data.frame(lower = frs_preds_df[,1], 
                                 mean = frs_preds_df[,2], upper = frs_preds_df[,3], 
                                 year = seq(1:36))

# Frost frequency days graph
(frost <- ggplot() +
    geom_point(data = qhi.frs.sum, aes(x = year, y = sum.ice), alpha = 0.8, colour = "#3cd0ea", size = 4) +
    scale_x_continuous(limits = c(1979, 2017), breaks = c(1979, 1989, 1999, 2009, 2017)) +
    geom_ribbon(data = frs_preds_df, aes(x = year + 1979, ymin = lower, ymax = upper), 
                fill = "#3cd0ea", alpha = 0.2) +
    geom_line(data = frs_preds_df, aes(x = year + 1979, y = mean), colour = "#3cd0ea") +
    theme_QHI() +
    labs(x = "", y = "Frost frequency\n", title = "(e) Frost frequency\n"))

# Snow melt MCMCglmm model
# Include 'SPP' as a random effect as it represents the transect where the data were collected
snow_m <- MCMCglmm(P1 ~ I(Year - 2000), random = ~ Spp + Year, data = snow, 
                   family = "gaussian", pr=TRUE, nitt = 100000, burnin = 20000, prior = prior2)
summary(snow_m)
# plot(snow_m$VCV)
# plot(snow_m$Sol)
autocorr(snow_m$VCV)

# Calculating model predictions
nyears <- 17
niter <- length(snow_m$Sol[,"(Intercept)"])

snow_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    snow_preds[i,j] <- snow_m$Sol[i,"(Intercept)"] + snow_m$Sol[i,"I(Year - 2000)"]*j
  }
}

snow_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  snow_preds_df[i,] <- quantile(snow_preds[,i], c(0.025, 0.5, 0.975))
}

snow_preds_df <- cbind.data.frame(lower = snow_preds_df[,1], 
                                  mean = snow_preds_df[,2], upper = snow_preds_df[,3], 
                                  year = seq(1:17))

# Snow melt graph
(snow.melt <- ggplot() +
    geom_point(data = snow, aes(x = Year, y = P1), alpha = 0.8, colour = "#3cd0ea", size = 4) +
    scale_x_continuous(limits = c(2001, 2017), breaks = c(2001, 2006, 2011, 2017)) +
    geom_ribbon(data = snow_preds_df, aes(x = year + 2000, ymin = lower, ymax = upper), 
                fill = "#3cd0ea", alpha = 0.2) +
    geom_line(data = snow_preds_df, aes(x = year + 2000, y = mean), colour = "#3cd0ea") +
    theme_QHI() +
    labs(x = "", y = "Day of the year\n", title = "(f) Snow melt date\n"))

# Sea ice concentration MCMCglmm model

# Creating a column of trial hits and misses for the binomial model
cwa.sub$failures <- round(100 - 100*cwa.sub$Total.Concentration)

ice_m <- MCMCglmm(cbind(round(100*(Total.Concentration)), failures) ~ I(Year-1979), 
                  data = cwa.sub, family="multinomial2", pr=TRUE, nitt = 100000, burnin = 20000)
summary(ice_m)
# plot(ice_m$VCV)
# plot(ice_m$Sol)
autocorr(ice_m$VCV)

# Calculating model predictions
nyears <- 38
niter <- length(ice_m$Sol[,"(Intercept)"])
ice_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    ice_preds[i,j] <- 100*plogis(ice_m$Sol[i,"(Intercept)"] + ice_m$Sol[i,"I(Year - 1979)"]*j)
  }
}

ice_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  ice_preds_df[i,] <- quantile(ice_preds[,i], c(0.025, 0.5, 0.975))
}

ice_preds_df <- cbind.data.frame(lower = ice_preds_df[,1], 
                                 mean = ice_preds_df[,2], upper = ice_preds_df[,3], 
                                 year = seq(1:38))

# Minimum Sea Ice Concentration graph
(sea.ice <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "darkgrey") +
    geom_line(data = cwa.sub, aes(x = Year, y = Total.Concentration), size = 2, alpha = 0.4, colour = "#3cd0ea") +
    scale_x_continuous(limits = c(1979, 2017), breaks = c(1979, 1989, 1999, 2009, 2017)) +
    geom_ribbon(data = ice_preds_df, aes(x = year + 1979, ymin = lower/100, ymax = upper/100), 
                fill = "#3cd0ea", alpha = 0.2) +
    geom_line(data = ice_preds_df, aes(x = year + 1979, y = mean/100), colour = "#3cd0ea") +
    theme_QHI() +
    labs(x = "", y = "Minimum sea ice concentration\n", title = "(g) Sea ice concentration\n"))

# Temperature at 12 m depth model for 15 m borehole

# Making a column to use as a legend in the figure
soil.temp.qhi$id <- paste(soil.temp.qhi$Bore.hole, soil.temp.qhi$Depth)
soil.temp.qhi$id <- recode(soil.temp.qhi$id, "15 m 12" = "First borehole 12 m")
soil.temp.qhi$id <- recode(soil.temp.qhi$id, "15 m 15" = "First borehole 15 m")
soil.temp.qhi$id <- recode(soil.temp.qhi$id, "43 m 12" = "Second borehole 12 m")
soil.temp.qhi$id <- recode(soil.temp.qhi$id, "43 m 16" = "Second borehole 16 m")
soil.temp.qhi$id <- as.factor(soil.temp.qhi$id)
soil.temp.qhi$id <- factor(soil.temp.qhi$id, 
                                 levels = c("First borehole 12 m", "First borehole 15 m",
                                            "Second borehole 12 m", "Second borehole 16 m"),
                                 labels=c("First borehole 12 m", "First borehole 15 m", 
                                          "Second borehole 12 m", "Second borehole 16 m"))
soil.temp.qhi <- soil.temp.qhi %>% select(Bore.hole, Depth, date, 
                                          temperature, year, id)
soil.temp.qhi <- na.omit(soil.temp.qhi)

soil.temp.qhi.sum <- soil.temp.qhi %>% group_by(id, year) %>% 
  summarise(minTemp = min(temperature)) %>% as.data.frame()

soil.temp.qhi15.12 <- filter(soil.temp.qhi.sum, id == "First borehole 12 m")
soil.temp.qhi15.12$squared <- (soil.temp.qhi15.12$year - 2000)^2
soil.temp.qhi15.15 <- filter(soil.temp.qhi.sum, id == "First borehole 15 m")
soil.temp.qhi15.15$squared <- (soil.temp.qhi15.15$year - 2000)^2
soil.temp.qhi43.12 <- filter(soil.temp.qhi.sum, id == "Second borehole 12 m")
soil.temp.qhi43.12$squared <- (soil.temp.qhi43.12$year - 2004)^2
soil.temp.qhi43.16 <- filter(soil.temp.qhi.sum, id == "Second borehole 16 m")
soil.temp.qhi43.16$squared <- (soil.temp.qhi43.16$year - 2004)^2

# Temperature at 12 m depth model for 15 m borehole
soil15_12_m <- MCMCglmm(minTemp ~ squared, data = soil.temp.qhi15.12, family = "gaussian", 
                        pr = TRUE, nitt = 100000, burnin = 20000)
summary(soil15_12_m)
# plot(soil15_12_m$VCV)
# plot(soil15_12_m$Sol)
autocorr(soil15_12_m$VCV)

# Calculating model predictions
nyears <- 17
niter <- length(soil15_12_m$Sol[,"(Intercept)"])
soil15_12_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    soil15_12_preds[i,j] <- soil15_12_m$Sol[i,"(Intercept)"] + soil15_12_m$Sol[i,"squared"]*j^2
  }
}

soil15_12_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  soil15_12_preds_df[i,] <- quantile(soil15_12_preds[,i], c(0.025, 0.5, 0.975))
}

soil15_12_preds_df <- cbind.data.frame(lower = soil15_12_preds_df[,1], 
                                       mean = soil15_12_preds_df[,2], 
                                       upper = soil15_12_preds_df[,3], 
                                       year = seq(1:17))

# Temperature at 15 m depth model for 15 m borehole
soil15_15_m <- MCMCglmm(minTemp ~ squared, data = soil.temp.qhi15.15, family = "gaussian", pr = TRUE, nitt = 100000, burnin = 20000)
summary(soil15_15_m)
# plot(soil15_15_m$VCV)
# plot(soil15_15_m$Sol)
autocorr(soil15_15_m$VCV)

# Calculating model predictions
nyears <- 17
niter <- length(soil15_15_m$Sol[,"(Intercept)"])
soil15_15_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    soil15_15_preds[i,j] <- soil15_15_m$Sol[i,"(Intercept)"] + soil15_15_m$Sol[i,"squared"]*j^2
  }
}

soil15_15_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  soil15_15_preds_df[i,] <- quantile(soil15_15_preds[,i], c(0.025, 0.5, 0.975))
}

soil15_15_preds_df <- cbind.data.frame(lower = soil15_15_preds_df[,1], 
                                       mean = soil15_15_preds_df[,2], 
                                       upper = soil15_15_preds_df[,3], 
                                       year = seq(1:17))

# Temperature at 12 m depth model for 43 m borehole
soil43_12_m <- MCMCglmm(minTemp ~ squared, data = soil.temp.qhi43.12, family = "gaussian", pr = TRUE, nitt = 100000, burnin = 20000)
summary(soil43_12_m)
# plot(soil43_12_m$VCV)
# plot(soil43_12_m$Sol)
autocorr(soil43_12_m$VCV)

# Calculating model predictions
nyears <- 13
niter <- length(soil43_12_m$Sol[,"(Intercept)"])
soil43_12_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    soil43_12_preds[i,j] <- soil43_12_m$Sol[i,"(Intercept)"] + soil43_12_m$Sol[i,"squared"]*j^2
  }
}

soil43_12_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  soil43_12_preds_df[i,] <- quantile(soil43_12_preds[,i], c(0.025, 0.5, 0.975))
}

soil43_12_preds_df <- cbind.data.frame(lower = soil43_12_preds_df[,1], 
                                       mean = soil43_12_preds_df[,2], 
                                       upper = soil43_12_preds_df[,3], 
                                       year = seq(5:17))
soil43_12_preds_df$year <- soil43_12_preds_df$year + 4

# Temperature at 16 m depth model for 43 m borehole
soil43_16_m <- MCMCglmm(minTemp ~ squared, data = soil.temp.qhi43.16, family = "gaussian", pr = TRUE, nitt = 100000, burnin = 20000)
summary(soil43_16_m)
# plot(soil43_16_m$VCV)
# plot(soil43_16_m$Sol)
autocorr(soil43_16_m$VCV)

# Calculating model predictions
nyears <- 13
niter <- length(soil43_16_m$Sol[,"(Intercept)"])
soil43_16_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    soil43_16_preds[i,j] <- soil43_16_m$Sol[i,"(Intercept)"] + soil43_16_m$Sol[i,"squared"]*j^2
  }
}

soil43_16_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  soil43_16_preds_df[i,] <- quantile(soil43_16_preds[,i], c(0.025, 0.5, 0.975))
}

soil43_16_preds_df <- cbind.data.frame(lower = soil43_16_preds_df[,1], 
                                       mean = soil43_16_preds_df[,2], 
                                       upper = soil43_16_preds_df[,3], 
                                       year = seq(5:17))
soil43_16_preds_df$year <- soil43_16_preds_df$year + 4

# Temperature at 12 m depth graph
(soil.qhi <- ggplot() +
    geom_point(data = soil.temp.qhi.sum, aes(x = year, y = minTemp, colour = id), alpha = 0.8, size = 4) +
    scale_colour_manual(values = c("#3cd0ea", "#005463", "#4886d6", "#002e6b")) +
    geom_vline(xintercept = 2005, linetype = 2, colour = "gray20") +
    geom_ribbon(data = soil15_12_preds_df, aes(x = year + 2000, ymin = lower, ymax = upper), 
                fill = "#3cd0ea", alpha = 0.2) +
    geom_line(data = soil15_12_preds_df, aes(x = year + 2000, y = mean), colour = "#3cd0ea") +
    geom_ribbon(data = soil15_15_preds_df, aes(x = year + 2000, ymin = lower, ymax = upper), 
                fill = "#005463", alpha = 0.2) +
    geom_line(data = soil15_15_preds_df, aes(x = year + 2000, y = mean), colour = "#005463") +
    geom_ribbon(data = soil43_12_preds_df, aes(x = year + 2000, ymin = lower, ymax = upper), 
                fill = "#4886d6", alpha = 0.2) +
    geom_line(data = soil43_12_preds_df, aes(x = year + 2000, y = mean), colour = "#4886d6") +
    geom_ribbon(data = soil43_16_preds_df, aes(x = year + 2000, ymin = lower, ymax = upper), 
                fill = "#002e6b", alpha = 0.2) +
    geom_line(data = soil43_16_preds_df, aes(x = year + 2000, y = mean), colour = "#002e6b") +
    scale_x_continuous(limits = c(2001, 2017), breaks = c(2001, 2006, 2011, 2017)) +
    theme_QHI() +
    theme(legend.position = c(0.35, 0.90),
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5),
          legend.background = element_rect(fill = "white")) +
    labs(x = "", y = "Temperature (°C)\n", title = "(h) Minimum soil temperature\n"))

# Arranging in a panel and saving the file
Figure3 <- grid.arrange(spring, summer, fall, winter, frost, snow.melt, sea.ice, soil.qhi, ncol=2)

ggsave("Qikiqtaruk_manuscript/figures/Figure3_temp_frost_snow_ice_soil.pdf", 
       Figure3, width = 30, height = 60, units = "cm")

# Figure 4. Active layer depth ----

# Models
# Active layer depth across years
act_layer_m_year <- MCMCglmm(value ~ I(year - 1985), random = ~ year + veg_type, data = act.layer, 
                             family = "gaussian", pr = TRUE, nitt = 100000, burnin = 20000)
summary(act_layer_m_year)

# Calculating model predictions
nyears <- 33
niter <- length(act_layer_m_year$Sol[,"(Intercept)"])
depth_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    depth_preds[i,j] <- act_layer_m_year$Sol[i,"(Intercept)"] + act_layer_m_year$Sol[i,"I(year - 1985)"]*j
  }
}

depth_preds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  depth_preds_df[i,] <- quantile(depth_preds[,i], c(0.025, 0.5, 0.975))
}

depth_preds_df <- cbind.data.frame(lower = depth_preds_df[,1], 
                                   mean = depth_preds_df[,2], 
                                   upper = depth_preds_df[,3], year = seq(1:33))

# Active layer depth across the growing season in Herschel veg type
# Calculate summary statistics
act.layer.sum <- act.layer.2017 %>% group_by(plot, veg.type, day) %>% 
  summarise(mean.depth = mean(depth), min.depth = min(depth), 
         max.depth = max(depth))

# Filter by site (Community composition plots vs Phenology ridge)
act.layer.2017.HE <- filter (act.layer.sum, veg.type == "HE")
act.layer.2017.KO <- filter (act.layer.sum, veg.type == "KO")

act_layer_m_HE <- MCMCglmm(mean.depth ~ I(day - 174), data = act.layer.2017.HE, 
                           family = "gaussian", pr = TRUE, nitt = 100000, burnin = 20000)
summary(act_layer_m_HE)

# Calculating model predictions
nyears <- 54
niter <- length(act_layer_m_HE$Sol[,"(Intercept)"])
depth_preds_HE <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    depth_preds_HE[i,j] <- act_layer_m_HE$Sol[i,"(Intercept)"] + act_layer_m_HE$Sol[i,"I(day - 174)"]*j
  }
}

depth_preds_HE_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  depth_preds_HE_df[i,] <- quantile(depth_preds_HE[,i], c(0.025, 0.5, 0.975))
}

depth_preds_HE_df <- cbind.data.frame(lower = depth_preds_HE_df[,1], 
                                      mean = depth_preds_HE_df[,2], 
                                      upper = depth_preds_HE_df[,3], day = seq(1:54))

# Active layer depth across the growing season in Komakuk veg type
act_layer_m_KO <- MCMCglmm(mean.depth ~ I(day - 174), data = act.layer.2017.KO, 
                           family = "gaussian", pr = TRUE, nitt = 100000, burnin = 20000)
summary(act_layer_m_KO)

nyears <- 54
niter <- length(act_layer_m_HE$Sol[,"(Intercept)"])
depth_preds_KO <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    depth_preds_KO[i,j] <- act_layer_m_KO$Sol[i,"(Intercept)"] + act_layer_m_KO$Sol[i,"I(day - 174)"]*j
  }
}

depth_preds_KO_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  depth_preds_KO_df[i,] <- quantile(depth_preds_KO[,i], c(0.025, 0.5, 0.975))
}

depth_preds_KO_df <- cbind.data.frame(lower = depth_preds_KO_df[,1], 
                                      mean = depth_preds_KO_df[,2], 
                                      upper = depth_preds_KO_df[,3], day = seq(1:54))


# Active layer depth across the years
act.layer.y$veg.type <- factor(act.layer.y$veg.type, levels = c("HE", "KO", "HE and/or KO", "HE and KO"))
act.layer.old <- filter(act.layer.y, year == "1985")

(act.layer.year <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbar(data = act.layer, aes(x = year, y = value, ymin = min.depth, 
                                        ymax = max.depth, colour = veg_type), size = 0.5, width = 0.5) +
    geom_point(data = act.layer, aes(x = year, y = value, colour = veg_type), size = 6, alpha = 0.8) +
    geom_ribbon(data = depth_preds_df, aes(x = year + 1985, ymin = lower, ymax = upper), 
                fill = "gray65", alpha = 0.2) +
    geom_line(data = depth_preds_df, aes(x = year + 1985, y = mean), colour = "gray65") +
    scale_color_manual(values = c("#ffa544", "#2b299b", "gray65", "black"), name = "", 
                       labels = c("Her.", "Kom.", "Her./Kom.")) +
    scale_fill_manual(values = c("#ffa544","#2b299b","gray65", "black")) +
    scale_x_continuous(breaks = c(1987, 2003, 2007, 2017)) +
    scale_y_reverse(lim=c(90,0), breaks = c(0, 20, 40, 60, 80)) +
    theme_QHI() + 
    labs(y = "Active layer depth (cm)\n", x = "\nYear", title = "(a) Across years\n") +
    theme(legend.position = c(0.9, 0.85), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5),
          axis.text = element_text(size = 24),
          axis.title = element_text(size = 30),
          plot.title = element_text(size = 30)) +
    guides(shape = FALSE, colour = guide_legend(override.aes = list(linetype = 0))))

act.layer.HE <- na.omit(act.layer.2017.HE)
act.layer.KO <- na.omit(act.layer.2017.KO)

# Active layer depth near community composition plots
(act.layer.season <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbar(data = act.layer.HE, aes(x = day, ymin = min.depth, ymax = max.depth),
               colour = "#ffa544", size = 0.5, width = 0.8) +
    geom_point(data = act.layer.HE, aes(x = day, y = mean.depth), 
               colour = "#ffa544", size = 6, alpha = 0.8) +
    geom_errorbar(data = act.layer.KO, aes(x = day, ymin = min.depth, ymax = max.depth),
               colour = "#2b299b", size = 0.5, width = 0.8) +
    geom_point(data = act.layer.KO, aes(x = day, y = mean.depth),
               colour = "#2b299b", size = 6, alpha = 0.8) +
    geom_ribbon(data = depth_preds_HE_df, aes(x = day + 174, ymin = lower, ymax = upper), 
                fill = "#ffa544", alpha = 0.2) +
    geom_line(data = depth_preds_HE_df, aes(x = day + 174, y = mean), colour = "#ffa544") +
    geom_ribbon(data = depth_preds_KO_df, aes(x = day + 174, ymin = lower, ymax = upper), 
                fill = "#2b299b", alpha = 0.2) +
    geom_line(data = depth_preds_KO_df, aes(x = day + 174, y = mean), colour = "#2b299b") +
    scale_color_manual(values = c("#ffa544", "#2b299b"), name="", labels = c("Her.", "Kom.")) +
    scale_fill_manual(values = c("#ffa544","#2b299b")) +
    scale_y_reverse(lim=c(90,0), breaks = c(0, 20, 40, 60, 80)) +
    theme_QHI() + 
    labs(y = "Mean active layer depth (cm)\n", x = "\nDay of year", title = "(b) Across the growing season\n") +
    theme(legend.position = c(0.1, 0.9), 
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5),
          axis.text = element_text(size = 24),
          axis.title = element_text(size = 30),
          plot.title = element_text(size = 30)) +
    guides(colour = guide_legend(override.aes = list(linetype = 0))))

# Arranging in a panel and saving the file
Figure4 <- grid.arrange(act.layer.year, act.layer.season, ncol = 2)

ggsave("Qikiqtaruk_manuscript/figures/Figure4_act_layer.pdf", 
       Figure4, width = 60, height = 30, units = "cm")

# Figure 6. Phenology change ----
# ** Phenology models ----

# Change column name Spp to SPP to match with script (column name was changes in 2 Nov. 2017 data update)
colnames(qiphen)[which(colnames(qiphen)=="Spp")] <- "SPP"

# P2 Model
# Keep only observations for which our variable of interest ("P2") was observed
qi <- qiphen[!is.na(qiphen$P2),]

# In order for the interval censoring analysis to work, we must have two columns of "response" data - one contains the "lower bound" (i.e., the last day on which an event was observed to have NOT yet occurred) and an "upper bound" (i.e., the first day on which an event was observed to have occurred). So we know that the event of interest happened somewhere between the upper bound and the lower bound.
head(qi[,c("P2_before","P2")])

# Make sure that P2 is always a larger number than P2_before, or else JAGS will be unhappy (NA's not accepted)
which(qi$P2_before>qi$P2)
qi[is.na(qi$P2_before),]

# Can set P2_before values that don't meet this criteria to 10 days before the "P2" date
# qi$P2_before[qi$P2_before > qi$P2 | is.na(qi$P2_before)] <- qi$P2[qi$P2_before > qi$P2 | is.na(qi$P2_before)] - 10

# Subset to SALARC
qi <- subset(qi, SPP == "SALARC")

# PREPARE DATA FOR JAGS IMPORT

qi <- qi[order(qi$Year, qi$Plot.ID),]
qi$PlotNum <- as.numeric(factor(qi$Plot.ID, levels=unique(qi$Plot.ID)))

head(qi)

jags.dat <- list(
  lim = cbind(qi$P2_before,qi$P2), # ANNE REMOVED -1 (subtracted from before date)
  y1s = ifelse(qi$P2_before < qi$P2,1,0), #column of 1's indicates iterval censoring
  year = qi$Year - min(qi$Year) + 1, #year must start with 1
  nyear = length(unique(qi$Year)),
  plot = qi$PlotNum,
  n = nrow(qi),
  nplot = length(unique(qi$PlotNum)),
  xhat = c(sort(unique(qi$Year - min(qi$Year) + 1))),
  nxhat = length(unique(qi$Year - min(qi$Year) + 1))
)

str(jags.dat) #make sure things that should be the same length are the same length and everything is numeric

write("
      model{
      # PRIORS
      a ~ dunif(0, 365)
      b ~ dnorm(0, 0.001)
      for (i in 1:nplot){
      aplot[i] ~ dnorm(0, tau.plot)
      }
      sigma.plot ~ dunif(0, 100)
      tau.plot <- 1/(sigma.plot*sigma.plot)
      sigma.year ~ dunif(0, 100)
      tau.year <- 1/(sigma.year*sigma.year)
      sigma ~ dunif(0, 100)
      tau <- 1/(sigma*sigma)
      #LIKELIHOOD
      for (j in 1:nyear){
      a.year[j]~dnorm(muyear[j], tau.year)
      muyear[j] <- a+b*j
      }
      for (i in 1:n){
      y1s[i]~dinterval(t[i],lim[i,])
      t[i]~dnorm(pred[i],tau)
      pred[i] <- aplot[plot[i]]+a.year[year[i]]
      }
      # DERIVED QUANTITIES
      for (j in 1:nxhat){
      phat[j] <- a + b*xhat[j]
      }
      }
      ","Qikiqtaruk_manuscript/models/int_cens_P2.jags")

# INITIAL VALUES ("t" is mandatory)
inits <- function() list(aplot = rnorm(jags.dat$nplot,0,2), sigma.plot = runif(1,0,1), sigma = runif(1,0,1), t = as.vector(apply(jags.dat$lim, 1, mean)))

# PARAMETERS TO MONITOR
params <- c("a", "b", "phat","aplot", "sigma", "sigma.plot", "sigma.year") ## ADD PARAMETERS TO MONITOR

modoutP2 <- jags(jags.dat, inits, params, model.file="Qikiqtaruk_manuscript/models/int_cens_P2.jags", n.chains=3, n.iter=20000, n.burnin=10000, n.thin=5, DIC=FALSE, working.directory=NULL, progress.bar = "text")

# plot(modoutP2) # check convergence

# EXTRACT OUTPUT
coefs <- as.data.frame(modoutP2$BUGSoutput$summary[,c('mean', 'sd', '2.5%', '97.5%', 'Rhat', 'n.eff')])
head(coefs,10)
hist(coefs$Rhat) # all should be below 1.1
coefs$Type <- as.vector(sapply(strsplit(rownames(coefs),"[[]",fixed=FALSE), "[", 1))

# GRAPH RESULTS
coefs.predP2 <- coefs[coefs$Type=="phat",]

# Add species and year information
coefs.predP2$SPP <- qi$SPP[1]
coefs.predP2$YearC <- c(1:17)
coefs.predP2$Year <- coefs.predP2$YearC + min(qi$Year) - 1

# Save coefs and subset of phenology data for later use
qiP2 <- qi
coefsP2 <- coefs.predP2
coefsP2.all <- coefs

# # 5% and 95% quantiles
# P2quants <- rbind(round(quantile(modoutP2$BUGSoutput$sims.list$a, probs = c(0.05, 0.95)), 2), round(quantile(modoutP2$BUGSoutput$sims.list$b, probs = c(0.05, 0.95)), 2), round(quantile(modoutP2$BUGSoutput$sims.list$sigma.year, probs = c(0.05, 0.95)), 2), round(quantile(modoutP2$BUGSoutput$sims.list$sigma.plot, probs = c(0.05, 0.95)), 2), round(quantile(modoutP2$BUGSoutput$sims.list$sigma, probs = c(0.05, 0.95)), 2))
# rownames(P2quants) <- c("a", "b", "sigma.year", "sigma.plot", "sigma")
# P2quants

# P3 Model
# Keep only observations for which our variable of interest ("P3") was observed
qi <- qiphen[!is.na(qiphen$P3),]

# In order for the interval censoring analysis to work, we must have two columns of "response" data - one contains the "lower bound" (i.e., the last day on which an event was observed to have NOT yet occurred) and an "upper bound" (i.e., the first day on which an even was observed to have occurred). So we know that the event of interest happened somewhere between the upper bound and the lower bound.
head(qi[,c("P3_before","P3")])

# Make sure that P3 is always a larger number than P3_before, or else JAGS will be unhappy (NA's not accepted)
which(qi$P3_before>qi$P3)
qi[is.na(qi$P3_before),]

# Can set P3_before values that don't meet this criteria to 10 days before the "P3" date
# qi$P3_before[qi$P3_before > qi$P3 | is.na(qi$P3_before)] <- qi$P3[qi$P3_before > qi$P3 | is.na(qi$P3_before)] - 10

# PREPARE DATA FOR JAGS IMPORT
qi$PlotNum <- as.numeric(as.factor(qi$Plot.ID))
qi$SpeciesNum <- as.numeric(as.factor(qi$SPP))
qi <- qi[order(qi$SPP, qi$Year, qi$PlotNum),]

head(qi)

jags.dat <- list(
  lim = cbind(qi$P3_before,qi$P3),
  y1s = ifelse(qi$P3_before<qi$P3,1,0), #column of 1's indicates iterval censoring
  year = qi$Year - min(qi$Year) + 1, #year must start with 1 (this becomes necessary in the more complex models)
  nyear = length(unique(qi$Year)),
  plot = qi$PlotNum,
  species = qi$SpeciesNum,
  n = nrow(qi),
  nspp = length(unique(qi$SpeciesNum)),
  nplot = length(unique(qi$PlotNum)),
  xhat = matrix(c(sort(unique(qi$Year - min(qi$Year) + 1))), nrow = 3,ncol = length(unique(qi$Year - min(qi$Year) + 1)),byrow = T),
  nxhat = length(unique(qi$Year - min(qi$Year) + 1))
)

str(jags.dat) #make sure things that should be the same length are the same length and everything is numeric

write("
      model{
      # PRIORS
      for (i in 1:nspp){
      aspp[i] ~ dunif(0, 365)
      bspp[i] ~ dnorm(0, 0.001)
      }
      for (i in 1:nplot){
      aplot[i] ~ dnorm(0, tau.plot)
      }
      sigma.plot ~ dunif(0, 100)
      tau.plot <- 1/(sigma.plot*sigma.plot)
      sigma.year ~ dunif(0, 100)
      tau.year <- 1/(sigma.year*sigma.year)
      sigma ~ dunif(0, 100)
      tau <- 1/(sigma*sigma)
      #LIKELIHOOD
      for (j in 1:nspp){
      for (k in 1:nyear){
      a.sppyear[j,k]~dnorm(musppyear[j,k], tau.year)
      musppyear[j,k] <- aspp[j]+bspp[j]*k
      }
      }
      for (i in 1:n){
      y1s[i]~dinterval(t[i],lim[i,])
      t[i]~dnorm(pred[i],tau)
      pred[i] <- aplot[plot[i]]+a.sppyear[species[i],year[i]]
      }
      # DERIVED QUANTITIES
      for (i in 1:nspp){
      for (j in 1:nxhat){
      phat[i,j] <- aspp[i] + bspp[i]*xhat[i,j]
      }
      }
      diff1to2 <- aspp[1] - aspp[2]
      diff1to3 <- aspp[1] - aspp[3]
      diff2to3 <- aspp[2] - aspp[3]
      }
      ","Qikiqtaruk_manuscript/models/int_cens_P3.jags")

# INITIAL VALUES ("t" is mandatory)
inits <- function() list(aplot = rnorm(jags.dat$nplot,0,2), sigma.plot = runif(1,0,1), sigma = runif(1,0,1), t = as.vector(apply(jags.dat$lim, 1, mean)))

# PARAMETERS TO MONITOR
params <- c("aspp", "bspp", "phat","aplot", "sigma", "sigma.plot", "sigma.year", "diff1to2", "diff1to3", "diff2to3") ## ADD PARAMETERS TO MONITOR

modoutP3 <- jags(jags.dat,inits, params, model.file="Qikiqtaruk_manuscript/models/int_cens_P3.jags", n.chains=3, n.iter=20000, n.burnin=10000, n.thin=5, DIC=FALSE, working.directory=NULL, progress.bar = "text")

# plot(modoutP3) #check convergence, etc.

# EXTRACT OUTPUT
coefs <- as.data.frame(modoutP3$BUGSoutput$summary[,c('mean', 'sd', '2.5%', '97.5%', 'Rhat', 'n.eff')])
head(coefs,10)
hist(coefs$Rhat)
coefs$Type <- as.vector(sapply(strsplit(rownames(coefs),"[[]",fixed=FALSE), "[", 1))

# GRAPH RESULTS
coefs.predP3 <- coefs[coefs$Type=="phat",] 

# Add species and year information
coefs.predP3$SpeciesNum <- rep(c(1:3), times = length(unique(jags.dat$year)))
coefs.predP3$SPP <- qi$SPP[match(coefs.predP3$SpeciesNum, qi$SpeciesNum)]
coefs.predP3$YearC <- rep(c(1:17), each = jags.dat$nspp)
coefs.predP3$Year <- coefs.predP3$YearC + min(qi$Year) - 1
# remove prediction for ERIVAG 2001 as no data is available for this year
coefs.predP3 <- coefs.predP3[!(coefs.predP3$Year == 2001 & coefs.predP3$SPP == "ERIVAG"),]

# Save coefs and subset of phenology data for later use
qiP3 <- qi
coefsP3 <- coefs.predP3
coefsP3.all <- coefs

# 5% and 95% quantiles
# P3quants <- rbind(round(quantile(modoutP3$BUGSoutput$sims.list$aspp[,1], probs = c(0.05, 0.95)), 2), round(quantile(modoutP3$BUGSoutput$sims.list$bspp[,1], probs = c(0.05, 0.95)), 2), round(quantile(modoutP3$BUGSoutput$sims.list$aspp[,2], probs = c(0.05, 0.95)), 2), round(quantile(modoutP3$BUGSoutput$sims.list$bspp[,2], probs = c(0.05, 0.95)), 2), round(quantile(modoutP3$BUGSoutput$sims.list$aspp[,3], probs = c(0.05, 0.95)), 2), round(quantile(modoutP3$BUGSoutput$sims.list$bspp[,3], probs = c(0.05, 0.95)), 2), round(quantile(modoutP3$BUGSoutput$sims.list$sigma.year, probs = c(0.05, 0.95)), 2), round(quantile(modoutP3$BUGSoutput$sims.list$sigma.plot, probs = c(0.05, 0.95)), 2), round(quantile(modoutP3$BUGSoutput$sims.list$sigma, probs = c(0.05, 0.95)), 2))
# rownames(P3quants) <- c("a1", "b1", "a2", "b2", "a3", "b3", "sigma.year", "sigma.plot", "sigma")
# P3quants

# P5 Model
# Keep only observations for which our variable of interest ("P5") was observed
qi <- qiphen[!is.na(qiphen$P5),]

# In order for the interval censoring analysis to work, we must have two columns of "response" data - one contains the "lower bound" (i.e., the last day on which an event was observed to have NOT yet occurred) and an "upper bound" (i.e., the first day on which an even was observed to have occurred). So we know that the event of interest happened somewhere between the upper bound and the lower bound.
head(qi[,c("P5_before","P5")])

# Make sure that P5 is always a larger number than P5_before, or else JAGS will be unhappy (NA's not accepted)
which(qi$P5_before>qi$P5)
qi[is.na(qi$P5_before),]

# Can set P5_before values that don't meet this criteria to 10 days before the "P5" date
# qi$P5_before[qi$P5_before > qi$P5 | is.na(qi$P5_before)] <- qi$P5[qi$P5_before > qi$P5 | is.na(qi$P5_before)] - 10

# Subset to SALARC
qi <- subset(qi, SPP == "SALARC")

# Thereis only one data point for 2001, remove that year from the dataset
qi <- qi[qi$Year > 2001,]

# PREPARE DATA FOR JAGS IMPORT
#qi$Plot.Num <- qi$Plot.ID
#qi$Plot.ID <- recode(qi$Plot.ID, F1 = "1", F2 = "2", F3 = "3", F4 = "4", F5 = "5", F6 = "6", F7 = "7", F8 = "8", F9 = "9", F10 = "10", M1 = "11", M2 = "12", M3 = "13", M4 = "14", M5 = "15", M6 = "16", M7 = "17", M8 = "18", M9 = "19", M10 = "20") #order doesn't matter here
#qi$PlotNum <- as.numeric(as.character(qi$Plot.ID))
qi$PlotNum <- as.numeric(as.factor(qi$Plot.ID))
qi <- qi[order(qi$Year, qi$PlotNum),]

jags.dat <- list(
  lim = cbind(qi$P5_before,qi$P5), #changed P5_after to P5 as this is actually the end of the window
  y1s = ifelse(qi$P5_before<qi$P5,1,0), #column of 1's indicates iterval censoring
  year = qi$Year - min(qi$Year) + 1, #year must start with 1 (this becomes necessary in the more complex models)
  nyear = length(unique(qi$Year)),
  plot = qi$PlotNum,
  n = nrow(qi),
  nplot = length(unique(qi$PlotNum)),
  xhat = c(sort(unique(qi$Year - min(qi$Year) + 1))),
  nxhat = length(unique(qi$Year - min(qi$Year) + 1))
)

str(jags.dat) #make sure things that should be the same length are the same length and everything is numeric

write("
      model{
      # PRIORS
      a ~ dunif(0, 365)
      b ~ dnorm(0, 0.001)
      for (i in 1:nplot){
      aplot[i] ~ dnorm(0, tau.plot)
      }
      sigma.plot ~ dunif(0, 100)
      tau.plot <- 1/(sigma.plot*sigma.plot)
      sigma.year ~ dunif(0, 100)
      tau.year <- 1/(sigma.year*sigma.year)
      sigma ~ dunif(0, 100)
      tau <- 1/(sigma*sigma)
      #LIKELIHOOD
      for (j in 1:nyear){
      a.year[j]~dnorm(muyear[j], tau.year)
      muyear[j] <- a+b*j
      }
      for (i in 1:n){
      y1s[i]~dinterval(t[i],lim[i,])
      t[i]~dnorm(pred[i],tau)
      pred[i] <- aplot[plot[i]]+a.year[year[i]]
      }
      # DERIVED QUANTITIES
      for (j in 1:nxhat){
      phat[j] <- a + b*xhat[j]
      }
      }
      ","Qikiqtaruk_manuscript/models/int_cens_P5.jags")

# INITIAL VALUES
inits <- function() list(aplot = rnorm(jags.dat$nplot,0,2), sigma.plot = runif(1,0,1), sigma = runif(1,0,1), t = as.vector(apply(jags.dat$lim, 1, mean)))

# PARAMETERS TO MONITOR
params <- c("a", "b", "phat","aplot", "sigma", "sigma.plot", "sigma.year") ## ADD PARAMETERS TO MONITOR

modoutP5 <- jags(jags.dat, inits, params, model.file="Qikiqtaruk_manuscript/models/int_cens_P5.jags", n.chains=3, n.iter=20000, n.burnin=10000, n.thin=5, DIC=FALSE, working.directory=NULL, progress.bar = "text")

# plot(modoutP5) #check convergence, etc.

# EXTRACT OUTPUT
coefs <- as.data.frame(modoutP5$BUGSoutput$summary[,c('mean', 'sd', '2.5%', '97.5%', 'Rhat', 'n.eff')])
head(coefs,10)
hist(coefs$Rhat)
coefs$Type <- as.vector(sapply(strsplit(rownames(coefs),"[[]",fixed=FALSE), "[", 1))

# GRAPH RESULTS
# Graph the trends over time for each species. Include both the raw data and the modeled predictions.
coefs.predP5 <- coefs[coefs$Type=="phat",] #change this name depending on what you called it

# Add species and year information
coefs.predP5$SPP <- qi$SPP[1]
coefs.predP5$YearC <- c(1:16) #change this if using all years!
coefs.predP5$Year <- coefs.predP5$YearC + min(qi$Year) - 1

# Save coefs and subset of phenology data for later use
qiP5 <- qi
coefsP5 <- coefs.predP5
coefsP5.all <- coefs

# 5% and 95% quantiles
# P5quants <- rbind(round(quantile(modoutP5$BUGSoutput$sims.list$a, probs = c(0.05, 0.95)), 2), round(quantile(modoutP5$BUGSoutput$sims.list$b, probs = c(0.05, 0.95)), 2), round(quantile(modoutP5$BUGSoutput$sims.list$sigma.year, probs = c(0.05, 0.95)), 2), round(quantile(modoutP5$BUGSoutput$sims.list$sigma.plot, probs = c(0.05, 0.95)), 2), round(quantile(modoutP5$BUGSoutput$sims.list$sigma, probs = c(0.05, 0.95)), 2))
# rownames(P5quants) <- c("a", "b", "sigma.year", "sigma.plot", "sigma")
# P5quants

# Growing Season Length

# DIFFERENCE between P2 and P5 for Salix arctica only
# length of time between green-up and senescence (Salix only)
# look at P6 instead? (date of last leaf turning yellow instead of first)
# NB: There is VIRTUALLY NO data for senescense in 2001 (1 data point) so the model below only runs from 2002 onwards
# Could also model this DIRECTLY - i.e. P5-P2 per plot per year - tried this but it doesn't work with censoring
# Need to add in right- and left- censoring

# Keep only observations for which our variables of interest ("P2" and "P5") were observed
qiphenP2 <- qiphen[!is.na(qiphen$P2) & qiphen$SPP=="SALARC",]
qiphenP5 <- qiphen[!is.na(qiphen$P5) & qiphen$SPP=="SALARC",]

# In order for the interval censoring analysis to work, we must have two columns of "response" data - one contains the "lower bound" (i.e., the last day on which an event was observed to have NOT yet occurred) and an "upper bound" (i.e., the first day on which an even was observed to have occurred). So we know that the event of interest happened somewhere between the upper bound and the lower bound.
head(qiphenP2[,c("P2_before","P2")])

# Make sure that P2_after is always a larger number than P2_before, or else JAGS will be unhappy (NA's not accepted)
which(qiphenP2$P2_before>qiphenP2$P2)
which(qiphenP5$P5_before>qiphenP5$P5)
qiphenP2[is.na(qiphenP2$P2_before) & !is.na(qiphenP2$P2),]
qiphenP5[is.na(qiphenP5$P5_before) & !is.na(qiphenP5$P5),]

qiphenP2gs <- qiphenP2 
qiphenP5gs <- qiphenP5[qiphenP5$Year > 2001,]

# PREPARE DATA FOR JAGS IMPORT

qiphenP2gs <- qiphenP2gs[order(qiphenP2gs$Year, qiphenP2gs$Plot.ID),]
qiphenP2gs$PlotNum <- as.numeric(factor(qiphenP2gs$Plot.ID, levels = unique(qiphenP2gs$Plot.ID)))

qiphenP5gs <- qiphenP5gs[order(qiphenP5gs$Year, qiphenP5gs$Plot.ID),]
qiphenP5gs$PlotNum <- as.numeric(factor(qiphenP5gs$Plot.ID, levels = unique(qiphenP5gs$Plot.ID)))

jags.dat <- list(
  limP2=cbind(qiphenP2gs$P2_before,qiphenP2gs$P2),
  limP5=cbind(qiphenP5gs$P5_before,qiphenP5gs$P5),
  y1sP2=rep(1,length(qiphenP2gs$P2_before)), #column of 1's indicates iterval censoring
  y1sP5=rep(1,length(qiphenP5gs$P5_before)),
  yearP2=qiphenP2gs$Year-min(qiphenP2gs$Year)+1, #year must start with 1
  yearP5=qiphenP5gs$Year-min(qiphenP5gs$Year)+1,
  nyearP2=length(unique(qiphenP2gs$Year)),
  nyearP5=length(unique(qiphenP5gs$Year)),
  plotP2=qiphenP2gs$PlotNum,
  plotP5=qiphenP5gs$PlotNum,
  nP2=nrow(qiphenP2gs),
  nP5=nrow(qiphenP5gs),
  nplotP2=length(unique(qiphenP2gs$PlotNum)),
  nplotP5=length(unique(qiphenP5gs$PlotNum)),
  xhatP2=c(sort(unique(qiphenP2gs$Year-min(qiphenP2gs$Year)+1))),
  xhatP5=c(sort(unique(qiphenP5gs$Year-min(qiphenP5gs$Year)+1))),
  nxhatP2=length(unique(qiphenP2gs$Year-min(qiphenP2gs$Year)+1)),
  nxhatP5=length(unique(qiphenP5gs$Year-min(qiphenP5gs$Year)+1))
)

str(jags.dat) #make sure things that should be the same length are the same length and everything is numeric

write("
      model{

      #PRIORS P2
      aP2~dunif(0,365)
      bP2~dnorm(0,0.001)
      #Priors for plot random effect
      for (i in 1:nplotP2){
      aplotP2[i]~dnorm(0,tau.plotP2)
      }
      sigma.plotP2~dunif(0,100)
      tau.plotP2 <- 1/(sigma.plotP2*sigma.plotP2)
      sigmaP2~dunif(0,100)
      tauP2 <- 1/(sigmaP2*sigmaP2)
      sigma.yearP2~dunif(0,100)
      tau.yearP2 <- 1/(sigma.yearP2*sigma.yearP2)

      #PRIORS P5
      aP5~dunif(0,365)
      bP5~dnorm(0,0.001)
      #Priors for plot random effect
      for (i in 1:nplotP5){
      aplotP5[i]~dnorm(0,tau.plotP5)
      }
      sigma.plotP5~dunif(0,100)
      tau.plotP5 <- 1/(sigma.plotP5*sigma.plotP5)
      sigmaP5~dunif(0,100)
      tauP5 <- 1/(sigmaP5*sigmaP5)
      sigma.yearP5~dunif(0,100)
      tau.yearP5 <- 1/(sigma.yearP5*sigma.yearP5)

      #LIKELIHOOD P2
      for (k in 1:nyearP2){
      #this adds in the random variation
      a.yearP2[k]~dnorm(muyearP2[k], tau.yearP2)
      #this is the predicted mean effect per year
      muyearP2[k] <- aP2+bP2*k
      }
      for (i in 1:nP2){
      y1sP2[i]~dinterval(tP2[i],limP2[i,])
      #y1s is a column of 1's (because t[i] always falls between lim[i,1] and lim[i,2])
      tP2[i]~dnorm(predP2[i],tauP2)
      predP2[i] <- aplotP2[plotP2[i]]+a.yearP2[yearP2[i]]
      }

      #LIKELIHOOD P5
      for (k in 1:nyearP5){
      #this adds in the random variation
      a.yearP5[k]~dnorm(muyearP5[k], tau.yearP5)
      #this is the predicted mean effect per year
      muyearP5[k] <- aP5+bP5*k
      }
      for (i in 1:nP5){
      y1sP5[i]~dinterval(tP5[i],limP5[i,])
      #y1s is a column of 1's (because t[i] always falls between lim[i,1] and lim[i,2])
      tP5[i]~dnorm(predP5[i],tauP5)
      predP5[i] <- aplotP5[plotP5[i]]+a.yearP5[yearP5[i]]
      }
      # DERIVED QUANTITIES
      for (j in 1:nxhatP2){
      yhatP2[j] <- aP2+bP2*xhatP2[j]
      }
      for (j in 1:nxhatP5){
      yhatP5[j] <- aP5+bP5*xhatP5[j]
      }
      slopeDiff <- bP5-bP2
      for (j in 1:nxhatP5){
      yhatDiff[j] <- yhatP5[j]-yhatP2[j+1] #add 1 to yhatP2 because it starts one year earlier
      }
      }
      ","Qikiqtaruk_manuscript/models/growing_season_pheno.jags")

# INITIAL VALUES

inits <- function() list(aplotP2=rnorm(jags.dat$nplotP2,0,2), tP2=as.vector(apply(jags.dat$limP2,1,mean)), sigma.plotP2=runif(1,0,1), sigma.yearP2=runif(1,0,1), sigmaP2=runif(1,0,1), aplotP5=rnorm(jags.dat$nplotP2,0,2), tP5=as.vector(apply(jags.dat$limP5,1,mean)), sigma.plotP5=runif(1,0,1), sigma.yearP5=runif(1,0,1), sigmaP5=runif(1,0,1))

# PARAMETERS TO MONITOR

params <- c("aplotP2","aplotP5","sigma.plotP2","sigma.plotP5","aP2","aP5","sigmaP2","sigmaP5","bP2","bP5","yhatP2","yhatP5","sigma.yearP2","sigma.yearP5","slopeDiff","yhatDiff")

modoutP2P5 <- jags(jags.dat,inits, params, model.file="Qikiqtaruk_manuscript/models/growing_season_pheno.jags", n.chains=3,n.iter=30000,n.burnin=15000, n.thin=2, DIC=FALSE, working.directory=NULL, progress.bar = "text")

# plot(modoutP2P5) #check convergence, etc.

# EXTRACT OUTPUT
coefs <- as.data.frame(modoutP2P5$BUGSoutput$summary)
head(coefs,10)
hist(coefs$Rhat)
coefs$Type <- as.vector(sapply(strsplit(rownames(coefs),"[[]",fixed=FALSE), "[", 1))

# GRAPH DIFFERENCE BETWEEN P2 and P5
coefs.diff <- coefs[coefs$Type=="yhatDiff",]
coefs.diff$YearC <- 1:16
coefs.diff$Year <- coefs.diff$YearC+min(qiphenP2gs$Year)-1+1 # +1 b/c started at year 2

# Save coefs for later use
coefsP2P5.all <- coefs
coefsP2P5 <- coefs.diff

# Get probability that the change in growing season length is greater than 0
slopeDiff <- modoutP2P5$BUGSoutput$sims.list$slopeDiff
pnorm(0, mean(slopeDiff), sd(slopeDiff), lower.tail=F)
#probability that a normally distributed random number with this mean and sd will be greater than 0

# # 5% and 95% quantiles
# P2P5quants <- round(quantile(modoutP2P5$BUGSoutput$sims.list$slopeDiff, probs = c(0.05, 0.95)), 2)
# P2P5quants

# Model coefficients summaries
coefsP2
coefsP3
coefsP5
coefsP2P5

# ** Phenology Figure ----

Figure6 <- grid.arrange(
  ggplot() +
    geom_point(data = qiP2,
               aes(x = Year, y = as.vector(apply(cbind(qiP2$P2_before, qiP2$P2), 1, mean)), colour = SPP),
               alpha = 0.5, size = 3) +
    geom_ribbon(data = coefs.predP2,
                aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, fill = SPP),
                alpha = 0.5, show.legend=FALSE) +
    geom_line(data = coefs.predP2,
              aes(x = Year, y = mean, colour = SPP), size = 1, show.legend=FALSE) +
    theme_QHI() +
    scale_colour_manual(values = c("#66c2a5", "#fc8d62", "#6157d1"), name = "SPECIES") +
    scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#6157d1"), name = "SPECIES") +
    scale_x_continuous(breaks = c(2001, 2006, 2011, 2017)) +
    scale_y_continuous(breaks = c(140, 160, 180, 200)) +
    coord_cartesian(ylim = c(130, 210)) +
    labs(x = " ", y = "Day of Year\n", title = expression(paste("(a) ", italic("Salix arctica"), " leaf out (P2)")), subtitle = "")  +
    theme(legend.position="none", axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5)),
  
  ggplot() +
    geom_point(data = qiP3,
               aes(x = Year, y = as.vector(apply(cbind(qiP3$P3before, qiP3$P3), 1, mean)), colour = SPP),
               alpha = 0.5, size = 3) +
    geom_ribbon(data = coefs.predP3, aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, fill = SPP),
                alpha = 0.5) +
    geom_line(data = coefs.predP3,
              aes(x = Year,y = mean, colour = SPP)) +
    theme(legend.position="none") +
    theme_QHI() +
    scale_colour_manual(values = c("#6157d1", "#fc8d62", "#66c2a5"), name = "", labels = c("D. integrifolia", "E. vaginatum", "S. arctica")) +
    scale_fill_manual(values = c("#6157d1", "#fc8d62", "#66c2a5"), name = "", labels = c("D. integrifolia", "E. vaginatum", "S. arctica")) +
    scale_x_continuous(breaks = c(2001, 2006, 2011, 2017))+
    scale_y_continuous(breaks=c(140, 160, 180, 200, 220))+
    coord_cartesian(ylim = c(130, 220)) +
    labs(x = " ", y = "Day of Year\n", title = "(b) Flowering (P3)\n", subtitle = ""),
  
  ggplot() +
    geom_point(data = qiphen[qiphen$SPP=="SALARC" & !is.na(qiphen$P5),], #anne changed this from qiP5 bc Isla wants 2001 in graph
               aes(x = Year, y = as.vector(apply(cbind(qiphen$P5before[qiphen$SPP=="SALARC" & !is.na(qiphen$P5)], qiphen$P5[qiphen$SPP=="SALARC" & !is.na(qiphen$P5)]), 1, mean)), colour = SPP),
               alpha = 0.5, size = 3) +
    geom_ribbon(data = coefs.predP5, aes(x = Year, ymin = `2.5%`, ymax = `97.5%`, fill = SPP),
                alpha = 0.5, show.legend = FALSE) +
    geom_line(data = coefs.predP5, aes(x = Year, y = mean, colour = SPP), show.legend = FALSE) +
    theme_QHI() +
    scale_colour_manual(values = c("#66c2a5", "#fc8d62", "#6157d1"), name = "SPECIES") +
    scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#6157d1"), name = "SPECIES") +
    scale_x_continuous(breaks=c(2001, 2006, 2011, 2017))+
    scale_y_continuous(breaks=c(190, 210, 230, 250))+
    coord_cartesian(ylim = c(180, 260)) +
    labs(x=" ", y = "Day of Year\n",  title = expression(paste("(c) \n", italic("Salix arctica\n"), " senescence (P5)\n")), subtitle = "") +
    theme(legend.position="none", axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5)),
  
  ggplot() +
    geom_point(data = qiphen[qiphen$SPP=="SALARC",],
               aes(x = Year, y = as.vector(apply(cbind(qiphen$P5_before[qiphen$SPP=="SALARC"],
                                                       qiphen$P5[qiphen$SPP=="SALARC"]), 1, mean) -
                                             apply(cbind(qiphen$P2_before[qiphen$SPP=="SALARC"],
                                                         qiphen$P2[qiphen$SPP=="SALARC"]), 1, mean))),
               position = position_dodge(width=0.2),
               alpha = 0.5, colour="#66c2a5", size = 3) +
    geom_ribbon(data = coefs.diff, aes(x = Year, ymin = `2.5%`, ymax = `97.5%`),
                alpha = 0.5, fill = "#66c2a5", show.legend = FALSE) +
    geom_line(data = coefs.diff, aes(x = Year, y = mean), colour = "#66c2a5", show.legend = FALSE) +
    theme(legend.position="none") +
    theme_QHI() +
    scale_x_continuous(breaks =  c(2001, 2006, 2011, 2017)) +
    scale_y_continuous(breaks=c(0, 25, 50, 75, 100)) +
    coord_cartesian(ylim = c(0, 110)) +
    labs(x = " ", y = "Days\n", title = expression(paste("(d) \n", italic("Salix arctica\n")," growing season (P5-P2)\n")), subtitle = "")
)

# Arranging in a panel and saving the file
ggsave("Qikiqtaruk_manuscript/figures/Figure6_phenology.pdf", 
       Figure6, width = 33, height = 33, units = "cm")

### Calculate mean and SD for interval censored data
# We will do so by calcualting the mean between the observed data and the data prior
# This is certainly not the most elegant method and a lot of informaiton will be lost,
# but for now it will have to do.

# Calculate mean between observed date and last visit to the plot
qiphen <- qiphen %>% mutate(P1_mean = (P1 + P1_before)/2,
                            P2_mean = (P2 + P2_before)/2,
                            P3_mean = (P3 + P3_before)/2,
                            P4_mean = (P4 + P4_before)/2,
                            P5_mean = (P5 + P5_before)/2,
                            P6_mean = (P6 + P6_before)/2,
                            P7_mean = (P7 + P7_before)/2)

# Calculate sd of difference in mean values of P5 and P2 (SALARC growing season length)
qiphen <- qiphen %>% mutate(P5_P2_mean = P5_mean-P2_mean)

# Calculate SD of mean values of intervals for each spp
qiphen_mean_sds <- qiphen %>% group_by(SPP) %>% summarise(P1_mean_sd = sd(P1_mean, na.rm =T),
                                                          P2_mean_sd = sd(P2_mean, na.rm =T),
                                                          P3_mean_sd = sd(P3_mean, na.rm =T),
                                                          P4_mean_sd = sd(P4_mean, na.rm =T),
                                                          P5_mean_sd = sd(P5_mean, na.rm =T),
                                                          P6_mean_sd = sd(P6_mean, na.rm =T),
                                                          P7_mean_sd = sd(P7_mean, na.rm =T),
                                                          P5_P2_mean_sd = sd(P5_P2_mean, na.rm = T))

# Output coefficients and effect sizes

table.P2 <- rbind(cbind.data.frame(Variable = "Intercept", coefsP2.all[coefsP2.all$Type == "a", c("mean","2.5%","97.5%","n.eff")]),
                     cbind(Variable = "Year", coefsP2.all[coefsP2.all$Type == "b", c("mean","2.5%","97.5%","n.eff")]),
                     cbind(Variable = "Sigma-Year", coefsP2.all[coefsP2.all$Type == "sigma.year", c("mean","2.5%","97.5%","n.eff")]),
                     cbind(Variable = "Sigma-Plot", coefsP2.all[coefsP2.all$Type == "sigma.plot", c("mean","2.5%","97.5%","n.eff")]),
                     cbind(Variable = "Sigma-Resid", coefsP2.all[coefsP2.all$Type == "sigma", c("mean","2.5%","97.5%","n.eff")]))
table.P2$Model <- "P2"

table.P3 <- rbind(cbind.data.frame(Variable = "Intercept", coefsP3.all[coefsP3.all$Type == "aspp", c("mean","2.5%","97.5%","n.eff")][1,]),
                  cbind(Variable = "Year", coefsP3.all[coefsP3.all$Type == "bspp", c("mean","2.5%","97.5%","n.eff")][1,]),
                  cbind.data.frame(Variable = "Intercept", coefsP3.all[coefsP3.all$Type == "aspp", c("mean","2.5%","97.5%","n.eff")][2,]),
                  cbind(Variable = "Year", coefsP3.all[coefsP3.all$Type == "bspp", c("mean","2.5%","97.5%","n.eff")][2,]),
                  cbind.data.frame(Variable = "Intercept", coefsP3.all[coefsP3.all$Type == "aspp", c("mean","2.5%","97.5%","n.eff")][3,]),
                  cbind(Variable = "Year", coefsP3.all[coefsP3.all$Type == "bspp", c("mean","2.5%","97.5%","n.eff")][3,]),
                  cbind(Variable = "Sigma-Year", coefsP3.all[coefsP3.all$Type == "sigma.year", c("mean","2.5%","97.5%","n.eff")]),
                  cbind(Variable = "Sigma-Plot", coefsP3.all[coefsP3.all$Type == "sigma.plot", c("mean","2.5%","97.5%","n.eff")]),
                  cbind(Variable = "Sigma-Resid", coefsP3.all[coefsP3.all$Type == "sigma", c("mean","2.5%","97.5%","n.eff")]))
table.P3$Model <- "P3"

table.P5 <- rbind(cbind.data.frame(Variable = "Intercept", coefsP5.all[coefsP5.all$Type == "a", c("mean","2.5%","97.5%","n.eff")]),
                  cbind(Variable = "Year", coefsP5.all[coefsP5.all$Type == "b", c("mean","2.5%","97.5%","n.eff")]),
                  cbind(Variable = "Sigma-Year", coefsP5.all[coefsP5.all$Type == "sigma.year", c("mean","2.5%","97.5%","n.eff")]),
                  cbind(Variable = "Sigma-Plot", coefsP5.all[coefsP5.all$Type == "sigma.plot", c("mean","2.5%","97.5%","n.eff")]),
                  cbind(Variable = "Sigma-Resid", coefsP5.all[coefsP5.all$Type == "sigma", c("mean","2.5%","97.5%","n.eff")]))
table.P5$Model <- "P5"

table.P2P5 <- rbind(cbind.data.frame(Variable = "Year", coefsP2P5.all[coefsP2P5.all$Type == "slopeDiff", c("mean","2.5%","97.5%","n.eff")]))
table.P2P5$Model <- "P2P5"

table.pheno.all <- rbind(table.P2, table.P3, table.P5, table.P2P5)

write.csv(table.pheno.all, "Qikiqtaruk_manuscript/model_outputs/Phenology_JAGS_coefficients_all_models.csv")

# standardized effect sizes

pheno.es <- table.pheno.all[table.pheno.all$Variable=="Year",]
pheno.es$Data_SD[pheno.es$Model=="P2"] <- as.numeric(qiphen_mean_sds[3,"P2_mean_sd"])
pheno.es$Data_SD[pheno.es$Model=="P3"][1] <- as.numeric(qiphen_mean_sds[1,"P3_mean_sd"])
pheno.es$Data_SD[pheno.es$Model=="P3"][2] <- as.numeric(qiphen_mean_sds[2,"P3_mean_sd"])
pheno.es$Data_SD[pheno.es$Model=="P3"][3] <- as.numeric(qiphen_mean_sds[3,"P3_mean_sd"])
pheno.es$Data_SD[pheno.es$Model=="P5"] <- as.numeric(qiphen_mean_sds[3,"P5_mean_sd"])
pheno.es$Data_SD[pheno.es$Model=="P2P5"] <- as.numeric(qiphen_mean_sds[3,"P5_P2_mean_sd"])

write.csv(pheno.es, file="Qikiqtaruk_manuscript/model_outputs/Phenology_JAGS_EffectSizes_all_models.csv")

# Figure 7. Canopy height and radial growth ----

# HE canopy height model
HE_canopy_m <- MCMCglmm(Mean.C.H ~ I(YEAR - 1998), random = ~ YEAR + PLOT, 
                        family = "gaussian", data = avg_heights[avg_heights$SUBSITE == "HE",], 
                        pr=TRUE, nitt = 100000, burnin = 20000, prior = prior2)
summary(HE_canopy_m)
# plot(HE_canopy_m$VCV)
# plot(HE_canopy_m$Sol)
autocorr(HE_canopy_m$VCV)

# Calculating model predictions
nyears <- 19
niter <- length(HE_canopy_m$Sol[,"(Intercept)"])
HE_canopy_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    HE_canopy_preds[i,j] <- HE_canopy_m$Sol[i,"(Intercept)"] + HE_canopy_m$Sol[i,"I(YEAR - 1998)"]*j
  }
}

HE_canopy_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  HE_canopy_preds_df[i,] <- quantile(HE_canopy_preds[,i], c(0.025, 0.5, 0.975))
}

HE_canopy_preds_df <- cbind.data.frame(lower = HE_canopy_preds_df[,1], 
                                       mean = HE_canopy_preds_df[,2], 
                                       upper = HE_canopy_preds_df[,3], year = seq(1:19))

# KO canopy height model
KO_canopy_m <- MCMCglmm(Mean.C.H ~ I(YEAR - 1998), random = ~ YEAR + PLOT, 
                        family = "gaussian", data = avg_heights[avg_heights$SUBSITE=="KO",], 
                        pr=TRUE, nitt = 100000, burnin = 20000, prior = prior2)
summary(KO_canopy_m)
# plot(KO_canopy_m$VCV)
# plot(KO_canopy_m$Sol)
autocorr(KO_canopy_m$VCV)

# Calculating model predictions
nyears <- 19
niter <- length(KO_canopy_m$Sol[,"(Intercept)"])
KO_canopy_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    KO_canopy_preds[i,j] <- KO_canopy_m$Sol[i,"(Intercept)"] + KO_canopy_m$Sol[i,"I(YEAR - 1998)"]*j
  }
}

KO_canopy_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  KO_canopy_preds_df[i,] <- quantile(KO_canopy_preds[,i], c(0.025, 0.5, 0.975))
}

KO_canopy_preds_df <- cbind.data.frame(lower = KO_canopy_preds_df[,1], 
                                       mean = KO_canopy_preds_df[,2], 
                                       upper = KO_canopy_preds_df[,3], year = seq(1:19))

# HE + KO canopy height graph
(canopy.height <- ggplot() +
    geom_point(data = avg_heights, aes(x = YEAR, y = Mean.C.H, colour = factor(SUBSITE), 
                                       fill = factor(SUBSITE)), alpha = 0.8, size = 4) +
    scale_color_manual(values = c("#ffa544", "#2b299b"), name = "", labels = c("Her.", "Kom.")) +
    scale_fill_manual(values = c("#ffa544","#2b299b")) +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2013, 2017)) +
    geom_line(data = HE_canopy_preds_df, aes(x = year + 1998, y = mean), colour = "#ffa544") +
    geom_ribbon(data = HE_canopy_preds_df, aes(x = year + 1998, ymin = lower, ymax = upper),
                fill = "#ffa544", alpha = 0.2) +
    geom_line(data = KO_canopy_preds_df, aes(x = year + 1998, y = mean), colour = "#2b299b") +
    geom_ribbon(data = KO_canopy_preds_df, aes(x = year + 1998, ymin = lower, ymax = upper),
                fill = "#2b299b", alpha = 0.2) +
    guides(fill = FALSE) +
    theme_QHI() +
    theme(legend.position = c(0.1, 0.95), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    labs(x = " ", y = "Mean canopy height (cm)\n", title = "(a) Mean canopy height\n"))

# Salix pulchra model
salpuls <- subset(abundance, SPP == "SALPUL"|SPP == "Salix pulchra")
salpuls <- subset(salpuls, !is.na(Height..cm.))
avg_salpuls <- ddply(salpuls, .(YEAR, SUBSITE, PLOT), summarise,
                     Mean.C.H = mean(Height..cm.))

Salix_canopy_m <- MCMCglmm(Height..cm. ~ I(YEAR - 1998), random = ~ YEAR + PLOT, 
                           family = "gaussian", data = salpuls, pr=TRUE, nitt = 100000, burnin = 20000, prior = prior2)
summary(Salix_canopy_m)
# plot(Salix_canopy_m$VCV)
# plot(Salix_canopy_m$Sol)
autocorr(Salix_canopy_m$VCV)

# Calculating model predictions
nyears <- 19
niter <- length(Salix_canopy_m$Sol[,"(Intercept)"])
Salix_canopy_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    Salix_canopy_preds[i,j] <- Salix_canopy_m$Sol[i,"(Intercept)"] + Salix_canopy_m$Sol[i,"I(YEAR - 1998)"]*j
  }
}

Salix_canopy_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  Salix_canopy_preds_df[i,] <- quantile(Salix_canopy_preds[,i], c(0.025, 0.5, 0.975))
}

Salix_canopy_preds_df <- cbind.data.frame(lower = Salix_canopy_preds_df[,1], 
                                          mean = Salix_canopy_preds_df[,2], 
                                          upper = Salix_canopy_preds_df[,3], 
                                          year = seq(1:19))

# Salix pulchra graph
(salix <- ggplot() +
    geom_point(data = avg_salpuls, aes(x = YEAR, y = Mean.C.H), 
               alpha = 0.8, colour = "#008c5f", size = 4) +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2013, 2017)) +
    geom_ribbon(data = Salix_canopy_preds_df, aes(x = year + 1998, 
                                                  ymin = lower, ymax = upper), 
                fill = "#008c5f", alpha = 0.2) +
    geom_line(data = Salix_canopy_preds_df, aes(x = year + 1998, y = mean), 
              colour = "#008c5f") +
    theme_QHI() +
    labs(x = "", y = "Height (cm)\n", 
         title = expression(paste("(b) ", italic("Salix pulchra"), " canopy height")),
         subtitle = ""))

# Dendro analysis
dendro.salric <- subset(dendro, Sp == "Salix richardsonii")
dendro.salpul <- subset(dendro, Sp == "Salix pulchra")
dendro.salarc <- subset(dendro, Sp == "Salix arctica")
dendro.salgla <- subset(dendro, Sp == "Salix glauca")

# Individual models for each species

# Salix richardsonii
salric_m <- MCMCglmm(rw ~ I(Year-1973), random = ~ Year + IndivUN, data = dendro.salric, 
                     family="gaussian", pr=TRUE, nitt = 100000, burnin = 20000, prior = prior2) 

summary(salric_m)
# plot(salric_m$VCV)
# plot(salric_m$Sol)
autocorr(salric_m$VCV)

nyears <- length(unique(dendro.salric$Year))
niter <- length(salric_m$Sol[,"(Intercept)"])
Salix_richardsonii_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    Salix_richardsonii_preds[i,j] <- salric_m$Sol[i,"(Intercept)"] + salric_m$Sol[i,"I(Year - 1973)"]*j
  }
}

Salix_richardsonii_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  Salix_richardsonii_preds_df[i,] <- quantile(Salix_richardsonii_preds[,i], c(0.025, 0.5, 0.975))
}

Salix_richardsonii_preds_df <- cbind.data.frame(lower = Salix_richardsonii_preds_df[,1], 
                                                mean = Salix_richardsonii_preds_df[,2], 
                                                upper = Salix_richardsonii_preds_df[,3], year = seq(1:nyears))

# Salix pulchra
salpul_m <- MCMCglmm(rw ~ I(Year-1973), random = ~ Year + IndivUN, data = dendro.salpul, 
                     family="gaussian", pr=TRUE, nitt = 100000, burnin = 20000, prior = prior2) 

summary(salpul_m)
# plot(salpul_m$VCV)
# plot(salpul_m$Sol)
autocorr(salpul_m$VCV)

nyears <- length(unique(dendro.salpul$Year))
niter <- length(salpul_m$Sol[,"(Intercept)"])
Salix_pulchra_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    Salix_pulchra_preds[i,j] <- salpul_m$Sol[i,"(Intercept)"] + salpul_m$Sol[i,"I(Year - 1973)"]*j
  }
}

Salix_pulchra_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  Salix_pulchra_preds_df[i,] <- quantile(Salix_pulchra_preds[,i], c(0.025, 0.5, 0.975))
}

Salix_pulchra_preds_df <- cbind.data.frame(lower = Salix_pulchra_preds_df[,1], 
                                           mean = Salix_pulchra_preds_df[,2], 
                                           upper = Salix_pulchra_preds_df[,3], year = seq(1:nyears))

# Salix arctica
salarc_m <- MCMCglmm(rw ~ I(Year-1973), random = ~ Year + IndivUN, data = dendro.salarc, 
                     family = "gaussian", pr=TRUE, nitt = 100000, burnin = 20000, prior = prior2) 

summary(salarc_m)
# plot(salarc_m$VCV)
# plot(salarc_m$Sol)
autocorr(salarc_m$VCV)

nyears <- length(unique(dendro.salarc$Year))
niter <- length(salarc_m$Sol[,"(Intercept)"])
Salix_arctica_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    Salix_arctica_preds[i,j] <- salarc_m$Sol[i,"(Intercept)"] + salarc_m$Sol[i,"I(Year - 1973)"]*j
  }
}

Salix_arctica_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  Salix_arctica_preds_df[i,] <- quantile(Salix_arctica_preds[,i], c(0.025, 0.5, 0.975))
}

Salix_arctica_preds_df <- cbind.data.frame(lower = Salix_arctica_preds_df[,1], 
                                           mean = Salix_arctica_preds_df[,2], 
                                           upper = Salix_arctica_preds_df[,3], 
                                           year = seq(1:nyears))

# Salix glauca
salgla_m <- MCMCglmm(rw ~ I(Year-1973), random = ~ Year + IndivUN, 
                     data = dendro.salgla, family="gaussian", pr=TRUE, 
                     nitt = 100000, burnin = 20000, prior = prior2) 

summary(salgla_m)
# plot(salgla_m$VCV)
# plot(salgla_m$Sol)
autocorr(salgla_m$VCV) 

nyears <- length(unique(dendro.salgla$Year))
niter <- length(salgla_m$Sol[,"(Intercept)"])
Salix_glauca_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    Salix_glauca_preds[i,j] <- salgla_m$Sol[i,"(Intercept)"] + salgla_m$Sol[i,"I(Year - 1973)"]*j
  }
}

Salix_glauca_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  Salix_glauca_preds_df[i,] <- quantile(Salix_glauca_preds[,i], c(0.025, 0.5, 0.975))
}

Salix_glauca_preds_df <- cbind.data.frame(lower = Salix_glauca_preds_df[,1], 
                                          mean = Salix_glauca_preds_df[,2], 
                                          upper = Salix_glauca_preds_df[,3], year = seq(1:nyears))

(radial.growth <- ggplot() +
    geom_point(data = dendro, aes(x = Year, y = rw, colour = factor(Sp), 
                                  fill = factor(Sp)), alpha = 0.8, size = 4) +
    scale_y_continuous(limits = c(0, 3)) +
    scale_color_manual(breaks=c('Salix arctica','Salix glauca', 'Salix pulchra','Salix richardsonii'),
                       labels=c("S. arctica","S. glauca", "S. pulchra","S. richardsonii"),
                       values=c("#66c2a5","#edb81c","#006b77","#7c337c")) +
    scale_fill_manual(values=c("#66c2a5","#edb81c","#006b77","#7c337c")) +
    geom_ribbon(data = Salix_pulchra_preds_df, aes(x = year + 1973, ymin = lower, ymax = upper),
                fill = "#008c5f", alpha = 0.2) +
    geom_ribbon(data = Salix_richardsonii_preds_df, aes(x = year + 1973, ymin = lower, ymax = upper),
                fill = "#CD96CD", alpha = 0.2) +
    geom_ribbon(data = Salix_arctica_preds_df, aes(x = year + 1973, ymin = lower, ymax = upper),
                fill = "#66c2a5", alpha = 0.2) +
    geom_ribbon(data = Salix_glauca_preds_df, aes(x = year + 1973, ymin = lower, ymax = upper),
                fill = "#EEDC82", alpha = 0.2) +
    geom_line(data = Salix_arctica_preds_df, aes(x = year + 1973, y = mean), colour = "#66c2a5") +
    geom_line(data = Salix_glauca_preds_df, aes(x = year + 1973, y = mean), colour = "#edb81c") +
    geom_line(data = Salix_pulchra_preds_df, aes(x = year + 1973, y = mean), colour = "#006b77") +
    geom_line(data = Salix_richardsonii_preds_df, aes (x = year + 1973, y = mean), colour = "#7c337c") +
    guides(fill = FALSE) +
    theme_QHI() +
    theme(legend.position = "top", 
          legend.text = element_text(size = 18, face="italic"),
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    labs(x = " ", y = "Radial growth (mm)\n", title = "(c) Radial growth"))

# Combined canopy height, salix height and radial growth figure
# Arranging in a panel and saving the file

Figure7 <- grid.arrange(grid.arrange(canopy.height, salix, ncol = 2, nrow = 1),
                        radial.growth, ncol = 1, nrow = 2)

ggsave("Qikiqtaruk_manuscript/figures/Figure7_height_dendro.pdf", 
       Figure7, width = 30, height = 30, units = "cm")

# Figure 8. Vegetation cover community composition changes ----

# Biomass (Vegetation cover index) model for HE
biomass_HE_m <- MCMCglmm(Biomass ~ I(YEAR-1998), random = ~ YEAR + PLOT, 
                         data = biomass_hits[biomass_hits$SUBSITE == "HE",], 
                         family="gaussian", pr=TRUE, nitt = 100000, burnin = 20000, prior=prior2)
summary(biomass_HE_m)
# plot(biomass_HE_m$VCV)
# plot(biomass_HE_m$Sol)
autocorr(biomass_HE_m$VCV)

# Calculate model predictions
nyears <- 19
niter <- length(biomass_HE_m$Sol[,"(Intercept)"])

biomass_HE_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    biomass_HE_preds[i,j] <- biomass_HE_m$Sol[i,"(Intercept)"] + biomass_HE_m$Sol[i,"I(YEAR - 1998)"]*j
  }
}

biomass_HE_preds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  biomass_HE_preds_df[i,] <- quantile(biomass_HE_preds[,i], c(0.025, 0.5, 0.975))
}

biomass_HE_preds_df <- cbind.data.frame(lower = biomass_HE_preds_df[,1], 
                                        mean = biomass_HE_preds_df[,2], upper = biomass_HE_preds_df[,3], year = seq(1:19))

# Biomass (Vegetation cover index) model for KO
biomass_KO_m <- MCMCglmm(Biomass ~ I(YEAR - 1998), random = ~ YEAR + PLOT, data = biomass_hits[biomass_hits$SUBSITE == "KO",], family = "gaussian", pr = TRUE, nitt = 100000, burnin = 20000, prior = prior2)
summary(biomass_KO_m)
# plot(biomass_KO_m$VCV)
# plot(biomass_KO_m$Sol)
autocorr(biomass_KO_m$VCV)

# Calculate model predictions
nyears <- 19
niter <- length(biomass_KO_m$Sol[,"(Intercept)"])

biomass_KO_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    biomass_KO_preds[i,j] <- biomass_KO_m$Sol[i,"(Intercept)"] + biomass_KO_m$Sol[i,"I(YEAR - 1998)"]*j
  }
}

biomass_KO_preds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  biomass_KO_preds_df[i,] <- quantile(biomass_KO_preds[,i], c(0.025, 0.5, 0.975))
}

biomass_KO_preds_df <- cbind.data.frame(lower = biomass_KO_preds_df[,1], 
                                        mean = biomass_KO_preds_df[,2], upper = biomass_KO_preds_df[,3], year = seq(1:19))

# Biomass (Vegetation cover index) graph
(veg.cover <- ggplot() +
    geom_point(data = biomass_hits, aes(x = YEAR, y = Biomass, colour = factor(SUBSITE)), alpha = 0.8, size = 4) +
    scale_color_manual(values = c("#ffa544", "#2b299b"), name = "", labels = c("Her.", "Kom.")) +
    scale_fill_manual(values = c("#ffa544","#2b299b")) +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2013, 2017)) +
    scale_y_continuous(breaks = c(0, 2.5, 5, 7.5, 10)) +
    geom_ribbon(data = biomass_HE_preds_df, aes(x = year + 1998, ymin = lower, ymax = upper), 
                fill = "#ffa544", alpha = 0.2) +
    geom_line(data = biomass_HE_preds_df, aes(x = year + 1998, y = mean), colour = "#ffa544") +
    geom_ribbon(data = biomass_KO_preds_df, aes(x = year + 1998, ymin = lower, ymax = upper), 
                fill = "#2b299b", alpha = 0.2) +
    geom_line(data = biomass_KO_preds_df, aes(x = year + 1998, y = mean), colour = "#2b299b") +
    theme_QHI() +
    theme(legend.position = c(0.1, 0.95), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    labs(x = "", y = "Vegetation cover index\n", title = "(a) Vegetation cover\n"))

# Bare ground HE model

# Creating a column for the trial hits and misses

bareground$failures <- 100 - bareground$Bareground
bareground_HE_m <- MCMCglmm(cbind(Bareground, failures) ~ I(YEAR - 1998), 
                            random = ~YEAR + PLOT, data = bareground[bareground$SUBSITE == "HE",], 
                            family = "multinomial2", pr = TRUE, nitt = 100000, burnin = 20000, prior = prior2)
summary(bareground_HE_m)
# plot(bareground_HE_m$VCV)
# plot(bareground_HE_m$Sol)
autocorr(bareground_HE_m$VCV)

# Calculate model predictions
nyears <- 19
niter <- length(bareground_HE_m$Sol[,"(Intercept)"])

bareground_HE_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    bareground_HE_preds[i,j] <- 100*plogis(bareground_HE_m$Sol[i,"(Intercept)"] + bareground_HE_m$Sol[i,"I(YEAR - 1998)"]*j)
  }
}

bareground_HE_preds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  bareground_HE_preds_df[i,] <- quantile(bareground_HE_preds[,i], c(0.025, 0.5, 0.975))
}

bareground_HE_preds_df <- cbind.data.frame(lower = bareground_HE_preds_df[,1], 
                                           mean = bareground_HE_preds_df[,2], upper = bareground_HE_preds_df[,3], year = seq(1:19))

# Bare ground KO model
bareground_KO_m <- MCMCglmm(cbind(Bareground, failures) ~ I(YEAR - 1998), 
                            random = ~YEAR + PLOT, data = bareground[bareground$SUBSITE == "KO",], 
                            family = "multinomial2", pr = TRUE, nitt = 100000, burnin = 20000, prior = prior2)
summary(bareground_KO_m)
# plot(bareground_KO_m$VCV)
# plot(bareground_KO_m$Sol)
autocorr(bareground_KO_m$VCV)

# Calculate model predictions
nyears <- 19
niter <- length(bareground_KO_m$Sol[,"(Intercept)"])

bareground_KO_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    bareground_KO_preds[i,j] <- 100*plogis(bareground_KO_m$Sol[i,"(Intercept)"] + bareground_KO_m$Sol[i,"I(YEAR - 1998)"]*j)
  }
}

bareground_KO_preds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  bareground_KO_preds_df[i,] <- quantile(bareground_KO_preds[,i], c(0.025, 0.5, 0.975))
}

bareground_KO_preds_df <- cbind.data.frame(lower = bareground_KO_preds_df[,1], 
                                           mean = bareground_KO_preds_df[,2], upper = bareground_KO_preds_df[,3], year = seq(1:19))

# Bare ground change plot
(bare.ground <- ggplot() +
    geom_point(data = bareground, aes(x = YEAR, y = Bareground, colour = factor(SUBSITE)), alpha = 0.8, size = 4) +
    scale_color_manual(values = c("#ffa544", "#2b299b"), name = "", labels = c("Her.", "Kom.")) +
    scale_fill_manual(values = c("#ffa544","#2b299b")) +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2013, 2017)) +
    geom_ribbon(data = bareground_HE_preds_df, aes(x = year + 1998, ymin = lower, ymax = upper), 
                fill = "#ffa544", alpha = 0.2) +
    geom_line(data = bareground_HE_preds_df, aes(x = year + 1998, y = mean), colour = "#ffa544") +
    geom_ribbon(data = bareground_KO_preds_df, aes(x = year + 1998, ymin = lower, ymax = upper), 
                fill = "#2b299b", alpha = 0.2) +
    geom_line(data = bareground_KO_preds_df, aes(x = year + 1998, y = mean), colour = "#2b299b") +
    theme_QHI() +
    theme(legend.position = c(0.9, 0.95), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    labs(x = "", y = "Bare ground (% cover per plot)\n", title = "(b) Bare ground\n"))

# Species richness model

#Site alpha diversity
alpha_site <- ddply(diversity,.(sub_name, year), summarise,
                    richness = length(unique(name)))
alpha_site$plot_unique <- paste(alpha_site$sub_name, alpha_site$plot, sep = "")

# Richness HE model
richness_HE_m <- MCMCglmm(richness ~ I(year - 1998), data = alpha_site[alpha_site$sub_name == "QHI:HE",], 
                          family = "gaussian", pr = TRUE, nitt = 100000, burnin = 20000)
summary(richness_HE_m)
# plot(richness_HE_m$VCV)
# plot(richness_HE_m$Sol)
autocorr(richness_HE_m$VCV)

# Calculate model predictions
nyears <- 19
niter <- length(richness_HE_m$Sol[,"(Intercept)"])

richness_HE_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    richness_HE_preds[i,j] <- richness_HE_m$Sol[i,"(Intercept)"] + richness_HE_m$Sol[i,"I(year - 1998)"]*j
  }
}

richness_HE_preds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  richness_HE_preds_df[i,] <- quantile(richness_HE_preds[,i], c(0.025, 0.5, 0.975))
}

richness_HE_preds_df <- cbind.data.frame(lower = richness_HE_preds_df[,1], 
                                         mean = richness_HE_preds_df[,2], upper = richness_HE_preds_df[,3], year = seq(1:19))

# Richness KO model
richness_KO_m <- MCMCglmm(richness ~ I(year - 1998), data = alpha_site[alpha_site$sub_name == "QHI:KO",], 
                          family = "gaussian", pr = TRUE, nitt = 100000, burnin = 20000)
summary(richness_KO_m)
# plot(richness_KO_m$VCV)
# plot(richness_KO_m$Sol)
autocorr(richness_KO_m$VCV)

# Calculate model predictions
nyears <- 19
niter <- length(richness_KO_m$Sol[,"(Intercept)"])

richness_KO_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    richness_KO_preds[i,j] <- richness_KO_m$Sol[i,"(Intercept)"] + richness_KO_m$Sol[i,"I(year - 1998)"]*j
  }
}

richness_KO_preds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  richness_KO_preds_df[i,] <- quantile(richness_KO_preds[,i], c(0.025, 0.5, 0.975))
}

richness_KO_preds_df <- cbind.data.frame(lower = richness_KO_preds_df[,1], 
                                         mean = richness_KO_preds_df[,2], upper = richness_KO_preds_df[,3], year = seq(1:19))

# Richness plot
(richness.plot <- ggplot() +
    geom_point(data = alpha_site, aes(x = year, y = richness, colour = factor(sub_name)), 
               alpha = 0.8, size = 4, position = position_jitter(height = 0.3, width = 0.3)) +
    scale_color_manual(values = c("#ffa544", "#2b299b"), name = "", labels = c("Her.", "Kom.")) +
    scale_fill_manual(values = c("#ffa544","#2b299b", labels = c("Her.", "Kom."))) +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2013, 2017)) +
    geom_ribbon(data = richness_HE_preds_df, aes(x = year + 1998, ymin = lower, ymax = upper), 
                fill = "#ffa544", alpha = 0.2) +
    geom_line(data = richness_HE_preds_df, aes(x = year + 1998, y = mean), colour = "#ffa544") +
    geom_ribbon(data = richness_KO_preds_df, aes(x = year + 1998, ymin = lower, ymax = upper), 
                fill = "#2b299b", alpha = 0.2) +
    geom_line(data = richness_KO_preds_df, aes(x = year + 1998, y = mean), colour = "#2b299b") +
    theme_QHI() +
    coord_cartesian(ylim = c(20, 45), xlim = c(1999, 2017)) +
    theme(legend.position = c(0.1, 0.95), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    labs(x = "", y = "Species richness\n", title = "(c) Species richness\n"))

# Evenness calculation
diversity2 <- ddply(diversity,.(sub_name, plot, year, name), summarise,
               plnp = mean(rel_cover/100 * (log(rel_cover/100))),
               p2 = mean((rel_cover/100)^2))
diversity2 <- subset(diversity2, !is.na(plnp))
Indices <- ddply(diversity2,.(sub_name, year, plot), summarise,
                 lnS = log(length(unique(name))),
                 Shannon = sum(plnp)*-1,
                 Simpson = 1-(sum(p2)),
                 Evenness = Shannon / lnS)

#Add back in unique plots
Indices$PLOT2 <- paste(Indices$sub_name, Indices$plot, sep="")

Indices_HE <- subset(Indices, sub_name=="QHI:HE")
Indices_KO <- subset(Indices, sub_name=="QHI:KO")

# Evenness HE model
# Creating a column for the trial hits and misses
Indices_HE$failures <- round(100 - 100*Indices_HE$Evenness)

evenness_HE_m <- MCMCglmm(cbind(round(100*Evenness), failures) ~ I(year - 1998), 
                          random = ~year + plot, data = Indices_HE, 
                          family = "multinomial2", nitt = 100000, burnin = 20000, 
                          pr = TRUE, prior = prior2)
summary(evenness_HE_m)
# plot(evenness_HE_m$VCV)
# plot(evenness_HE_m$Sol)
autocorr(evenness_HE_m$VCV)

# Calculate model predictions
nyears <- 19
niter <- length(evenness_HE_m$Sol[,"(Intercept)"])

evenness_HE_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    evenness_HE_preds[i,j] <- 100*plogis(evenness_HE_m$Sol[i,"(Intercept)"] + evenness_HE_m$Sol[i,"I(year - 1998)"]*j)/100
  }
}

evenness_HE_preds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  evenness_HE_preds_df[i,] <- quantile(evenness_HE_preds[,i], c(0.025, 0.5, 0.975))
}

evenness_HE_preds_df <- cbind.data.frame(lower = evenness_HE_preds_df[,1], 
                                         mean = evenness_HE_preds_df[,2], upper = evenness_HE_preds_df[,3], year = seq(1:19))

# Evenness KO model
# Creating a column for the trial hits and misses
Indices_KO$failures <- round(100 - 100*Indices_KO$Evenness)

evenness_KO_m <- MCMCglmm(cbind(round(100*Evenness), failures) ~ I(year - 1998), 
                          random = ~ year + plot, data = Indices_KO, 
                          family = "multinomial2", nitt = 100000, burnin = 20000, 
                          pr = TRUE, prior = prior2)
summary(evenness_KO_m)
# plot(evenness_KO_m$VCV)
# plot(evenness_KO_m$Sol)
autocorr(evenness_KO_m$VCV)

# Calculate model predictions
nyears <- 19
niter <- length(evenness_KO_m$Sol[,"(Intercept)"])

evenness_KO_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    evenness_KO_preds[i,j] <- 100*plogis(evenness_KO_m$Sol[i,"(Intercept)"] + evenness_KO_m$Sol[i,"I(year - 1998)"]*j)/100
  }
}

evenness_KO_preds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  evenness_KO_preds_df[i,] <- quantile(evenness_KO_preds[,i], c(0.025, 0.5, 0.975))
}

evenness_KO_preds_df <- cbind.data.frame(lower = evenness_KO_preds_df[,1], 
                                         mean = evenness_KO_preds_df[,2], upper = evenness_KO_preds_df[,3], year = seq(1:19))

# Evenness plot
(evenness.plot <- ggplot() +
    geom_point(data = Indices, aes(x = year, y = Evenness, colour = factor(sub_name)), alpha = 0.8, size = 4) +
    scale_color_manual(values = c("#ffa544", "#2b299b"), name = "", labels = c("Her.", "Kom.")) +
    scale_fill_manual(values = c("#ffa544","#2b299b")) +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2013, 2017)) +
    geom_ribbon(data = evenness_HE_preds_df, aes(x = year + 1998, ymin = lower, ymax = upper), 
                fill = "#ffa544", alpha = 0.2) +
    geom_line(data = evenness_HE_preds_df, aes(x = year + 1998, y = mean), colour = "#ffa544") +
    geom_ribbon(data = evenness_KO_preds_df, aes(x = year + 1998, ymin = lower, ymax = upper), 
                fill = "#2b299b", alpha = 0.2) +
    geom_line(data = evenness_KO_preds_df, aes(x = year + 1998, y = mean), colour = "#2b299b") +
    theme_QHI() +
    coord_cartesian(ylim = c(0.40, 1), xlim = c(1999, 2017)) +
    theme(legend.position = c(0.1, 0.95), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    labs(x = "", y = "Evenness\n", title = "(d) Evenness\n"))

# Herschel type model for E. vaginatum
HE_plots_ev <- subset(HEcoverHE, name == "Eriophorum vaginatum")

# Creating a column for the trial hits and misses
HE_plots_ev$failures <- round(100 - HE_plots_ev$cover)

HE_plot_evag <- MCMCglmm(cbind(round(cover), failures) ~ I(year - 1998), 
                         random = ~ year + plot, data = HE_plots_ev, 
                         family = "multinomial2", pr = TRUE, nitt = 100000, 
                         burnin = 20000, prior = prior2)
summary(HE_plot_evag)  
# plot(HE_plot_evag$VCV)  
# plot(HE_plot_evag$Sol)
autocorr(HE_plot_evag$VCV)

# Calculating model predictions
nyears <- 19
niter <- length(HE_plot_evag$Sol[,"(Intercept)"])

HE_plot_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    HE_plot_preds[i,j] <- 100*plogis(HE_plot_evag$Sol[i,"(Intercept)"] + HE_plot_evag$Sol[i,"I(year - 1998)"]*j)
  }
}

HE_plot_preds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  HE_plot_preds_df[i,] <- quantile(HE_plot_preds[,i], c(0.025, 0.5, 0.975))
}

HE_plot_preds_df <- cbind.data.frame(lower = HE_plot_preds_df[,1], 
                                     mean = HE_plot_preds_df[,2], upper = HE_plot_preds_df[,3], year = seq(1:19))

# Herschel type model for S. pulchra
HE_plots_sp <- subset(HEcoverHE, name == "Salix pulchra")

# Creating a column for the trial hits and misses
HE_plots_sp$failures <- round(100 - HE_plots_sp$cover)

HE_plot_salpul <- MCMCglmm(cbind(round(cover), failures) ~ I(year - 1998), 
                           random = ~ year + plot, data = HE_plots_sp, 
                           family = "multinomial2", pr = TRUE, nitt = 100000, 
                           burnin = 20000, prior = prior2)
summary(HE_plot_salpul)  
# plot(HE_plot_salpul$VCV)  
# plot(HE_plot_salpul$Sol)
autocorr(HE_plot_salpul$VCV)

# Calculating model predictions
nyears <- 19
niter <- length(HE_plot_salpul$Sol[,"(Intercept)"])

HE_plot_preds2 <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    HE_plot_preds2[i,j] <- 100*plogis(HE_plot_salpul$Sol[i,"(Intercept)"] + HE_plot_salpul$Sol[i,"I(year - 1998)"]*j)
  }
}

HE_plot_preds_df2 <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  HE_plot_preds_df2[i,] <- quantile(HE_plot_preds2[,i], c(0.025, 0.5, 0.975))
}

HE_plot_preds_df2 <- cbind.data.frame(lower = HE_plot_preds_df2[,1], 
                                      mean = HE_plot_preds_df2[,2], upper = HE_plot_preds_df2[,3], year = seq(1:19))

# Herschel type graph
(herschel <- ggplot() +
    geom_point(data = subset(HEcoverHE, name == "Eriophorum vaginatum" | name == "Salix pulchra"), 
               aes(x = year, y = cover, colour = factor(name)), alpha = 0.8, size = 4) +
    scale_color_manual(values = c("#fc8d62","#008c5f"), name = "", labels = c("Eriophorum vaginatum", "Salix pulchra")) +
    scale_fill_manual(values = c("#fc8d62","#008c5f"), name = "", labels = c("Eriophorum vaginatum", "Salix pulchra")) +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2013, 2017)) +
    geom_ribbon(data = HE_plot_preds_df, aes(x = year + 1998, ymin = lower, ymax = upper), fill = "#fc8d62", alpha = 0.2) +
    geom_line(data = HE_plot_preds_df, aes(x = year + 1998, y = mean), colour = "#fc8d62") +
    geom_ribbon(data = HE_plot_preds_df2, aes(x = year + 1998, ymin = lower, ymax = upper), fill = "#008c5f", alpha = 0.2) +
    geom_line(data = HE_plot_preds_df2, aes(x = year + 1998, y = mean), colour = "#008c5f") +
    theme_QHI() +
    theme(legend.position = c(0.315, 0.95), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    labs(x = "", y = "Species cover per plot\n", title = "(e) Herschel plots\n"))

# Komakuk type model

# Arctagrostis latifolia
KO_plots_al <- subset(HEcoverKO,name == "Arctagrostis latifolia")

# Creating a column for the trial hits and misses
KO_plots_al$failures <- round(100 - KO_plots_al$cover)

# Binomial model
KO_plot_arclat <- MCMCglmm(cbind(round(cover), failures) ~ I(year - 1998), 
                           random = ~ year + plot, data = KO_plots_al, 
                           family = "multinomial2", pr = TRUE, nitt = 100000, 
                           burnin = 20000, prior = prior2)
summary(KO_plot_arclat)
# plot(KO_plot_arclat$VCV)
# plot(KO_plot_arclat$Sol)
autocorr(KO_plot_arclat$VCV)

# Calculating model predictions - binomial
nyears <- 19
niter <- length(KO_plot_arclat$Sol[,"(Intercept)"])

KO_plot_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    KO_plot_preds[i,j] <- 100*plogis(KO_plot_arclat$Sol[i,"(Intercept)"] + KO_plot_arclat$Sol[i,"I(year - 1998)"]*j)
  }
}

KO_plot_preds_df <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  KO_plot_preds_df[i,] <- quantile(KO_plot_preds[,i], c(0.025, 0.5, 0.975))
}

KO_plot_preds_df <- cbind.data.frame(lower = KO_plot_preds_df[,1],
                                     mean = KO_plot_preds_df[,2], upper = KO_plot_preds_df[,3], year = seq(1:19))

# Komakuk type model for Alopecurus alpinus
KO_plots_aa <- subset(HEcoverKO, name == "Alopecurus alpinus")

# Creating a column for the trial hits and misses
KO_plots_aa$failures <- round(100 - KO_plots_aa$cover)

# Binomial model
KO_plot_aloalp <- MCMCglmm(cbind(round(cover), failures) ~ I(year - 1998), 
                           random = ~ year + plot, data = KO_plots_aa, 
                           family = "multinomial2", pr = TRUE, nitt = 100000, 
                           burnin = 20000, prior = prior2)
summary(KO_plot_aloalp)
# plot(KO_plot_aloalp$VCV)
# plot(KO_plot_aloalp$Sol)
autocorr(KO_plot_aloalp$VCV)

# Calculating model predictions - binomial
nyears <- 19
niter <- length(KO_plot_aloalp$Sol[,"(Intercept)"])

KO_plot_preds2 <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    KO_plot_preds2[i,j] <- 100*plogis(KO_plot_aloalp$Sol[i,"(Intercept)"] + KO_plot_aloalp$Sol[i,"I(year - 1998)"]*j)
  }
}

KO_plot_preds_df2 <- array(NA, dim = c(nyears,3))

for (i in 1:nyears){
  KO_plot_preds_df2[i,] <- quantile(KO_plot_preds2[,i], c(0.025, 0.5, 0.975))
}

KO_plot_preds_df2 <- cbind.data.frame(lower = KO_plot_preds_df2[,1],
                                      mean = KO_plot_preds_df2[,2], upper = KO_plot_preds_df2[,3], year = seq(1:19))

# Komakuk type graph
(komakuk <- ggplot() +
    geom_point(data = subset(HEcoverKO, name == "Arctagrostis latifolia" | name == "Alopecurus alpinus"), 
               aes(x = year, y = cover, colour = factor(name)), alpha = 0.8, size = 4) +
    scale_color_manual(values = c("#ffcd44","#1b74d3"), name = "", labels = c("Alopecurus alpinus", "Arctagrostis latifolia")) +
    scale_fill_manual(values = c("#ffcd44","#1b74d3"), name = "") +
    scale_x_continuous(breaks = c(1999, 2004, 2009, 2013, 2017)) +
    geom_ribbon(data = KO_plot_preds_df, aes(x = year + 1998, ymin = lower, ymax = upper), fill = "#1b74d3", alpha = 0.2) +
    geom_line(data = KO_plot_preds_df, aes(x = year + 1998, y = mean), colour = "#1b74d3") +
    geom_ribbon(data = KO_plot_preds_df2, aes(x = year + 1998, ymin = lower, ymax = upper), fill = "#ffcd44", alpha = 0.2) +
    geom_line(data = KO_plot_preds_df2, aes(x = year + 1998, y = mean), colour = "#ffcd44") +
    theme_QHI() +
    theme(legend.position = c(0.275, 0.95), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    labs(x = "", y = "Species cover per plot\n", title = "(f) Komakuk plots\n") +
    coord_cartesian(ylim = c(0, 100), xlim = c(1999, 2017)))

# Arranging in a panel and saving the file
Figure8 <- grid.arrange(veg.cover, bare.ground, richness.plot, evenness.plot, herschel, komakuk, ncol=2)

ggsave("Qikiqtaruk_manuscript/figures/Figure8_com_comp.pdf", Figure8,
       width = 30, height = 45, units = "cm")

# All species model Herschel Vegetation Type

# Creating a column for the trial hits and misses
HEcoverHE$failures <- round(100 - HEcoverHE$cover)

# Linear model
HE_plot_all_linear <- MCMCglmm(cover ~ I(year - 1998) * name - 1, 
                               random = ~ year + plot, data = HEcoverHE, 
                               family = "gaussian", pr = TRUE, nitt = 100000, 
                               burnin = 20000, prior = prior2)
summary(HE_plot_all_linear)

# Binomial model
HE_plot_all <- MCMCglmm(cbind(round(cover), failures) ~ I(year - 1998) * name, 
                        random = ~ year + plot, data = HEcoverHE, 
                        family = "multinomial2", pr = TRUE, nitt = 100000, 
                        burnin = 20000, prior = prior2)
summary(HE_plot_all)
# plot(HE_plot_all$VCV)
# plot(HE_plot_all$Sol)

# All species model Komakuk Vegetation Type

# Creating a column for the trial hits and misses
HEcoverKO$failures <- round(100 - HEcoverKO$cover)

# Linear model
KO_plot_all_linear <- MCMCglmm(cover ~ I(year - 1998) * name -1, random = ~ year + plot, data = HEcoverKO, 
                               family = "gaussian", pr = TRUE, nitt = 100000, burnin = 20000, prior = prior2)
summary(KO_plot_all_linear)

# Binomial model
KO_plot_all <- MCMCglmm(cbind(round(cover), failures) ~ I(year - 1998) * name, 
                        random = ~ year + plot, data = HEcoverKO, 
                        family = "multinomial2", pr = TRUE, nitt = 100000, 
                        burnin = 20000, prior = prior2)
summary(KO_plot_all)
# plot(HE_plot_all$VCV)
# plot(HE_plot_all$Sol)

# Figure 9. Species accumulation curves ----
# Species accumulation models

# Herschel
# The best convergence of the three! And autocorrelation is within the limits!
HE_SAC_m <- MCMCglmm(accumulated ~ log(distance.all), family = "gaussian", 
                      data = SAC.acc[SAC.acc$veg.type == "HE",], 
                     pr=TRUE, nitt = 100000, burnin = 20000)
summary(HE_SAC_m)
# plot(HE_SAC_m$VCV)
# plot(HE_SAC_m$Sol)
autocorr(HE_SAC_m$VCV)

# Calculating model predictions
ndist <- 90
niter <- length(HE_SAC_m$Sol[,"(Intercept)"])
HE_SAC_preds <- array(NA, dim = c(niter,ndist))

for (i in 1:niter){
  for (j in 1:ndist){
    HE_SAC_preds[i,j] <- HE_SAC_m$Sol[i,"(Intercept)"] + HE_SAC_m$Sol[i,"log(distance.all)"]*log(j)
  }
}

HE_SAC_preds_df <- array(NA, dim = c(ndist, 3))

for (i in 1:ndist){
  HE_SAC_preds_df[i,] <- quantile(HE_SAC_preds[,i], c(0.025, 0.5, 0.975))
}

HE_SAC_preds_df <- cbind.data.frame(lower = HE_SAC_preds_df[,1], 
                                       mean = HE_SAC_preds_df[,2], 
                                       upper = HE_SAC_preds_df[,3], 
                                    distance = seq(1:90))

# Adding on a predictin for 0.125 meters and 0 species at 0 meters
new_df <- array(NA, dim = c(2,4))
new_df[1,1] <- 0
new_df[1,2] <- 0
new_df[1,3] <- 0
new_df[1,4] <- 0
new_df[2,2] <- mean(HE_SAC_m$Sol[,"(Intercept)"] + HE_SAC_m$Sol[,"log(distance.all)"]*log(0.125))
new_df[2,4] <- 0.125
new_df[2,1] <- HPDinterval(HE_SAC_m$Sol[,"(Intercept)"] + HE_SAC_m$Sol[,"log(distance.all)"]*log(0.125))[,1]
new_df[2,3] <- HPDinterval(HE_SAC_m$Sol[,"(Intercept)"] + HE_SAC_m$Sol[,"log(distance.all)"]*log(0.125))[,2]
colnames(new_df) <- c("lower", "mean", "upper", "distance")
HE_SAC_preds_df <- rbind(HE_SAC_preds_df, new_df)

# Komakuk type
KO_SAC_m <- MCMCglmm(accumulated ~ log(distance.all), family = "gaussian", 
                     data = SAC.acc[SAC.acc$veg.type == "KO",], 
                     pr=TRUE, nitt = 100000, burnin = 20000)
summary(KO_SAC_m)
# plot(KO_SAC_m$VCV)
# plot(KO_SAC_m$Sol)
autocorr(KO_SAC_m$VCV)

# Calculating model predictions
ndist <- 105
niter <- length(KO_SAC_m$Sol[,"(Intercept)"])
KO_SAC_preds <- array(NA, dim = c(niter,ndist))

for (i in 1:niter){
  for (j in 1:ndist){
    KO_SAC_preds[i,j] <- KO_SAC_m$Sol[i,"(Intercept)"] + KO_SAC_m$Sol[i,"log(distance.all)"]*log(j)
  }
}

KO_SAC_preds_df <- array(NA, dim = c(ndist, 3))

for (i in 1:ndist){
  KO_SAC_preds_df[i,] <- quantile(KO_SAC_preds[,i], c(0.025, 0.5, 0.975))
}

KO_SAC_preds_df <- cbind.data.frame(lower = KO_SAC_preds_df[,1], 
                                    mean = KO_SAC_preds_df[,2], 
                                    upper = KO_SAC_preds_df[,3], 
                                    distance = seq(1:105))

# Adding on a predictin for 0.125 meters and 0 species at 0 meters
new_df <- array(NA, dim = c(2,4))
new_df[1,1] <- 0
new_df[1,2] <- 0
new_df[1,3] <- 0
new_df[1,4] <- 0
new_df[2,2] <- mean(KO_SAC_m$Sol[,"(Intercept)"] + KO_SAC_m$Sol[,"log(distance.all)"]*log(0.125))
new_df[2,4] <- 0.125
new_df[2,1] <- HPDinterval(KO_SAC_m$Sol[,"(Intercept)"] + KO_SAC_m$Sol[,"log(distance.all)"]*log(0.125))[,1]
new_df[2,3] <- HPDinterval(KO_SAC_m$Sol[,"(Intercept)"] + KO_SAC_m$Sol[,"log(distance.all)"]*log(0.125))[,2]
colnames(new_df) <- c("lower", "mean", "upper", "distance")
KO_SAC_preds_df <- rbind(KO_SAC_preds_df, new_df)

# Species pool figures with model predictions
(SAC.plot <- ggplot() +
    geom_point(data = SAC.acc, aes(x = distance.all, y = accumulated, 
                                   colour = factor(veg.type)), alpha = 0.8, size = 8) +
    geom_ribbon(data = HE_SAC_preds_df, aes(x = distance, ymin = lower, ymax = upper), 
               fill = "#ffa544", alpha = 0.2) +
    geom_line(data = HE_SAC_preds_df, aes(x = distance, y = mean), colour = "#ffa544") +
    geom_ribbon(data = KO_SAC_preds_df, aes(x = distance, ymin = lower, ymax = upper), 
              fill = "#2b299b", alpha = 0.2) +
    geom_line(data = KO_SAC_preds_df, aes(x = distance, y = mean), colour = "#2b299b") +
    scale_color_manual(values = c("#ffa544", "#2b299b"), name = "", labels = c("Her.", "Kom.")) +
    scale_fill_manual(values = c("#ffa544","#2b299b"), name = "", labels = c("Her.", "Kom.")) +
    theme_QHI() +
    coord_cartesian(ylim = c(0, 65), xlim = c(0, 110)) +
    scale_x_continuous(breaks = c(0, 30, 60, 90)) +
    scale_y_continuous(breaks = c(0, 20, 40, 60)) +
    theme(legend.position = c(0.1, 0.94), 
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5),
          axis.text = element_text(size = 30),
          axis.title = element_text(size = 36),
          legend.text = element_text(size = 26)) +
    labs(x = "Distance (m)", y = "Number of species\n") +
    guides(fill = FALSE))

# Saving the file
Figure9 <- grid.arrange(SAC.plot, ncol = 1)

ggsave("Qikiqtaruk_manuscript/figures/Figure9_SAC.pdf", 
       Figure9, width = 30, height = 30, units = "cm")

# Standardised effect sizes ----

all_effects <- data.frame(matrix("NA", nrow = 33, ncol = 9))
colnames(all_effects) <- c("Parameter", "Estimate", "SD", "Lower", "Upper", 
                           "Standardized", "Lower_std", "Upper_std", "Colour")
all_effects$Parameter <- c("", "Temperature change", "Temp. (Spring)", "Temp. (Summer)",
                        "Temp. (Autumn)", "Temp. (Winter)", "Frost frequency",
                        "Ice, snow and permafrost", "Sea ice", "Snow melt date",
                        "Soil temp.", "Active layer depth", "Phenology", 
                        "Leaf emergence (S. arctica)", "Flowering time (S. arctica)",
                        "Flowering time (D. integrifolia)", "Flowering time (E. vaginatum)",
                        "Leaf senescence (S. arctica)", "Growing season (S. arctica)",
                        "Growth", "Canopy height (Kom.)", "Canopy height (Her.)",
                        "Canopy height (S. pulchra)", "Shrub growthrings", "Cover change",
                        "Cover (Kom.)", "Cover (Her.)", "Bare ground (Kom.)", "Bare ground (Her.)",
                        "Species pool (Kom.)", "Species pool (Her.)", "Active layer depth (Kom.)", 
                        "Active layer depth (Her.)")

all_effects$Estimate <- as.numeric(all_effects$Estimate)
all_effects$SD <- as.numeric(all_effects$SD)
all_effects$Lower <- as.numeric(all_effects$Lower)
all_effects$Upper <- as.numeric(all_effects$Upper)
all_effects$Standardized <- as.numeric(all_effects$Standardized)
all_effects$Lower_std <- as.numeric(all_effects$Lower_std)
all_effects$Upper_std <- as.numeric(all_effects$Upper_std)
all_effects$Colour <- as.character(all_effects$Colour)

all_effects[1,2:8] <- 0
all_effects[1,9] <- as.character("#FFFFFF")
all_effects[2,2:8] <- 0
all_effects[2,9] <- as.character("#FFFFFF")

all_effects[3,2] <- mean(spring_m$Sol[,"I(season - 2000)"])
all_effects[3,3] <- sd(qhi.tmp.spring.2$meanT)
all_effects[3,4] <- HPDinterval(spring_m$Sol[,"I(season - 2000)"])[,1]
all_effects[3,5] <- HPDinterval(spring_m$Sol[,"I(season - 2000)"])[,2]
all_effects[3,6] <- mean(spring_m$Sol[,"I(season - 2000)"])/sd(qhi.tmp.spring.2$meanT)
all_effects[3,7] <- HPDinterval(spring_m$Sol[,"I(season - 2000)"])[,1]/sd(qhi.tmp.spring.2$meanT)
all_effects[3,8] <- HPDinterval(spring_m$Sol[,"I(season - 2000)"])[,2]/sd(qhi.tmp.spring.2$meanT)
all_effects[3,9] <- as.character("#b70000")

all_effects[4,2] <- mean(summer_m$Sol[,"I(season - 2000)"])
all_effects[4,3] <- sd(qhi.tmp.summer.2$meanT)
all_effects[4,4] <- HPDinterval(summer_m$Sol[,"I(season - 2000)"])[,1]
all_effects[4,5] <- HPDinterval(summer_m$Sol[,"I(season - 2000)"])[,2]
all_effects[4,6] <- mean(summer_m$Sol[,"I(season - 2000)"])/sd(qhi.tmp.summer.2$meanT)
all_effects[4,7] <- HPDinterval(summer_m$Sol[,"I(season - 2000)"])[,1]/sd(qhi.tmp.summer.2$meanT)
all_effects[4,8] <- HPDinterval(summer_m$Sol[,"I(season - 2000)"])[,2]/sd(qhi.tmp.summer.2$meanT)
all_effects[4,9] <- as.character("#b70000")

all_effects[5,2] <- mean(fall_m$Sol[,"I(season - 2000)"])
all_effects[5,3] <- sd(qhi.tmp.fall.2$meanT)
all_effects[5,4] <- HPDinterval(fall_m$Sol[,"I(season - 2000)"])[,1]
all_effects[5,5] <- HPDinterval(fall_m$Sol[,"I(season - 2000)"])[,2]
all_effects[5,6] <- mean(fall_m$Sol[,"I(season - 2000)"])/sd(qhi.tmp.fall.2$meanT)
all_effects[5,7] <- HPDinterval(fall_m$Sol[,"I(season - 2000)"])[,1]/sd(qhi.tmp.fall.2$meanT)
all_effects[5,8] <- HPDinterval(fall_m$Sol[,"I(season - 2000)"])[,2]/sd(qhi.tmp.fall.2$meanT)
all_effects[5,9] <- as.character("#b70000")

all_effects[6,2] <- mean(winter_m$Sol[,"I(season - 2000)"])
all_effects[6,3] <- sd(qhi.tmp.winter.2$meanT)
all_effects[6,4] <- HPDinterval(winter_m$Sol[,"I(season - 2000)"])[,1]
all_effects[6,5] <- HPDinterval(winter_m$Sol[,"I(season - 2000)"])[,2]
all_effects[6,6] <- mean(winter_m$Sol[,"I(season - 2000)"])/sd(qhi.tmp.winter.2$meanT)
all_effects[6,7] <- HPDinterval(winter_m$Sol[,"I(season - 2000)"])[,1]/sd(qhi.tmp.winter.2$meanT)
all_effects[6,8] <- HPDinterval(winter_m$Sol[,"I(season - 2000)"])[,2]/sd(qhi.tmp.winter.2$meanT)
all_effects[6, 9] <- as.character("#b70000")

all_effects[7,2] <- mean(frs_m$Sol[,"I(year - 1979)"])
all_effects[7,3] <- sd(qhi.frs.sum$sum.ice)
all_effects[7,4] <- HPDinterval(frs_m$Sol[,"I(year - 1979)"])[,1]
all_effects[7,5] <- HPDinterval(frs_m$Sol[,"I(year - 1979)"])[,2]
all_effects[7,6] <- mean(frs_m$Sol[,"I(year - 1979)"])/sd(qhi.frs.sum$sum.ice)
all_effects[7,7] <- HPDinterval(frs_m$Sol[,"I(year - 1979)"])[,1]/sd(qhi.frs.sum$sum.ice)
all_effects[7,8] <- HPDinterval(frs_m$Sol[,"I(year - 1979)"])[,2]/sd(qhi.frs.sum$sum.ice)
all_effects[7,9] <- as.character("#b70000")

all_effects[8,2:8] <- 0
all_effects[8,9] <- as.character("#FFFFFF")

all_effects[9,2] <- mean(ice_m$Sol[,"I(Year - 1979)"])
all_effects[9,3] <- sd(qhi.frs.sum$sum.ice)
all_effects[9,4] <- HPDinterval(ice_m$Sol[,"I(Year - 1979)"])[,1]
all_effects[9,5] <- HPDinterval(ice_m$Sol[,"I(Year - 1979)"])[,2]
all_effects[9,6] <- mean(ice_m$Sol[,"I(Year - 1979)"])/sd(qhi.frs.sum$sum.ice)
all_effects[9,7] <- HPDinterval(ice_m$Sol[,"I(Year - 1979)"])[,1]/sd(qhi.frs.sum$sum.ice)
all_effects[9,8] <- HPDinterval(ice_m$Sol[,"I(Year - 1979)"])[,2]/sd(qhi.frs.sum$sum.ice)
all_effects[9,9] <- as.character("#3cd0ea")

all_effects[10,2] <- mean(snow_m$Sol[,"I(Year - 2000)"])
all_effects[10,3] <- sd(snow$P1)
all_effects[10,4] <- HPDinterval(snow_m$Sol[,"I(Year - 2000)"])[,1]
all_effects[10,5] <- HPDinterval(snow_m$Sol[,"I(Year - 2000)"])[,2]
all_effects[10,6] <- mean(snow_m$Sol[,"I(Year - 2000)"])/sd(snow$P1)
all_effects[10,7] <- HPDinterval(snow_m$Sol[,"I(Year - 2000)"])[,1]/sd(snow$P1)
all_effects[10,8] <- HPDinterval(snow_m$Sol[,"I(Year - 2000)"])[,2]/sd(snow$P1)
all_effects[10,9] <- as.character("#3cd0ea")

all_effects[11,2] <- mean(soil15_15_m$Sol[,"squared"])
all_effects[11,3] <- sd(soil.temp.qhi15.15$minTemp)
all_effects[11,4] <- HPDinterval(soil15_15_m$Sol[,"squared"])[,1]
all_effects[11,5] <- HPDinterval(soil15_15_m$Sol[,"squared"])[,2]
all_effects[11,6] <- mean(soil15_15_m$Sol[,"squared"])/sd(soil.temp.qhi15.15$minTemp)
all_effects[11,7] <- HPDinterval(soil15_15_m$Sol[,"squared"])[,1]/sd(soil.temp.qhi15.15$minTemp)
all_effects[11,8] <- HPDinterval(soil15_15_m$Sol[,"squared"])[,2]/sd(soil.temp.qhi15.15$minTemp)
all_effects[11,9] <- as.character("#3cd0ea")

all_effects[12,2] <- mean(act_layer_m_year$Sol[,"I(year - 1985)"])
all_effects[12,3] <- sd(act.layer$value, na.rm = TRUE)
all_effects[12,4] <- HPDinterval(act_layer_m_year$Sol[,"I(year - 1985)"])[,1]
all_effects[12,5] <- HPDinterval(act_layer_m_year$Sol[,"I(year - 1985)"])[,2]
all_effects[12,6] <- mean(act_layer_m_year$Sol[,"I(year - 1985)"])/sd(act.layer$value, na.rm = TRUE)
all_effects[12,7] <- HPDinterval(act_layer_m_year$Sol[,"I(year - 1985)"])[,1]/sd(act.layer$value, na.rm = TRUE)
all_effects[12,8] <- HPDinterval(act_layer_m_year$Sol[,"I(year - 1985)"])[,2]/sd(act.layer$value, na.rm = TRUE)
all_effects[12,9] <- as.character("#3cd0ea")

all_effects[13,2:8] <- 0
all_effects[13,9] <- as.character("#FFFFFF")

# Effect sizes for leaf emergence S arctica
# Read in output saved in phenology section
pheno.es <- read.csv("Qikiqtaruk_manuscript/model_outputs/Phenology_JAGS_EffectSizes_all_models.csv", stringsAsFactors = F)
all_effects[14,2] <- pheno.es$mean[pheno.es$Model=="P2"]
all_effects[14,3] <- pheno.es$Data_SD[pheno.es$Model=="P2"]
all_effects[14,4] <- pheno.es$X2.5.[pheno.es$Model=="P2"]
all_effects[14,5] <- pheno.es$X97.5.[pheno.es$Model=="P2"]
all_effects[14,6] <- pheno.es$mean[pheno.es$Model=="P2"]/pheno.es$Data_SD[pheno.es$Model=="P2"]
all_effects[14,7] <- pheno.es$X2.5.[pheno.es$Model=="P2"]/pheno.es$Data_SD[pheno.es$Model=="P2"]
all_effects[14,8] <- pheno.es$X97.5.[pheno.es$Model=="P2"]/pheno.es$Data_SD[pheno.es$Model=="P2"]
all_effects[14,9] <- as.character("#66c2a5")

# Effect sizes for flowering time S arctica
all_effects[15,2] <- pheno.es$mean[pheno.es$Model=="P3"][3]
all_effects[15,3] <- pheno.es$Data_SD[pheno.es$Model=="P3"][3]
all_effects[15,4] <- pheno.es$X2.5.[pheno.es$Model=="P3"][3]
all_effects[15,5] <- pheno.es$X97.5.[pheno.es$Model=="P3"][3]
all_effects[15,6] <- pheno.es$mean[pheno.es$Model=="P3"][3]/pheno.es$Data_SD[pheno.es$Model=="P3"][3]
all_effects[15,7] <- pheno.es$X2.5.[pheno.es$Model=="P3"][3]/pheno.es$Data_SD[pheno.es$Model=="P3"][3]
all_effects[15,8] <- pheno.es$X97.5.[pheno.es$Model=="P3"][3]/pheno.es$Data_SD[pheno.es$Model=="P3"][3]
all_effects[15,9] <- as.character("#66c2a5")

# Effect sizes for flowering time D integrifolia
all_effects[16,2] <- pheno.es$mean[pheno.es$Model=="P3"][1]
all_effects[16,3] <- pheno.es$Data_SD[pheno.es$Model=="P3"][1]
all_effects[16,4] <- pheno.es$X2.5.[pheno.es$Model=="P3"][1]
all_effects[16,5] <- pheno.es$X97.5.[pheno.es$Model=="P3"][1]
all_effects[16,6] <- pheno.es$mean[pheno.es$Model=="P3"][1]/pheno.es$Data_SD[pheno.es$Model=="P3"][1]
all_effects[16,7] <- pheno.es$X2.5.[pheno.es$Model=="P3"][1]/pheno.es$Data_SD[pheno.es$Model=="P3"][1]
all_effects[16,8] <- pheno.es$X97.5.[pheno.es$Model=="P3"][1]/pheno.es$Data_SD[pheno.es$Model=="P3"][1]
all_effects[16,9] <- as.character("#6157d1")

# Effect sizes for flowering time E vaginatum
all_effects[17,2] <- pheno.es$mean[pheno.es$Model=="P3"][2]
all_effects[17,3] <- pheno.es$Data_SD[pheno.es$Model=="P3"][2]
all_effects[17,4] <- pheno.es$X2.5.[pheno.es$Model=="P3"][2]
all_effects[17,5] <- pheno.es$X97.5.[pheno.es$Model=="P3"][2]
all_effects[17,6] <- pheno.es$mean[pheno.es$Model=="P3"][2]/pheno.es$Data_SD[pheno.es$Model=="P3"][2]
all_effects[17,7] <- pheno.es$X2.5.[pheno.es$Model=="P3"][2]/pheno.es$Data_SD[pheno.es$Model=="P3"][2]
all_effects[17,8] <- pheno.es$X97.5.[pheno.es$Model=="P3"][2]/pheno.es$Data_SD[pheno.es$Model=="P3"][2]
all_effects[17,9] <- as.character("#fc8d62")

# Effect sizes for leaf senescence S arctica
all_effects[18,2] <- pheno.es$mean[pheno.es$Model=="P5"]
all_effects[18,3] <- pheno.es$Data_SD[pheno.es$Model=="P5"]
all_effects[18,4] <- pheno.es$X2.5.[pheno.es$Model=="P5"]
all_effects[18,5] <- pheno.es$X97.5.[pheno.es$Model=="P5"]
all_effects[18,6] <- pheno.es$mean[pheno.es$Model=="P5"]/pheno.es$Data_SD[pheno.es$Model=="P5"]
all_effects[18,7] <- pheno.es$X2.5.[pheno.es$Model=="P5"]/pheno.es$Data_SD[pheno.es$Model=="P5"]
all_effects[18,8] <- pheno.es$X97.5.[pheno.es$Model=="P5"]/pheno.es$Data_SD[pheno.es$Model=="P5"]
all_effects[18,9] <- as.character("#66c2a5")

# Effect sizes for growing season S arctica
all_effects[19,2] <- pheno.es$mean[pheno.es$Model=="P2P5"]
all_effects[19,3] <- pheno.es$Data_SD[pheno.es$Model=="P2P5"]
all_effects[19,4] <- pheno.es$X2.5.[pheno.es$Model=="P2P5"]
all_effects[19,5] <- pheno.es$X97.5.[pheno.es$Model=="P2P5"]
all_effects[19,6] <- pheno.es$mean[pheno.es$Model=="P2P5"]/pheno.es$Data_SD[pheno.es$Model=="P2P5"]
all_effects[19,7] <- pheno.es$X2.5.[pheno.es$Model=="P2P5"]/pheno.es$Data_SD[pheno.es$Model=="P2P5"]
all_effects[19,8] <- pheno.es$X97.5.[pheno.es$Model=="P2P5"]/pheno.es$Data_SD[pheno.es$Model=="P2P5"]
all_effects[19,9] <- as.character("#66c2a5")

all_effects[20,2:8] <- 0
all_effects[20,9] <- as.character("#FFFFFF")

all_effects[21,2] <- mean(KO_canopy_m$Sol[,"I(YEAR - 1998)"])
all_effects[21,3] <- sd(avg_heights[avg_heights$SUBSITE == "KO",]$Mean.C.H)
all_effects[21,4] <- HPDinterval(KO_canopy_m$Sol[,"I(YEAR - 1998)"])[,1]
all_effects[21,5] <- HPDinterval(KO_canopy_m$Sol[,"I(YEAR - 1998)"])[,2]
all_effects[21,6] <- mean(KO_canopy_m$Sol[,"I(YEAR - 1998)"])/sd(avg_heights[avg_heights$SUBSITE == "KO",]$Mean.C.H)
all_effects[21,7] <- HPDinterval(KO_canopy_m$Sol[,"I(YEAR - 1998)"])[,1]/sd(avg_heights[avg_heights$SUBSITE == "KO",]$Mean.C.H)
all_effects[21,8] <- HPDinterval(KO_canopy_m$Sol[,"I(YEAR - 1998)"])[,2]/sd(avg_heights[avg_heights$SUBSITE == "KO",]$Mean.C.H)
all_effects[21,9] <- as.character("#ffa544")

all_effects[22,2] <- mean(HE_canopy_m$Sol[,"I(YEAR - 1998)"])
all_effects[22,3] <- sd(avg_heights[avg_heights$SUBSITE == "HE",]$Mean.C.H)
all_effects[22,4] <- HPDinterval(HE_canopy_m$Sol[,"I(YEAR - 1998)"])[,1]
all_effects[22,5] <- HPDinterval(HE_canopy_m$Sol[,"I(YEAR - 1998)"])[,2]
all_effects[22,6] <- mean(HE_canopy_m$Sol[,"I(YEAR - 1998)"])/sd(avg_heights[avg_heights$SUBSITE == "HE",]$Mean.C.H)
all_effects[22,7] <- HPDinterval(HE_canopy_m$Sol[,"I(YEAR - 1998)"])[,1]/sd(avg_heights[avg_heights$SUBSITE == "HE",]$Mean.C.H)
all_effects[22,8] <- HPDinterval(HE_canopy_m$Sol[,"I(YEAR - 1998)"])[,2]/sd(avg_heights[avg_heights$SUBSITE == "HE",]$Mean.C.H)
all_effects[22,9] <- as.character("#2b299b")

all_effects[23,2] <- mean(Salix_canopy_m$Sol[,"I(YEAR - 1998)"])
all_effects[23,3] <- sd(salpuls$Height..cm.)
all_effects[23,4] <- HPDinterval(Salix_canopy_m$Sol[,"I(YEAR - 1998)"])[,1]
all_effects[23,5] <- HPDinterval(Salix_canopy_m$Sol[,"I(YEAR - 1998)"])[,2]
all_effects[23,6] <- mean(Salix_canopy_m$Sol[,"I(YEAR - 1998)"])/sd(salpuls$Height..cm.)
all_effects[23,7] <- HPDinterval(Salix_canopy_m$Sol[,"I(YEAR - 1998)"])[,1]/sd(salpuls$Height..cm.)
all_effects[23,8] <- HPDinterval(Salix_canopy_m$Sol[,"I(YEAR - 1998)"])[,2]/sd(salpuls$Height..cm.)
all_effects[23,9] <- as.character("#008c5f")

# Salix pulchra model
all_effects[24,2] <- mean(salpul_m$Sol[,"I(Year - 1973)"])
all_effects[24,3] <- sd(dendro.salpul$rw)
all_effects[24,4] <- HPDinterval(salpul_m$Sol[,"I(Year - 1973)"])[,1]
all_effects[24,5] <- HPDinterval(salpul_m$Sol[,"I(Year - 1973)"])[,2]
all_effects[24,6] <- mean(salpul_m$Sol[,"I(Year - 1973)"])/sd(dendro.salpul$rw)
all_effects[24,7] <- HPDinterval(salpul_m$Sol[,"I(Year - 1973)"])[,1]/sd(dendro.salpul$rw)
all_effects[24,8] <- HPDinterval(salpul_m$Sol[,"I(Year - 1973)"])[,2]/sd(dendro.salpul$rw)
all_effects[24,9] <- as.character("#008c5f")

all_effects[25,2:8] <- 0
all_effects[25,9] <- as.character("#FFFFFF")

all_effects[26,2] <- mean(biomass_KO_m$Sol[,"I(YEAR - 1998)"])
all_effects[26,3] <- sd(biomass_hits[biomass_hits$SUBSITE == "KO",]$Biomass)
all_effects[26,4] <- HPDinterval(biomass_KO_m$Sol[,"I(YEAR - 1998)"])[,1]
all_effects[26,5] <- HPDinterval(biomass_KO_m$Sol[,"I(YEAR - 1998)"])[,2]
all_effects[26,6] <- mean(biomass_KO_m$Sol[,"I(YEAR - 1998)"])/sd(biomass_hits[biomass_hits$SUBSITE == "KO",]$Biomass)
all_effects[26,7] <- HPDinterval(biomass_KO_m$Sol[,"I(YEAR - 1998)"])[,1]/sd(biomass_hits[biomass_hits$SUBSITE == "KO",]$Biomass)
all_effects[26,8] <- HPDinterval(biomass_KO_m$Sol[,"I(YEAR - 1998)"])[,2]/sd(biomass_hits[biomass_hits$SUBSITE == "KO",]$Biomass)
all_effects[26,9] <- as.character("#2b299b")

all_effects[27,2] <- mean(biomass_HE_m$Sol[,"I(YEAR - 1998)"])
all_effects[27,3] <- sd(biomass_hits[biomass_hits$SUBSITE == "HE",]$Biomass)
all_effects[27,4] <- HPDinterval(biomass_HE_m$Sol[,"I(YEAR - 1998)"])[,1]
all_effects[27,5] <- HPDinterval(biomass_HE_m$Sol[,"I(YEAR - 1998)"])[,2]
all_effects[27,6] <- mean(biomass_HE_m$Sol[,"I(YEAR - 1998)"])/sd(biomass_hits[biomass_hits$SUBSITE == "HE",]$Biomass)
all_effects[27,7] <- HPDinterval(biomass_HE_m$Sol[,"I(YEAR - 1998)"])[,1]/sd(biomass_hits[biomass_hits$SUBSITE == "HE",]$Biomass)
all_effects[27,8] <- HPDinterval(biomass_HE_m$Sol[,"I(YEAR - 1998)"])[,2]/sd(biomass_hits[biomass_hits$SUBSITE == "HE",]$Biomass)
all_effects[27,9] <- as.character("#ffa544")

all_effects[28,2] <- mean(bareground_KO_m$Sol[,"I(YEAR - 1998)"])
all_effects[28,3] <- sd(bareground[bareground$SUBSITE == "KO",]$Bareground)
all_effects[28,4] <- HPDinterval(bareground_KO_m$Sol[,"I(YEAR - 1998)"])[,1]
all_effects[28,5] <- HPDinterval(bareground_KO_m$Sol[,"I(YEAR - 1998)"])[,2]
all_effects[28,6] <- mean(bareground_KO_m$Sol[,"I(YEAR - 1998)"])/sd(bareground[bareground$SUBSITE == "KO",]$Bareground)
all_effects[28,7] <- HPDinterval(bareground_KO_m$Sol[,"I(YEAR - 1998)"])[,1]/sd(bareground[bareground$SUBSITE == "KO",]$Bareground)
all_effects[28,8] <- HPDinterval(bareground_KO_m$Sol[,"I(YEAR - 1998)"])[,2]/sd(bareground[bareground$SUBSITE == "KO",]$Bareground)
all_effects[28,9] <- as.character("#2b299b")

all_effects[29,2] <- mean(bareground_HE_m$Sol[,"I(YEAR - 1998)"])
all_effects[29,3] <- sd(bareground[bareground$SUBSITE == "HE",]$Bareground)
all_effects[29,4] <- HPDinterval(bareground_HE_m$Sol[,"I(YEAR - 1998)"])[,1]
all_effects[29,5] <- HPDinterval(bareground_HE_m$Sol[,"I(YEAR - 1998)"])[,2]
all_effects[29,6] <- mean(bareground_HE_m$Sol[,"I(YEAR - 1998)"])/sd(bareground[bareground$SUBSITE == "KO",]$Bareground)
all_effects[29,7] <- HPDinterval(bareground_HE_m$Sol[,"I(YEAR - 1998)"])[,1]/sd(bareground[bareground$SUBSITE == "KO",]$Bareground)
all_effects[29,8] <- HPDinterval(bareground_HE_m$Sol[,"I(YEAR - 1998)"])[,2]/sd(bareground[bareground$SUBSITE == "KO",]$Bareground)
all_effects[29,9] <- as.character("#ffa544")

all_effects[30,2] <- mean(KO_SAC_m$Sol[,"log(distance.all)"])
all_effects[30,3] <- sd(SAC.acc[SAC.acc$veg.type == "KO",]$accumulated)
all_effects[30,4] <- HPDinterval(KO_SAC_m$Sol[,"log(distance.all)"])[,1]
all_effects[30,5] <- HPDinterval(KO_SAC_m$Sol[,"log(distance.all)"])[,2]
all_effects[30,6] <- mean(KO_SAC_m$Sol[,"log(distance.all)"])/sd(SAC.acc[SAC.acc$veg.type == "KO",]$accumulated)
all_effects[30,7] <- HPDinterval(KO_SAC_m$Sol[,"log(distance.all)"])[,1]/sd(SAC.acc[SAC.acc$veg.type == "KO",]$accumulated)
all_effects[30,8] <- HPDinterval(KO_SAC_m$Sol[,"log(distance.all)"])[,2]/sd(SAC.acc[SAC.acc$veg.type == "KO",]$accumulated)
all_effects[30,9] <- as.character("#ffa544")

all_effects[31,2] <- mean(HE_SAC_m$Sol[,"log(distance.all)"])
all_effects[31,3] <- sd(SAC.acc[SAC.acc$veg.type == "HE",]$accumulated)
all_effects[31,4] <- HPDinterval(HE_SAC_m$Sol[,"log(distance.all)"])[,1]
all_effects[31,5] <- HPDinterval(HE_SAC_m$Sol[,"log(distance.all)"])[,2]
all_effects[31,6] <- mean(HE_SAC_m$Sol[,"log(distance.all)"])/sd(SAC.acc[SAC.acc$veg.type == "HE",]$accumulated)
all_effects[31,7] <- HPDinterval(HE_SAC_m$Sol[,"log(distance.all)"])[,1]/sd(SAC.acc[SAC.acc$veg.type == "HE",]$accumulated)
all_effects[31,8] <- HPDinterval(HE_SAC_m$Sol[,"log(distance.all)"])[,2]/sd(SAC.acc[SAC.acc$veg.type == "HE",]$accumulated)
all_effects[31,9] <- as.character("#2b299b")

all_effects[32,2] <- mean(act_layer_m_HE$Sol[,"I(day - 174)"])
all_effects[32,3] <- sd(act.layer.2017.HE$mean.depth)
all_effects[32,4] <- HPDinterval(act_layer_m_HE$Sol[,"I(day - 174)"])[,1]
all_effects[32,5] <- HPDinterval(act_layer_m_HE$Sol[,"I(day - 174)"])[,2]
all_effects[32,6] <- mean(act_layer_m_HE$Sol[,"I(day - 174)"])/sd(act.layer.2017.HE$mean.depth)
all_effects[32,7] <- HPDinterval(act_layer_m_HE$Sol[,"I(day - 174)"])[,1]/sd(act.layer.2017.HE$mean.depth)
all_effects[32,8] <- HPDinterval(act_layer_m_HE$Sol[,"I(day - 174)"])[,2]/sd(act.layer.2017.HE$mean.depth)
all_effects[32,9] <- as.character("#3cd0ea")

all_effects[33,2] <- mean(act_layer_m_KO$Sol[,"I(day - 174)"])
all_effects[33,3] <- sd(act.layer.2017.KO$mean.depth, na.rm = TRUE)
all_effects[33,4] <- HPDinterval(act_layer_m_KO$Sol[,"I(day - 174)"])[,1]
all_effects[33,5] <- HPDinterval(act_layer_m_KO$Sol[,"I(day - 174)"])[,2]
all_effects[33,6] <- mean(act_layer_m_KO$Sol[,"I(day - 174)"])/sd(act.layer.2017.KO$mean.depth, na.rm = TRUE)
all_effects[33,7] <- HPDinterval(act_layer_m_KO$Sol[,"I(day - 174)"])[,1]/sd(act.layer.2017.KO$mean.depth, na.rm = TRUE)
all_effects[33,8] <- HPDinterval(act_layer_m_KO$Sol[,"I(day - 174)"])[,2]/sd(act.layer.2017.KO$mean.depth, na.rm = TRUE)
all_effects[33,9] <- as.character("#3cd0ea")

write.csv(all_effects, file = "Qikiqtaruk_manuscript/model_outputs/All_standardised_effects_2017.csv")

# All MCMC model outputs for effect size figure and SI table ----

# Creating a summary table of model outputs for models with no random effects
dataList <- list(spring_m, summer_m, fall_m, winter_m, frs_m, snow_m, ice_m, 
                 soil15_15_m, act_layer_m_year, richness_HE_m, richness_KO_m)

# Create a list of input model names
dataListNames <- list("Spring", "Summer", "Fall", "Winter", "Frost frequency", 
                      "Snow melt days", "Sea ice concentration", "Minimum soil temperature (15 m)",
                      "Active layer depth", "Species richness HE", "Species richness KO")

# Get the clean.MCMC outputs and add modelName columns to each element for ID purposes
readyList <- mapply(cbind, lapply(dataList, clean.MCMC.2), "modelName" = dataListNames, SIMPLIFY = F)

# Turn the list of data.frames into one big data.frame
mcmc.outputs.simple <- as.data.frame(do.call(rbind, readyList))

# Creating a summary table of model outputs for models with random effects
dataList <- list(HE_canopy_m, KO_canopy_m, Salix_canopy_m, salric_m, salpul_m,
                 salarc_m, salgla_m, biomass_HE_m, biomass_KO_m,
                 bareground_HE_m, bareground_KO_m, evenness_HE_m, 
                 evenness_KO_m, HE_plot_evag, HE_plot_salpul,
                 KO_plot_arclat, KO_plot_aloalp)

# Create a list of input model names
dataListNames <- list("Herschel canopy height", "Komakuk canopy height", 
                      "Salix canopy height", "Salric radial growth", 
                      "Salpul radial growth", "Salarc radial growth", 
                      "Salgla radial growth", "Herschel biomass(veg cover index)", 
                      "Komakuk biomass(veg cover index)", "Herschel bareground", 
                      "Komakuk bareground", "Herschel evenness", "Komakuk evenness", 
                      "Herschel plots evag", "Herschel plots salpul", 
                      "Komakuk plots arclat", "Komakuk plots aloalp")

# Get the clean.MCMC outputs and add modelName columns to each element for ID purposes
readyList <- mapply(cbind, lapply(dataList, clean.MCMC), "modelName" = dataListNames, SIMPLIFY = F)

# Turn the list of data.frames into one big data.frame
mcmc.outputs.random <- as.data.frame(do.call(rbind, readyList))

# Combine outputs from models with and without random effects
mcmc.outputs <- rbind(mcmc.outputs.random, mcmc.outputs.simple)

# Write csv
write.csv(mcmc.outputs, file = "Qikiqtaruk_manuscript/model_outputs/All_mcmc_outputs_2017.csv")

# Generate html table
stargazer(mcmc.outputs, type = "html", summary = FALSE, digits = 2)

# Copy the html code, click on File/New/R HTML file
# Delete everything in the newly generated file except the <html> and </html> tags
# Click Knit, save the html file
# The file can now be opened with Word

# Tables for all KO and HE species
HE_plot_all_linear_outputs <- clean.MCMC(HE_plot_all_linear)
KO_plot_all_linear_outputs <- clean.MCMC(KO_plot_all_linear)

# Write csv
write.csv(HE_plot_all_linear_outputs, file = "Qikiqtaruk_manuscript/model_outputs/HE_plot_all_linear_outputs.csv")
write.csv(KO_plot_all_linear_outputs, file = "Qikiqtaruk_manuscript/model_outputs/KO_plot_all_linear_outputs.csv")

# Herschel species table
stargazer(HE_plot_all_linear_outputs, type = "html", summary = FALSE, digits = 2)

# Komakuk species table
stargazer(KO_plot_all_linear_outputs, type = "html", summary = FALSE, digits = 2)

# Figure 10. Effect Sizes ----

effects <- all_effects
effects$Parameter <- factor(effects$Parameter, levels=unique(effects$Parameter))
effects$Colour <- as.character(effects$Colour)
effects <- effects %>% filter(Parameter != "Species pool (Her.)", Parameter != "Species pool (Kom.)", Parameter != "Active layer depth (Kom.)", Parameter != "Active layer depth (Her.)")

Figure10 <- grid.arrange(
  ggplot(effects) +
    geom_hline(yintercept = 0) +
    geom_bar(aes(x = Parameter, y = Standardized, fill = Parameter),
             colour = "#000000", stat = "identity", position = "dodge") +
    geom_errorbar(aes(x = Parameter, ymin = Lower_std, ymax = Upper_std), width = 0.4) +
    theme_bw() +
    scale_fill_manual(values = effects$Colour, name = "") +
    scale_x_discrete(labels = effects$Parameter) +
    ylab("Standardized Effect Size") + xlab("") +
    theme(axis.line = element_line(colour = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(), 
          legend.position = "none", 
          axis.text.x = element_text(size = 20, angle = 45,
                                     vjust = 1, hjust = 1,
                                     colour = "black"),
          axis.text.y = element_text(hjust = 1, size = 20,
                                     colour = "black"),
          axis.title.x = element_text(size = 20, 
                                      margin = margin(20, 0, 0, 0),
                                      colour = "black"),
          axis.title.y = element_text(angle = 90, size = 24, face = "plain",
                                    margin = margin(0, 20, 0, 0), 
                                    colour = "black"),
          axis.ticks = element_blank(), 
          panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank(), 
          panel.grid.minor.y = element_blank(), 
          panel.grid.major.y = element_blank(), 
          legend.key = element_rect(colour = "white"), 
          legend.text = element_text(size = 18))
)

# Saving the file
ggsave("Qikiqtaruk_manuscript/figures/Figure10_effect_sizes.pdf", 
       Figure10, width = 45, height = 30, units = "cm")

# Models and figures not used in manuscript ----

# Seasonal temp plot all in one
(all_seasons <- ggplot() +
   geom_point(data = qhi.tmp.all[qhi.tmp.all$season>cutoff,], aes(x = season, y = meanT,
                                                                  colour = factor(Time)), 
              alpha = 0.8, size = 4) +
   scale_colour_manual(values = c("#EE7621", "#CD2990", "#458B00", "#00B2EE")) +
   scale_fill_manual(values = c("#EE7621", "#CD2990", "#458B00", "#00B2EE")) +
   geom_ribbon(data = spring_preds_df, aes(x = year + 2000, ymin = lower, ymax = upper), 
               fill = "#CD2990", alpha = 0.2) +
   geom_line(data = spring_preds_df, aes(x = year + 2000, y = mean), colour = "#CD2990") +
   geom_ribbon(data = summer_preds_df, aes(x = year + 2000, ymin = lower, ymax = upper), 
               fill = "#458B00", alpha = 0.2) +
   geom_line(data = summer_preds_df, aes(x = year + 2000, y = mean), colour = "#458B00") +
   geom_ribbon(data = fall_preds_df, aes(x = year + 2000, ymin = lower, ymax = upper), 
               fill = "#EE7621", alpha = 0.2) +
   geom_line(data = fall_preds_df, aes(x = year + 2000, y = mean), colour = "#EE7621") +
   geom_ribbon(data = winter_preds_df, aes(x = year + 2000, ymin = lower, ymax = upper), 
               fill = "#00B2EE", alpha = 0.2) +
   geom_line(data = winter_preds_df, aes(x = year + 2000, y = mean), colour = "#00B2EE") +
   theme_QHI() +
   theme(legend.position = "top") +
   theme(axis.line.x = element_line(color="black", size = 0.5),
         axis.line.y = element_line(color="black", size = 0.5)) +
   labs(x=" ", y = "Temperature (°C)\n", title = "(a) Mean temperature"))

# Herbivory

# Load data
QikHerbDF <- read.csv("Qikiqtaruk_manuscript/data/qhi_herbivores.csv")

# Filtering original data frame per species
caribou <- filter(QikHerbDF, Species == "Caribou")
muskox <- filter(QikHerbDF, Species == "Muskox")

# Caribou group size model
caribou_size_m <- MCMCglmm(Size ~ I(Year - 1987), data = caribou, 
                           family = "gaussian", pr=TRUE, nitt = 100000, 
                           burnin = 20000)
summary(caribou_size_m)
# plot(caribou_size_m$VCV)
# plot(caribou_size_m$Sol)
autocorr(caribou_size_m$VCV)

# Calculating model predictions
mean1 <- mean(caribou_size_m$Sol[,"(Intercept)"] + caribou_size_m$Sol[,"I(Year - 1987)"])
mean2 <- mean(caribou_size_m$Sol[,"(Intercept)"] + (23*caribou_size_m$Sol[,"I(Year - 1987)"]))
lower1 <- HPDinterval(caribou_size_m$Sol[,"(Intercept)"] + caribou_size_m$Sol[,"I(Year - 1987)"])[1]
upper1 <- HPDinterval(caribou_size_m$Sol[,"(Intercept)"] + caribou_size_m$Sol[,"I(Year - 1987)"])[2]
lower2 <- HPDinterval(caribou_size_m$Sol[,"(Intercept)"] + (23*caribou_size_m$Sol[,"I(Year - 1987)"]))[1] 
upper2 <- HPDinterval(caribou_size_m$Sol[,"(Intercept)"] + (23*caribou_size_m$Sol[,"I(Year - 1987)"]))[2] 

# Calculating model predictions
nyears <- 23
niter <- length(caribou_size_m$Sol[,"(Intercept)"])

caribou_size_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    caribou_size_preds[i,j] <- caribou_size_m$Sol[i,"(Intercept)"] + caribou_size_m$Sol[i,"I(Year - 1987)"]*j
  }
}

caribou_size_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  caribou_size_preds_df[i,] <- quantile(caribou_size_preds[,i], c(0.025, 0.5, 0.975))
}

caribou_size_preds_df <- cbind.data.frame(lower = caribou_size_preds_df[,1], 
                                          mean = caribou_size_preds_df[,2], upper = caribou_size_preds_df[,3], year = seq(1:23))

# Muskox group size model
muskox_size_m <- MCMCglmm(Size ~ I(Year - 1987), data = muskox, family = "gaussian", pr=TRUE, nitt = 100000, burnin = 20000)
summary(muskox_size_m)
# plot(muskox_size_m$VCV)
# plot(muskox_size_m$Sol)
autocorr(muskox_size_m$VCV)

# Calculating model predictions
mean1 <- mean(muskox_size_m$Sol[,"(Intercept)"] + muskox_size_m$Sol[,"I(Year - 1987)"])
mean2 <- mean(muskox_size_m$Sol[,"(Intercept)"] + (23*muskox_size_m$Sol[,"I(Year - 1987)"]))
lower1 <- HPDinterval(muskox_size_m$Sol[,"(Intercept)"] + muskox_size_m$Sol[,"I(Year - 1987)"])[1]
upper1 <- HPDinterval(muskox_size_m$Sol[,"(Intercept)"] + muskox_size_m$Sol[,"I(Year - 1987)"])[2]
lower2 <- HPDinterval(muskox_size_m$Sol[,"(Intercept)"] + (23*muskox_size_m$Sol[,"I(Year - 1987)"]))[1] 
upper2 <- HPDinterval(muskox_size_m$Sol[,"(Intercept)"] + (23*muskox_size_m$Sol[,"I(Year - 1987)"]))[2]

# Calculating model predictions
nyears <- 23
niter <- length(muskox_size_m$Sol[,"(Intercept)"])

muskox_size_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    muskox_size_preds[i,j] <- muskox_size_m$Sol[i,"(Intercept)"] + muskox_size_m$Sol[i,"I(Year - 1987)"]*j
  }
}

muskox_size_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  muskox_size_preds_df[i,] <- quantile(muskox_size_preds[,i], c(0.025, 0.5, 0.975))
}

muskox_size_preds_df <- cbind.data.frame(lower = muskox_size_preds_df[,1], 
                                         mean = muskox_size_preds_df[,2], 
                                         upper = muskox_size_preds_df[,3], year = seq(1:23))

# Max size of group spotted, caribou/muskox, 1988-2010
(group_size <- ggplot() + 
    geom_point(data = QikHerbDF, aes(Year, Size, colour = Species), alpha = 0.8, size = 4) +
    scale_color_manual(values=c("#ffcd44", "#1b74d3")) +
    scale_fill_manual(values=c("#ffcd44", "#1b74d3"), labels=c("Caribou", "Muskox")) +
    geom_ribbon(data = caribou_size_preds_df, aes(x = year + 1987, ymin = lower, ymax = upper),
                fill = "#ffcd44", alpha = 0.2) +
    geom_line(data = caribou_size_preds_df, aes(x = year + 1987, y = mean), colour = "#ffcd44") +
    geom_ribbon(data = muskox_size_preds_df, aes(x = year + 1987, ymin = lower, ymax = upper), 
                fill = "#1b74d3", alpha = 0.2) +
    geom_line(data = muskox_size_preds_df, aes(x = year + 1987, y = mean), colour = "#1b74d3") +
    labs(x=" ", y="Largest group size observed\n", title="(c) Herbivore group size") +
    theme_QHI() +
    theme(legend.position = c(0.15, 0.9), 
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.line.y = element_line(color = "black", size = 0.5)) +
    coord_cartesian(ylim = c(0, 80), xlim = c(1988, 2010)))

# Caribou number of groups model
caribou_number_m <- MCMCglmm(Groups ~ I(Year - 1987), data = caribou, 
                             family = "gaussian", pr=TRUE, nitt = 100000, 
                             burnin = 20000)
summary(caribou_number_m)
# plot(caribou_number_m$VCV)
# plot(caribou_number_m$Sol)
autocorr(caribou_number_m$VCV)

# Calculating model predictions
mean1 <- mean(caribou_number_m$Sol[,"(Intercept)"] + caribou_number_m$Sol[,"I(Year - 1987)"])
mean2 <- mean(caribou_number_m$Sol[,"(Intercept)"] + (23*caribou_number_m$Sol[,"I(Year - 1987)"]))
lower1 <- HPDinterval(caribou_number_m$Sol[,"(Intercept)"] + caribou_number_m$Sol[,"I(Year - 1987)"])[1]
upper1 <- HPDinterval(caribou_number_m$Sol[,"(Intercept)"] + caribou_number_m$Sol[,"I(Year - 1987)"])[2]
lower2 <- HPDinterval(caribou_number_m$Sol[,"(Intercept)"] + (23*caribou_number_m$Sol[,"I(Year - 1987)"]))[1] 
upper2 <- HPDinterval(caribou_number_m$Sol[,"(Intercept)"] + (23*caribou_number_m$Sol[,"I(Year - 1987)"]))[2] 

# Calculating model predictions
nyears <- 23
niter <- length(caribou_number_m$Sol[,"(Intercept)"])

caribou_number_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    caribou_number_preds[i,j] <- caribou_number_m$Sol[i,"(Intercept)"] + caribou_number_m$Sol[i,"I(Year - 1987)"]*j
  }
}

caribou_number_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  caribou_number_preds_df[i,] <- quantile(caribou_number_preds[,i], c(0.025, 0.5, 0.975))
}

caribou_number_preds_df <- cbind.data.frame(lower = caribou_number_preds_df[,1], 
                                            mean = caribou_number_preds_df[,2], upper = caribou_number_preds_df[,3], year = seq(1:23))

# Muskox number of groups model
muskox_number_m <- MCMCglmm(Groups ~ I(Year - 1987), data = muskox, family = "gaussian", pr=TRUE, nitt = 100000, burnin = 20000)
summary(muskox_number_m)
# plot(muskox_number_m$VCV)
# plot(muskox_number_m$Sol)
autocorr(muskox_number_m$VCV)

# Calculating model predictions
mean1 <- mean(muskox_number_m$Sol[,"(Intercept)"] + muskox_number_m$Sol[,"I(Year - 1987)"])
mean2 <- mean(muskox_number_m$Sol[,"(Intercept)"] + (23*muskox_number_m$Sol[,"I(Year - 1987)"]))
lower1 <- HPDinterval(muskox_number_m$Sol[,"(Intercept)"] + muskox_number_m$Sol[,"I(Year - 1987)"])[1]
upper1 <- HPDinterval(muskox_number_m$Sol[,"(Intercept)"] + muskox_number_m$Sol[,"I(Year - 1987)"])[2]
lower2 <- HPDinterval(muskox_number_m$Sol[,"(Intercept)"] + (23*muskox_number_m$Sol[,"I(Year - 1987)"]))[1] 
upper2 <- HPDinterval(muskox_number_m$Sol[,"(Intercept)"] + (23*muskox_number_m$Sol[,"I(Year - 1987)"]))[2] 

# Calculating model predictions
nyears <- 23
niter <- length(muskox_number_m$Sol[,"(Intercept)"])

muskox_number_preds <- array(NA, dim = c(niter,nyears))

for (i in 1:niter){
  for (j in 1:nyears){
    muskox_number_preds[i,j] <- muskox_number_m$Sol[i,"(Intercept)"] + muskox_number_m$Sol[i,"I(Year - 1987)"]*j
  }
}

muskox_number_preds_df <- array(NA, dim = c(nyears, 3))

for (i in 1:nyears){
  muskox_number_preds_df[i,] <- quantile(muskox_number_preds[,i], c(0.025, 0.5, 0.975))
}

muskox_number_preds_df <- cbind.data.frame(lower = muskox_number_preds_df[,1], 
                                           mean = muskox_number_preds_df[,2], upper = muskox_number_preds_df[,3], year = seq(1:23))

# Number of groups spotted, caribou/muskox, 1988-2010
(group_number <- ggplot() + 
    geom_point(data = QikHerbDF, aes(Year, Groups, colour = Species), alpha = 0.8, size = 4) +
    scale_color_manual(values = c("#ffcd44", "#1b74d3")) +
    scale_fill_manual(values = c("#ffcd44", "#1b74d3"), labels = c("Caribou", "Muskox")) +
    geom_ribbon(data = caribou_number_preds_df, aes(x = year + 1987, ymin = lower, ymax = upper), 
                fill = "#ffcd44", alpha = 0.2) +
    geom_line(data = caribou_number_preds_df, aes(x = year + 1987, y = mean), colour = "#ffcd44") +
    geom_ribbon(data = muskox_number_preds_df, aes(x = year + 1987, ymin = lower, ymax = upper), 
                fill = "#1b74d3", alpha = 0.2) +
    geom_line(data = muskox_number_preds_df, aes(x = year + 1987, y = mean), colour = "#1b74d3") +
    labs(x = " ", y = "Largest group number observed\n", title = "(d) Herbivore group number") +
    theme_QHI() +
    theme(legend.position = c(0.15, 0.9), 
          axis.line.x = element_line(color="black", size = 0.5),
          axis.line.y = element_line(color="black", size = 0.5)) +
    coord_cartesian(ylim = c(0, 155), xlim = c(1988, 2010)))

# Arranging in a panel and saving the file
SI_temp_herb <- grid.arrange(all_seasons, group_size, group_number, ncol = 2)

ggsave("Qikiqtaruk_manuscript/figures/SI_temp_herb.pdf", 
       SI_temp_herb, width = 30, height = 30, units = "cm")
