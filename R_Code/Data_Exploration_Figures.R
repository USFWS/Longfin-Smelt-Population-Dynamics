#### 1. INTRODUCTION ####
# 2024-03-14
# Vanessa Tobias <vanessa_tobias@fws.gov>

# The purpose of this code is to document the population dynamics model
#    appendix/technical note, as it appears in the Longfin Smelt Species
#    Status Assessment.

# This code is part of LFS Technical Note 4 (v. 3.2).
# It creates figures 2, 3, and 4, which are the data exploration graphs.

# The code file called "PopModel_BH_LN.R" fits the population model
#   and creates the later figures in the technical note.


#### 2. SETUP ####
## Packages ####
# This code does not require any additional packages


## Functions ####
## Residual Plotting Function
plotXresids <- function(dataset, model, x, xlab, main, loess = FALSE){
  # PLOTS a variable against residuals from a model that you specify
  # --- dataset = name of the dataset you want to plot from. No quotes.
  # --- model = name of the model that you want to extract residuals from
  # --- x = the variable name (in quotes) that you want to plot on the x axis
  # --- xlab = the label for the x axis (in quotes)
  # --- main = label for the top of the plot; should describe the model
  
  # rename the dataset for convenience
  dat <- dataset
  # extract the variable names from the model
  z <- all.vars(formula(model))[2]
  y <- all.vars(formula(model))[1]
  # remove any NAs from the dataset (to make the vectors the same length for plotting)
  dat <- dat[which(!is.na(eval(parse(text = paste0("dat$", y)))) & !is.na(eval(parse(text = paste0("dat$", z))))),]
  # draw the residual plot
  plot(eval(parse(text = paste0("dat$", x))),
       residuals(model),
       pch = 16, col = "black", bty = "l",
       xlab = xlab,
       ylab = "Residual",
       main = main)
  # add a horizontal line at residual = 0
  abline(h = 0, lty = 2, col = "grey")
  # add a smooth line, if requested
  if (loess == TRUE) lines(loess.smooth(eval(parse(text = paste0("dat$", x))),
                                        residuals(model), family = "gaussian"),
                           lwd = 2, col = "grey")
}


## Data ####
# Get Bay Study data from a CSV file
sfbs <- read.csv("./Data_Original/SFBS_just_LFS_Indices.csv", header = TRUE)
sfbs$log.ot.0 <- log(sfbs$ot.0)
sfbs$log.ot.1 <- log(sfbs$ot.1)
sfbs$log.ot.2 <- log(sfbs$ot.2)
sfbs$log.ot.2[which(is.infinite(sfbs$log.ot.2))] <- NA

sfbs$log.mwt.0 <- log(sfbs$mwt.0)
sfbs$log.mwt.1 <- log(sfbs$mwt.1)
sfbs$log.mwt.2 <- log(sfbs$mwt.2)
sfbs$log.mwt.2[which(is.infinite(sfbs$log.mwt.2))] <- NA
# Reorganize Bay Study data into year classes

yrClass <- data.frame(brood.year = sfbs$year,
                      ot.0 = sfbs$ot.0,
                      ot.1 = c(sfbs$ot.1[-1], NA),
                      ot.2 = c(sfbs$ot.2[-(1:2)], NA, NA),
                      mwt.0 = sfbs$mwt.0,
                      mwt.1 = c(sfbs$mwt.1[-1], NA),
                      mwt.2 = c(sfbs$mwt.2[-(1:2)], NA, NA))

# Calculate log of each index
yrClass$log.ot.0 <- log(yrClass$ot.0)
yrClass$log.ot.1 <- log(yrClass$ot.1)
yrClass$log.ot.2 <- log(yrClass$ot.2)
yrClass$log.ot.2[which(is.infinite(yrClass$log.ot.2))] <- NA
yrClass$log.mwt.0 <- log(yrClass$mwt.0)
yrClass$log.mwt.1 <- log(yrClass$mwt.1)
yrClass$log.mwt.2 <- log(yrClass$mwt.2)
yrClass$log.mwt.2[which(is.infinite(yrClass$log.mwt.2))] <- NA

# calculate survival ratios
yrClass$log.ot.01 <- yrClass$log.ot.1/yrClass$log.ot.0
yrClass$log.ot.02 <- yrClass$log.ot.2/yrClass$log.ot.0
yrClass$log.ot.12 <- yrClass$log.ot.2/yrClass$log.ot.1
yrClass$log.mwt.01 <- yrClass$log.mwt.1/yrClass$log.mwt.0
yrClass$log.mwt.02 <- yrClass$log.mwt.2/yrClass$log.mwt.0
yrClass$log.mwt.12 <- yrClass$log.mwt.2/yrClass$log.mwt.1

#### 3. ANALYSIS ####
## Correlations ####
### Otter Trawl ####
corrAge0_1.ot <- lm(log.ot.1 ~ log.ot.0, data = yrClass)
corrAge1_2.ot <- lm(log.ot.2 ~ log.ot.1, data = yrClass)
corrAge2_0.ot <- lm(log.ot.0 ~ log.ot.2, data = sfbs,
                    na.action = na.omit)

### Midwater Trawl ####
corrAge0_1.mwt <- lm(log.mwt.1 ~ log.mwt.0 - 1, data = yrClass)
corrAge1_2.mwt <- lm(log.mwt.2 ~ log.mwt.1 - 1, data = yrClass)
corrAge2_0.mwt <- lm(log.mwt.0 ~ log.mwt.2 -1, data = sfbs,
                     na.action = na.omit)

### Otter Trawl v. Midwater Trawl ####
corrAge0gears <- lm(log.mwt.0 ~ log.ot.0, data = yrClass)
corrAge1gears <- lm(log.mwt.1 ~ log.ot.1, data = yrClass)
corrAge2gears <- lm(log.mwt.2 ~ log.ot.2, data = yrClass)


#### 4. PLOTS ####

### Fig. 2 Otter Trawl ####
summary(corrAge0_1.ot)
summary(corrAge1_2.ot)
summary(corrAge2_0.ot)

png("./Figures/Fig2_OT.png",
    height = 10.125, width = 9.45, units = "in", res = 300)
par(mfrow = c(3, 2)) # create a 3x2 grid, fill by rows

# Row 1: age 0 v. age 1
plot(yrClass$log.ot.0, yrClass$log.ot.1,
     pch = 16, col = "black", bty = "l",
     xlab = "log(Age0)",
     ylab = "log(Age1)",
     main = "Otter Trawl")
abline(coefficients(corrAge0_1.ot))

plotXresids(dataset = yrClass, model = corrAge0_1.ot,
            x = "brood.year", xlab = "Brood Year",
            main = "Otter Trawl Age 0 - Age 1",
            loess = TRUE)

# Row 2: age 1 v. age 2
plot(yrClass$log.ot.1, yrClass$log.ot.2,
     pch = 16, col = "black", bty = "l",
     xlab = "log(Age1)",
     ylab = "log(Age2)",
     main = "Otter Trawl")
abline(coefficients(corrAge1_2.ot))

plotXresids(dataset = yrClass, model = corrAge1_2.ot,
            x = "brood.year", xlab = "Brood Year",
            main = "Otter Trawl Age 1 - Age 2",
            loess = TRUE)

# Row 3: age 2 v. age 0
plot(yrClass$log.ot.2, yrClass$log.ot.0,
     pch = 16, col = "black", bty = "l",
     xlab = "log(Age2)",
     ylab = "log(Age0)",
     main = "Otter Trawl")
abline(coefficients(corrAge2_0.ot))

plotXresids(dataset = sfbs,
            model = corrAge2_0.ot, 
            x = "year", xlab = "Brood Year",
            main = "Otter Trawl Age 2 - Age 0",
            loess = TRUE)
dev.off()

### Fig. 3 Midwater Trawl ####
summary(corrAge0_1.mwt)
summary(corrAge1_2.mwt)
summary(corrAge2_0.mwt)

png("./Figures/Fig3_MWT.png",
    height = 10.125, width = 9.45, units = "in", res = 300)
par(mfrow = c(3, 2)) # create a 3x2 grid, fill by rows

# Row 1: age 0 v. age 1
plot(yrClass$log.mwt.0, yrClass$log.mwt.1,
     pch = 16, col = "black", bty = "l",
     xlab = "log(Age0)",
     ylab = "log(Age1)",
     main = "Midwater Trawl")
abline(b=coefficients(corrAge0_1.mwt), a = 0)

plotXresids(dataset = yrClass, model = corrAge0_1.mwt, 
            x = "brood.year", xlab = "Brood Year",
            main = "Midwater Trawl Age 0 - Age 1",
            loess = TRUE)

# Row 2: age 1 v. age 2
plot(yrClass$log.mwt.1, yrClass$log.mwt.2,
     pch = 16, col = "black", bty = "l",
     xlab = "log(Age1)",
     ylab = "log(Age2)",
     main = "Midwater Trawl")
abline(b = coefficients(corrAge1_2.mwt), a = 0)

plotXresids(dataset = yrClass, model = corrAge1_2.mwt, 
            x = "brood.year", xlab = "Brood Year",
            main = "Midwater Trawl Age 1 - Age 2",
            loess = TRUE)

# Row 3: age 2 v. age 0
plot(yrClass$log.mwt.2, yrClass$log.mwt.0,
     pch = 16, col = "black", bty = "l",
     xlab = "log(Age2)",
     ylab = "log(Age0)",
     main = "Midwater Trawl")
abline(b = coefficients(corrAge2_0.mwt), a = 0)

plotXresids(dataset = sfbs, model = corrAge2_0.mwt, 
            x = "year", xlab = "Brood Year",
            main = "Midwater Trawl Age 2 - Age 0",
            loess = TRUE)

dev.off()


### Fig. 4 Gear Comparison ####
summary(corrAge0gears)
summary(corrAge1gears)
summary(corrAge2gears)

png("./Figures/Fig4_gear_comp.png",
    height = 10.125, width = 9.45, units = "in", res = 300)
par(mfrow = c(3, 2)) # create a 3x2 grid, fill by rows

# Row 1: Age-0
plot(yrClass$log.ot.0, yrClass$log.mwt.0,
     pch = 16, col = "black", bty = "l",
     ylab = "log(Midwater Trawl)",
     xlab = "log(Otter Trawl)",
     main = "Age-0")
abline(coefficients(corrAge0gears))

plotXresids(dataset = yrClass, model = corrAge0gears, 
            x = "brood.year", xlab = "Brood Year",
            main = "Age-0 OT v.MWT",
            loess = TRUE)

# Row 2: Age-1
plot(yrClass$log.ot.1, yrClass$log.mwt.1,
     pch = 16, col = "black", bty = "l",
     ylab = "log(Midwater Trawl)",
     xlab = "log(Otter Trawl)",
     main = "Age-1")
abline(coefficients(corrAge1gears))

plotXresids(dataset = yrClass, model = corrAge1gears, 
            x = "brood.year", xlab = "Brood Year",
            main = "Age-1 OT v.MWT",
            loess = TRUE)

# Row 3: Age-2
plot(yrClass$log.ot.2, yrClass$log.mwt.2,
     pch = 16, col = "black", bty = "l",
     ylab = "log(Midwater Trawl)",
     xlab = "log(Otter Trawl)",
     main = "Age-2")
abline(coefficients(corrAge2gears))

plotXresids(dataset = yrClass, model = corrAge2gears, 
            x = "brood.year", xlab = "Brood Year",
            main = "Age-2 OT v.MWT",
            loess = TRUE)

dev.off()