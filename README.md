# Longfin Smelt Population Dynamics Model

## Overview

This repository contains data and code to support the development of a population dynamics model for the Bay-Delta DPS of [Longfin Smelt (*Spirinchus thaleichthys*)](https://www.fws.gov/species/longfin-smelt-spirinchus-thaleichthys). This model uses indices of Longfin Smelt abundance for three age-classes to estimate survival and reproduction parameters. This code was used in Appendix E of the [Species Status Assessment (SSA)](https://www.fws.gov/node/4531791) for the Bay-Delta DPS of Longfin Smelt. Further development of this project may also occur.

There are two main R scripts in this repository. The first (Data_Exploration_Figures.R) creates graphs that are used to visualize the relationships among age-classes and gear types. The second (PopModel_BH_LN.R) fits a state space model to estimate age-0 survival, age-1 survival, and two parameters for reproduction based on the Beverton-Holt model for density dependent populationg growth.

Data for this study were derived from the Interagency Ecological Program's [long-term monitoring data](https://iep.ca.gov/Data/IEP-Survey-Data). Specifically, this model uses data and abundance indices from the San Francisco Bay Study program.


## Installation

The code contained in this repository produced the figures and derived data that are also contained in this repository using R version 4.3.1 ("Beagle Scouts").
This repository is set up to work as an RStudio project; as such, file paths contained in the code are relative paths. This code will also work as a set of stand-alone R scripts, but the user may need to specify the locations of data files in a more explicity way.

## Getting help

Contact the project maintainer for help with this code.  

## Contribute

Contact the [project maintainer](vanessa_tobias@fws.gov) for information about contributing to this project. Submit a [GitHub Issue](https://github.com/USFWS/Longfin-Smelt-Length-at-Date/issues) to report a bug or request a feature or enhancement.

