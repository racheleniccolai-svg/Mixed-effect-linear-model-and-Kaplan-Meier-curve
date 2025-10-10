# Mixed effect linear model and Kaplan Meier curve
Tumor-growth data were analyzed using a linear mixed-effects model (http://dx.doi.org/10.1007/b98882) fitted by Restricted Maximum Likelihood (REML) to evaluate the treatment effects over time. Tumor volumes were normalized relative to baseline measurements to correct for initial tumor-size differences and then log-transformed to model the natural exponential growth of tumor. The model incorporated fixed effects for the time, the interaction between time and treatment/cell line as well as the interaction between the three to assess differential temporal dynamics between groups. Individual animal variability was accounted by adding a random intercept for each mouse. The significance level was set at 0.05. All statistical analyses were performed using R (version 4.3.2), using the nlme package (version 3.1-164) to fit the mixed-effect model. After the models were fitted, the slopes of the tumor growth curves and their respectives confidence intervals were computed using the emmeans package (version 1Lenth R (2025). emmeans: Estimated Marginal Means, aka Least-Squares Means. R package version 1.11.1-00001.10.2, ).

## Installation / Requirements
```
library(dplyr)
library(emmeans)
library(forcats)
library(ggplot2)
library(multcomp)
library(nlme)
library(readxl)
library(stringr)
library(survival)
library(survminer)
library(tidyr)
```

# Script
For more information about the script open the file "In vivo statistical analysis"
