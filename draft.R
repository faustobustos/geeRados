library(roxygen2)
library(devtools)
library(geepack)
library(boot)
setwd('C:/Users/cjcar/Documents/Github/geeRados')

gee.boot <- function(formula, df, varname, n.boot,
                     fam='binomial', scalefix=TRUE) {

  gee_no2 <- geeglm(formula, family=fam, data=df, id=as.expression(varname), 
                    corstr = "exchangeable", std.err="san.se", scale.fix = scalefix)
  summary(gee_no2)
  
}

gee.boot(everinf ~ age, practice, varname = practice$homeid_cohort)

