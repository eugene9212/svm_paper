# Practice on functional data example #
rm(list = ls())

# load packages
library(fda); library(fda.usc)

# load data
data(tecator)
x=tecator$absorp.fdata
y.Fat=tecator$y$Fat

# -------------------------------- For Gaussian ---------------------------------------- #
# data trimming(1)
tt=x[["argvals"]] # extract time range from x(fdata)
dataf=as.data.frame(tecator$y) # transform y as data.frame structure
nbasis.x=11 # create basis used for fdata or fd covariates.
nbasis.b=4  # create basis used for beta parameter estimation.

# formula
f=Fat~Protein+x 

# create basis
basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)

# Fitting (train) _ gaussian
basis.x=list("x"=basis1) # has to be the same name
basis.b=list("x"=basis2)
ldata=list("df"=dataf,"x"=x)
res=fregre.glm(f,family=gaussian(), data=ldata, basis.x=basis.x, basis.b=basis.b)
summary(res)
fres$basis.b
res$beta.l
?fregre.glm

# summary
summary(res)

# -------------------------------- For Binomial ---------------------------------------- #
# data trimming(1)
tt=x[["argvals"]] # extract time range from x(fdata)

## Change y into 1 or 0
y1 <- as.vector(y.Fat)
y1[which(y1 < 13)] = 0
y1[which(y1 >= 13)] = 1

## formula
f=Fat~Protein+x 
# y : Fat, Water, Protein
# x : fdata of tecator$absorp

# Create basis
nbasis.x=11 # create basis used for fdata or fd covariates.
nbasis.b=10  # create basis used for beta parameter estimation.

basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)

# Fitting (train) _ gaussian
basis.x=list("x"=basis1) # has to be the same name
basis.b=list("x"=basis2)

# Appropriate data form
dataf=as.data.frame(tecator$y)
dataf$Fat = y1
ldata=list("df"=dataf,"x"=x)

# Fitting (train) _ binomial
res=fregre.glm(f,binomial(link = "logit"),data=ldata, basis.x=basis.x, basis.b=basis.b)
summary(res)
