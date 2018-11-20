# Practice on functional data example #

# load packages
library(fda); library(fda.usc)
# load data
data(tecator)
x=tecator$absorp.fdata
str(x)
y=tecator$y$Fat
str(y)

# data trimming(1)
tt=x[["argvals"]] # extract time range from x(fdata)
dataf=as.data.frame(tecator$y) # transform y as data.frame structure
nbasis.x=11 # create basis used for fdata or fd covariates.
nbasis.b=7  # create basis used for beta parameter estimation.

# create basis
basis1=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.x)
basis2=create.bspline.basis(rangeval=range(tt),nbasis=nbasis.b)

# formula
f=Fat~Protein+x 
# y : Fat, Water, Protein
# x : fdata of tecator$absorp

# data trimming(2) : basis
basis.x=list("x"=basis1)
basis.b=list("x"=basis2)
ldata=list("df"=dataf,"x"=x) 
# dataf : no fdata
# x : basis1 and basis2

# Fitting (train)
res=fregre.glm(f,family=gaussian(),data=ldata,basis.x=basis.x,basis.b=basis.b)

# summary
summary(res)