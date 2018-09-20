# for generating fourier basis simulation data
par(mfrow = c(1,1))
n <- 30
k <- 9
T <- 12
basisobj <- create.fourier.basis(c(0,12), k, T) # 월별데이터를 9개의 basis로
plot(basisobj)

D <- c(runif(n*k,0,2))
coefmat <- matrix(D, k, n)

tempfd <- fd(coefmat, basisobj)
plot(tempfd) # real data는 maybe 환경(지리,기온) 가능??










# Create a minimal Fourier basis for annual data
#  using 3 basis functions
yearbasis3 <- create.fourier.basis(c(0,365),
                                   axes=list("axesIntervals") )
#  plot the basis
plot(yearbasis3)

# Identify the months with letters
plot(yearbasis3, axes=list('axesIntervals', labels=monthLetters))

# The same labels as part of the basis object
yearbasis3. <- create.fourier.basis(c(0,365),
                                    axes=list("axesIntervals", labels=monthLetters) )
plot(yearbasis3.)

# set up the Fourier basis for the monthly temperature data,
#  using 9 basis functions with period 12 months.
monthbasis <- create.fourier.basis(c(0,12), 9, 12.0)

#  plot the basis
plot(monthbasis)

# Create a false Fourier basis using 1 basis function.
falseFourierBasis <- create.fourier.basis(nbasis=1)
#  plot the basis:  constant
plot(falseFourierBasis)
