wsvmQP <- function(Rmat, cvec, wvec, y)  {
    n <- length(cvec)
    zsmall <- c(Rmat, cvec, 0)
    wvec <- c(wvec)
    lenz <- length(zsmall)
    sol <- .Fortran("wsvmqp",
                    xn = as.double(rep(0,n+1)),
                    n = as.integer(n),
                    wvec = as.double(wvec),
                    zsmall = as.double(zsmall),
                    y = as.double(y),
                    lenz = as.integer(lenz),
                    inform = as.integer(0)#,
                    #PACKAGE="wsvmsurf"
                    )
    if (sol$inform!=0) print("convergence warning in initialization\n")
    list(alpha=sol$xn[1:n], obj=sol$zsmall[lenz])
  }