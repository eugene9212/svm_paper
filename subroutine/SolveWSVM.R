"SolveWSVM" <- 
function(weight, lambda, x, y, kernel.function = poly.kernel, param.kernel = 1) {
   K <- kernel.function(x, x, param.kernel = param.kernel)
   Kscript <- (K * outer(y,y))
   n <- length(y)
   cvec <- -rep(1,n)
   obj <- wsvmQP(Kscript/lambda, cvec, weight, y)$alpha
obj
}