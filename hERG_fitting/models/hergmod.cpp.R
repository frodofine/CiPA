`.sourceCpp_1_DLLInfo` <- dyn.load('models/hergmod.dll')

derivs <- Rcpp:::sourceCppFunction(function(t, y, parms) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_derivs')

rm(`.sourceCpp_1_DLLInfo`)
