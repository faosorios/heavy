# ID: control.R, last updated 2019/09/06, F.Osorio */

heavy.control <-
function(maxIter = 2000, tolerance = 1e-6, fix.shape = FALSE, ndraws = 500, algorithm = c("EM", "NEM"), ncycles = 5)
{
    algorithm <- match.arg(algorithm)
    choice <- switch(algorithm, "EM" = 0, "NEM" = 1)
    list(maxIter = maxIter, tolerance = tolerance, fix.shape = fix.shape, ndraws = ndraws,
         algorithm = choice, ncycles = ncycles)
}
