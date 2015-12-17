pkgname <- "ars"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('ars')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("ars")
### * ars

flush(stderr()); flush(stdout())

### Name: ars
### Title: Adaptive rejection sampler The main ars() function
### Aliases: ars

### ** Examples

samples <- ars( function(x) log( dnorm(x, 0, 1) ), 10000, c(-Inf, Inf))



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
