"Conditional inference in cis-Mendelian randomization using weak genetic factors"

by Ashish Patel, Dipender Gill, Paul Newcombe, and Stephen Burgess

load R code to use the S-LIML method: source(S-LIML.R)
load R code to use the F-LIML, F-AR, F-LM, and F-CLR methods: source(F-CLR.R)

Required data: two-sample summary data and an LD (genetic variant correlation) matrix.

bx: p-vector of genetic variant-exposure associations 
by: p-vector of genetic variant-outcome associations
sx: p-vector of standard errors of genetic variant-exposure associations 
sy: p-vector of standard errors of genetic variant-exposure associations
ld: (p x p) genetic variant correlation matrix

(the LD matrix can be from based on the same sample as one of the genetic association studies, or can be from a separate reference panel).
