# BBDG
BBDG_fun.R are R functions for Beta Binomial modeling with dynamic overdispersion rates on base pairs and its based differential expression analysis.
BBDG_MLE.o and BBDG_MLE.so are C object code and the dll.

------------------------------------

First, parameters (D, beta, gamma, pi, pn)
will be estimated and second, DE testing will be performed.

------------------------------------

#estimate D, beta and gamma

D, beta coefficients of local sequence around a particular
nucleotide of a particular gene and the overdispersion parameter gamma will be
estimated using the least-squares estimation method. The implement of
least-squares estimation is ordinary and is not described here.

------------------------------------

Next, the proportion of reads counts from two samples, pi, will be
set as the neutral proportion pn and pairwise D will be estimated by maximizing
the beta-binomial log likelihood.

The data matrix, “dat”, has four columns. The first column can be
the position indicator of current nucleotide, the second column is the local
sequence effect calculated based on estimated beta, and columns three and four
store the counts of reads beginning from the current nucleotide for two samples
separately.

The strand-specific D will be initiated by the least-squares
estimation in “a0”, which is a vector with two values for sense and antisense strands.

The fixed parameters, pi and gamma, will be set in “afix”, which is
a vector with three values (pi, gamma, 0).

“D.index” is a vector of indicators of strands, which is 1 for
sense and 2 for antisense.

#estimate D

D <- bb.mle (dat, which="D_SAS", a0=a0, afix=afix, d.index=D.index)

------------------------------------

Further, D, beta (for the full model), and gamma will be set as known
parameters and pi will be updated by maximizing the beta-binomial log
likelihood. The parameters will be submitted in “afix”, which is a vector with
three values (D, gamma, 0). For the full model estimation, beta should be
provided to “coe”.

#estimate pi

Pi <- bb.mle (dat, which="mu",  coe=NULL,  a0=afix)

------------------------------------

Given estimated pi, pn, D, beta (for the full model), and gamma,
the likelihood of beta binomial model will be calculated and the likelihood
ratio test will be performed to asset the significance of differential
expression. Those estimated parameters will be set in “afix”, which is a vector
with three values (D, gamma, 0). For the full model estimation, beta should be
provided to “coe”.

#DE test

fval.bb.pi <- .opl.u.i.c (pi, coe=NULL, afix, dat) * (-1)

fval.bb.pn <- .opl.u.i.c (pn, coe=NULL, afix, dat) * (-1)

p.b <- .bb.chis.test (fval.bb.pi, fval.bb.pn)

