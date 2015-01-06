ews_fli1_and_e2f_nar_article
============================

This repository provides the MATLAB code used for doing the statistical analysis of the research article
"EWS-FLI1 employs an E2F switch to drive target gene expression", 
R. Schwentner, T. Papamarkou2, M.O. Kauer, V. Stathopoulos, F. Yang, Sv. Bilke, P.S. Meltzer, M. Girolami, H. Kovar,
*Nucleic Acids Research*
(to appear soon).

The statistical analysis employs Bayesian model selection in order to identify the most suitable ODE system for the
the regulation of E2F3 target genes among a pool of four candidate mechanistic models. From a technical point of view,
population MCMC has been used in conjunction with the simplified manifold Metropolis-adjusted Langevin algorithm
(SMMALA), see the article for details. The population MCMC scheme was run in parallel on a cluster using the MATLAB
inteface to MPI.

It is noted that this public code has not been polished in the form of a MATLAB package. However, it is provided as a
guide to the interested researcher who would wish to replicate the statistical analysis after some configuration of
the provided repository or to recode the MCMC simulations in a programming language of their preference by using the
existing code as a guide.

To bring the code to an operational state, the user would need to set up some absolute paths pointing to the data and
the files holding the specified models. Furthermore, the following dependencies need to be installed on the cluster:

* cln-1.3.2
* ginac-1.6.2
* mxml-2.6
* sundials-2.5.0
* vfgen-2.4.1 (use the patched vfgen folder in dependencies/vfgen-2.4.1)
