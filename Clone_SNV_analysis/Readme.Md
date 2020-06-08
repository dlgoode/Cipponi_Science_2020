R code for analysis of subclonal SNVs from WGS data, including filtering and tree building steps.

collect_private_mutations_each_sample.r	- create a list of private variants for each sample

identify subclonal SNVs.r	- identify sub-clonal variants according to sample-specific cutoff

make_clone_identity_matrix.r	- calculate fraction of shared SNVs between samples, for tree building
