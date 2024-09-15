## Version 1.1.2 September 15, 2024

This update improves the functions in `plotting.R`, adds input checking to `scDist()`, and adds some missing documentation. See PR 8 and 10. 

## Version 1.1.1 September 4, 2024

This update primarily addresses bugs in `distGenes`.

## Version 1.1.0 July 10, 2024

This update contains new functionality for post-hoc analysis and simulated data.

**New functionality:** 

* The `distGenes` function can be used to visualize the genes associated with
the perturbation.
* The generation of simulated data is simplified with the `simData` function.
* The `scDist_tree` function performs the analysis at different resolutions
of a given clustering. 
* `condition_v_celltype_distance` compares the condition perturbation to the 
inter-cell-type distances, providing a possible way to diagnose issues in 
annotation. 

**Changes to existing functionality:**

* The `scDist` function now has a `weights` argument to allow different genes
to be weighted differently in the calculation of the distance. 

## Version 1.0.0 May 26, 2023

Initial release 
