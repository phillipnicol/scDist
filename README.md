scDist: Robust identification of perturbed cell types in single-cell
RNA-seq data
================
R package version 1.1.2

## System Requirements

`R` is required to use `scDist`. In development, `R` version 4.0.0 and
greater were used, but there may be compatibility with previous
versions.

## Installation

From the R console, `devtools::install_github("phillipnicol/scDist")`.
Installation should take less than a minute on a standard machine.

## Demo

The input to `scDist` is a normalized count matrix and correpsonding
metadata that describes what condition and patient each cell belongs to.
In this demo, we create a simulated dataset with 10 cell types. The demo
should take less than a minute to run on a standard machine. The code is
also modifiable to see how `scDist` performs for different parameter
values.

``` r
library(scDist)
set.seed(1126490984)
```

Generate simulated data with 10 cell types and 5 patients in each group:

``` r
sim <- simData(nct=10,N1=5,N2=5)

dim(sim$Y) #Normalized counts
```

    ## [1] 1000 4970

``` r
rownames(sim$Y) <- 1:1000
head(sim$meta.data)
```

    ##   response patient clusters
    ## 1        1       1        a
    ## 2        1       1        a
    ## 3        1       1        a
    ## 4        1       1        a
    ## 5        1       1        a
    ## 6        1       1        a

Now we apply scDist:

``` r
out <- scDist(sim$Y,sim$meta.data,fixed.effects = "response",
              random.effects="patient",
              clusters="clusters")
```

    ## ================================================================================

The results data frame gives a summary of the estimated distance and
uncertainty for each cell type

``` r
out$results
```

    ##       Dist. 95% CI (low) 95% CI (upper)        p.val
    ## a 32.959502     32.58274       33.34466 0.0000099999
    ## b 38.491548     37.58184       39.41957 0.0000099999
    ## c 23.595509     22.59581       24.63245 0.0000099999
    ## d  8.848442      0.00000       12.49404 0.7220227798
    ## e  0.000000      0.00000        0.00000 0.9785502145
    ## f 13.480642     11.90159       15.15108 0.0000099999
    ## g  0.000000      0.00000        0.00000 0.9504704953
    ## h  0.000000      0.00000        0.00000 0.9495605044
    ## i 11.991456     10.12235       13.93573 0.0000599994
    ## j 88.566794     87.61600       89.50476 0.0000099999

The true distances are

``` r
names(sim$D.true) <- letters[1:length(sim$D.true)]
sim$D.true
```

    ##         a         b         c         d         e         f         g         h 
    ## 31.497609 36.291212 20.618938  3.602791  1.823269  7.908243  1.228762  1.030859 
    ##         i         j 
    ##  6.250907 88.944895

We can also plot the results

``` r
DistPlot(out)
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

To get a plot of genes that are associated with the perturbation use
`distGenes`:

``` r
distGenes(out, cluster = "a")
```

![](README_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

## Reference

If you use `scDist` in your work, please cite:

Nicol, P.B., Paulson, D., Qian, G., Liu, X.S., Irizarry, R.A., and Sahu,
A.D. (2024). Robust identification of perturbed cell types in
single-cell RNA-seq data. *Nature Communications*. Vol 15(7610).
<https://doi.org/10.1038/s41467-024-51649-3>.
