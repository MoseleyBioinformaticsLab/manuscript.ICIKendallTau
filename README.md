
# ICI-Kendall-Tau Manuscript

[![DOI](https://zenodo.org/badge/429873882.svg)](https://zenodo.org/badge/latestdoi/429873882)

All the code for generating the manuscript:

Information-Content-Informed Kendall-tau Correlation: Utilizing Missing Values
Robert M Flight, Praneeth S Bhatt, Hunter NB Moseley
[biorxiv](https://doi.org/10.1101/2022.02.24.481854)

is contained within this repository.

## R Packages Needed

The repository used R 4.1.0, and {renv} 1.0.0
For all of the other packages needed, see the file `renv.lock`.

To setup to be able to rerun everything here, you can clone the repo from github or download it from Zenodo, and then from within that folder:

```r
# make sure renv is installed
install.packages("renv")
# restore the packages
renv::restore()
```

## Rerun All Underlying Analyses

```r
targets::tar_make()
```


## Obtaining targets

To make it easier to at least generate the manuscript, there is a copy of the {drake} cache on Zenodo at [here]().
You can download it using `wget`.

```
# make sure you are wherever you want the drake cache, like wherever
# you cloned the github repo to, at the top level of the directory.
wget https://zenodo.org/record/6310898/files/manuscript.ICIKendallTau.drake_cache.tgz?download=1
mv manuscript.ICIKendallTau.drake_cache.tgz?download=1 manuscript.ICIKendallTau.drake_cache.tgz
tar -xzf manuscript.ICIKendallTau.targets_cache.tgz
```
