
# ICI-Kendall-Tau Manuscript

All the code for generating the manuscript:

Information-Content-Informed Kendall-tau Correlation: Utilizing Missing Values
Robert M Flight, Praneeth S Bhatt, Hunter NB Moseley

is contained within this repository.

## R Packages Needed

The repository used R 4.1.0, and {renv} 0.14.0.
For all of the other packages needed, see the file `renv.lock`.

To setup to be able to rerun everything here, you can clone the repo from github or download it from figshare, and then from within that folder:

```r
# make sure renv is installed
install.packages("renv")
# restore the packages
renv::restore()
```

## Rerun All Underlying Analyses

```r
drake::r_make()
```

## Generate Manuscript

Because of the way {drake} handles dependencies, I didn't include the manuscript generation in the {drake} run.
Therefore the manuscript needs to be done separately, all in one go.

```r
source("generate_manuscript.R")
```

## Running on HPC

The way I ran everything, was to do all the heavy work on a remote machine with 

```r
drake::r_make()
```

Then I would `rsync` the .drake directory from the remote machine to my laptop, and then run the manuscript generation.
Hopefully that will work for you as well.
