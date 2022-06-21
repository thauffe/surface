### Parallelized version of of the R package surface

The main functions fit a series of Hansen models, then identify cases of convergent evolution where multiple lineages have shifted to the same adaptive peak. It enables to uses multiple CPUs through the *doParallel* package and is faster than the original version of surface. For more information on surface see [Ingram and Mahler (2013)](http://www.sciencemag.org/content/341/6143/292).

##### Usage
The package can be directly installed from this github repository.

```{r, warning = F, echo = F}
library(remotes)
install_github("thauffe/surface")
library(simDES)
```


The argument *ncores* allows to specify the number of CPU cores.

```{r, warning = F, echo = F}
data(surfaceDemo)
Fit <- runSurface(tree = surfaceDemo$tree, dat = surfaceDemo$sim$dat, verbose = TRUE, ncores = 2)
```


