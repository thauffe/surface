### Parallelized version of the R package surface

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
Fit <- runSurface(tree = surfaceDemo$tree, dat = surfaceDemo$sim$dat, ncores = 2)
```

The time needed to identify cases of convergent evolution decreases with the number of CPU cores.

```{r, warning = F, echo = F}
Time <- matrix(NA_real_, ncol = 4, nrow = 10)
for (i in 1:4) {
  for (y in 1:10) {
    Tmp <- system.time(Fit <- runSurface(tree = surfaceDemo$tree,
                                         dat = surfaceDemo$sim$dat,
                                         ncores = i))
    Time[y, i] <- Tmp[3]
  }
}
par(mar = c(4, 4, 0.1, 0.1), las = 1)
plot(1:4, colMeans(Time), type = "l", ylim = range(Time),
     xlab = "CPU cores (N)", ylab = "Time (s)", xaxt = "n")
axis(side = 1, 1:4)
points(rep(1:4, each = nrow(Time)), c(Time), 
       pch = 19, col = adjustcolor("black", alpha = 0.3))
```

![wall-time](https://github.com/thauffe/surface/blob/master/walltime.png)
