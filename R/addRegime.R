addRegime <- function(otree, 
                      odata, 
                      oldshifts, 
                      oldaic, 
                      oldfit, 
                      alloldaic=NULL, 
                      exclude=NULL, 
                      aic_threshold=0, 
                      verbose=FALSE, 
                      plotaic=FALSE, 
                      error_skip=FALSE, 
                      sample_shifts=FALSE, 
                      sample_threshold=2,
                      ncores = 1) {
  if (is.null(oldshifts)) {
    oldshifts<-c("1"="a")
  }
  Letters <- c(letters, paste("z",letters,sep=""),paste("zz",letters,sep="")) 
  Letters <- Letters[-which(Letters %in% oldshifts)]
  n <- otree@nterm
  nt<-dim(odata)[2]
  nodes <- otree@nodes
  uniqueshifts <- unique(oldshifts)
  k <- length(oldshifts)+1
  kk <- length(uniqueshifts)+1
  aics <- LnLs <- rep(NA, length(nodes))
  names(aics) <- names(LnLs) <- nodes
  shifts <- character(length(nodes))
  fits <- list()
  if (class(oldfit) != "list") {
    oldfit <- list(oldfit)
  }
  oldalphas <- sapply(oldfit, function(x) summary(x)$alpha)
  oldsigmas <- sapply(oldfit, function(x) summary(x)$sigma)
  oldparms <- data.frame(matrix(c(oldalphas, oldsigmas), nrow=2, byrow=T, dimnames=list(c("a","s"), names(odata))))
  odata2 <- rbind(oldparms, odata)
  if (!is.null(exclude)&!is.null(alloldaic)) {
    old <- sort(alloldaic,decreasing=TRUE)
    Nexcluded <- floor(exclude * length(old))
    excluded <- names(old[0:Nexcluded])
  }
  else {		
    excluded <- NULL	
  }	
  skip <- c(excluded, names(oldshifts))
  if (verbose) {
    print(paste("placing regime", k), quote=F)
    print(paste("testing ", length(nodes) - length(skip), "candidate models"), quote=F)
  }
  iter <- 2:length(nodes)
  iter <- iter[iter %in% skip == FALSE]
  if (ncores == 1) {
    for (i in iter) {
      Fitted <- fitHansen(i, error_skip, k, kk, n, nt, Letters, nodes, oldshifts, odata2, otree, oldaic, verbose)
      fits[[i]] <- Fitted$fit
      shifts[i] <- Fitted$letter
      names(shifts)[i] <- Fitted$node
      LnLs[i] <- Fitted$LnL
      aics[i] <- Fitted$aic
    }
  }
  else {
    registerDoParallel(ncores)
    Fitted <- foreach(i = iter, .packages = c('ouch')) %dopar% fitHansen(i, 
                                                                         error_skip, 
                                                                         k, 
                                                                         kk, 
                                                                         n, 
                                                                         nt, 
                                                                         Letters, 
                                                                         nodes, 
                                                                         oldshifts,
                                                                         odata2, 
                                                                         otree,
                                                                         oldaic, 
                                                                         verbose = FALSE)
    stopImplicitCluster()
    for (i in 1:length(Fitted)) {
      fits[[i]] <- Fitted[[i]]$fit
      shifts[i] <- Fitted[[i]]$letter
      names(shifts)[i] <- Fitted[[i]]$node
      LnLs[i] <- Fitted[[i]]$LnL
      aics[i] <- Fitted[[i]]$aic
    }
  }
  best <- names(sort(aics))[1]
  if (sample_shifts&(aics[best] - oldaic) < (aic_threshold)) {
    candidates <- aics[which((aics - min(aics, na.rm=TRUE)) <= sample_threshold & (aics-oldaic) < (aic_threshold))]
    if (verbose) {
      print(paste("sampling 1 of", length(candidates), "models within", sample_threshold, "units of best AICc"))
    }
    if (length(candidates) > 1) {
      best <- names(sample(candidates, 1))
    }
  }
  newshifts <- c(oldshifts, shifts[as.numeric(best)])
  xx<-summary(factor(newshifts))
  n_regimes <- c(k=k,kprime=kk,deltak=k-kk,c=sum(xx[xx>1]),kprime_conv=sum(xx>1),kprime_nonconv=sum(xx==1))
  if (plotaic) {
    plot(aics, ylim=c(min(oldaic, aics[best]),oldaic+20),xlim=c(0,dim(odata)[1]),main=k);abline(h=oldaic-seq(0,30,by=2));abline(h=aics[best],col="red")
  }
  if (verbose ) {
    print(paste("old AIC =",round(oldaic,2)),quote=F)
    print(paste("new AIC =",round(as.numeric(aics[best]),2)),quote=F)
    if ((aics[best]-oldaic) < (aic_threshold)) {
      print(paste("adding regime shift at node", names(aics[best])), quote=F)
    }
  }
  return(list(fit = fits[[as.numeric(best)]], all_aic = aics, aic = aics[best], savedshifts=newshifts, n_regimes=n_regimes))	
}

fitHansen <- function(i, error_skip, k, kk, n, nt, Letters, nodes, oldshifts, odata2, otree, oldaic, verbose) {
  # shifts[i] <- Letters[1]
  # names(shifts)[i] <- nodes[i]
  tempshifts <- c(oldshifts, Letters[1])
  names(tempshifts)[k] <- i
  tempregs <- repaint(otree, regshifts = tempshifts)
  if (error_skip) {
    te <- try( fit <- apply(odata2, 2, function(x) hansen(x[-c(1,2)], otree, regimes = tempregs, sqrt.alpha = sqrt(x[1]), sigma = sqrt(x[2]))) )
    if (class(te) == "try-error") {
      if (verbose) {
        print(paste("error fitting regime",i), quote=F)
      }
      LnL <- NA
      aic <- aic_threshold + 9999
    }
    else {
      LnL <- sum(sapply(fits[[i]], function(x)summary(x)$loglik))
      aic <- getAIC(LnLs[i], k+nt*(2+kk), n*nt,TRUE)
    }
  }
  else {
    fit <- apply(odata2, 2, function(x) hansen(x[-c(1,2)], otree, regimes = tempregs, sqrt.alpha = sqrt(x[1]), sigma = sqrt(x[2])))
    LnL <- sum(sapply(fit, function(x) summary(x)$loglik))
    aic <- getAIC(LnL, k+nt*(2+kk), n*nt, TRUE)
  }
  if (verbose) {
    print(c(names(aic), round(as.numeric(aic - oldaic), 2)), quote=F)
  }
  Out <- vector(mode = 'list', length = 5)
  Out[[1]] <- fit
  Out[[2]] <- LnL
  Out[[3]] <- aic
  Out[[4]] <- Letters[1]
  Out[[5]] <- nodes[i]
  names(Out) <- c('fit', 'LnL', 'aic', 'letter', 'node')
  return(Out)
}
