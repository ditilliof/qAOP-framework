# Phase plane analysis in R: grind.R 
# Rob de Boer, Utrecht University
grind_version <- "21-01-2024"

library(deSolve)     # run() calls the ode() function
library(rootSolve)   # newton() and continue() call steady()
library(coda)        # required by FME
library(FME)         # fit() calls modFit() and modCost()

options(stringsAsFactors=FALSE)
colors <- c("red","blue","darkgreen","darkorange","darkmagenta","gold","darkorchid","aquamarine","deeppink","gray",seq(2,991))
#colors <- seq(2,101)  # Use standard R colors
ncolors <- length(colors)
sizeLegend <- 0.1    # legend size is 75% of the default value in R
font.main <- 1        # plain (default is bold: 2)
font.sub  <- 1        # plain (default is bold: 2)

x_plane <- 1; xmin_plane <- 0; xmax_plane <- 1.1
y_plane <- 2; ymin_plane <- 0; ymax_plane <- 1.1
log_plane <- ""; eps_plane <- 0; addone_plane <- FALSE

plane <- function(xmin=0, xmax=1.1, ymin=0, ymax=1.1, xlab="", ylab="", log="", npixels=500, state=inistate, parms=pars, odes=feedbackqAOP, x=1, y=2, time=0, grid=5, eps=0, show=NULL, addone=FALSE, portrait=FALSE, vector=FALSE, add=FALSE, legend=TRUE, zero=TRUE, lwd=2, col="black", pch=20, ...) {
  # Make a phase plane with nullclines and/or phase portrait
  dots <- list(...)
  if (!is.null(dots)) {
    unknown <- names(dots[!names(dots) %in% c(args_run,args_plot)])
    if (length(unknown)>0) warning(paste("Unknown argument(s):",unknown,sep=" "))
  }
  if (!is.null(dots)) dots_run <- dots[names(dots) %in% args_run]
  else dots_run <- NULL
  if (add) {
    x <- x_plane
    y <- y_plane
    xmin <- xmin_plane; xmax <- xmax_plane
    ymin <- ymin_plane; ymax <- ymax_plane
    log <- log_plane; eps <- eps_plane; addone <- addone_plane
  } else {
    if (!is.numeric(x)) x <- index(x,names(parms))
    if (!is.numeric(y)) y <- index(y,names(state))
    x_plane <<- x
    y_plane <<- y
    xmin_plane <<- xmin; xmax_plane <<- xmax
    ymin_plane <<- ymin; ymax_plane <<- ymax 
    log_plane <<- log; eps_plane <<- eps; addone_plane <<- addone
  }
  if (!is.null(show)) ishows <- index(show, names(state))
  else ishows <- c(x, y)
  nvar <- length(state)
  if (zero) state[1:nvar] <- rep(0,nvar)
  lvec <- 50                         # length of vector
  logx <- ifelse(grepl('x',log), TRUE, FALSE)
  logy <- ifelse(grepl('y',log), TRUE, FALSE)
  
  if (logx) xc <- 10^seq(log10(xmin),log10(xmax),length.out=npixels)
  else {
    if (eps != 0 & !logx & xmin == 0) xmin <- xmin + eps
    xc <- seq(xmin,xmax,length.out=npixels)
  }
  if (logy) yc <- 10^seq(log10(ymin),log10(ymax),length.out=npixels)
  else {
    if (eps != 0 & !logy & ymin == 0) ymin <- ymin + eps
    yc <- seq(ymin,ymax,length.out=npixels)
  }
  if (xlab == "") xlab <- names(state)[x]
  if (ylab == "") ylab <- names(state)[y]
  if (addone) {
    if (logx) xlab <- paste(xlab,"+ 1")
    if (logy) ylab <- paste(ylab,"+ 1")
  }
  if (!add) {
    do.call('plot',c(list(1,1,type='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=xlab,ylab=ylab,log=log,font.main=font.main,font.sub=font.sub),dots[names(dots) %in% args_plot]))
    if (legend)
      legend("topright",legend=names(state)[ishows],col=colors[ishows],lty=1,lwd=lwd,cex=sizeLegend)
  }
  
  npixels2 <- npixels^2
  vstate <- as.list(state)
  vparms <- as.list(parms)
  vparms <- lapply(vparms,rep,vparms,npixels2)
  vstate <- lapply(vstate,rep,vstate,npixels2)
  #for (j in seq(1,nvar)) if (j!=x & j!=y) vstate[[j]]<-rep.int(vstate[[j]],npixels2)
  vstate[[x]] <- rep.int(xc, npixels)
  vstate[[y]] <- rep.int(yc, rep.int(npixels, npixels))
  if (addone & logx) vstate[[x]] <- vstate[[x]] - 1
  if (addone & logy) vstate[[y]] <- vstate[[y]] - 1
  dvstate <- odes(time,vstate,vparms)[[1]]
  dim(dvstate) <- c(npixels,npixels,nvar)
  for (i in ishows) 
    contour(xc,yc,dvstate[,,i],levels=0,drawlabels=FALSE,add=TRUE,col=colors[i],lwd=lwd)
  
  if (portrait | vector) {
    if (logx) {dx <- (log10(xmax)-log10(xmin))/grid; vx <- 1+3.32*grid*dx/lvec}
    else {dx <- (xmax-xmin)/grid; vx = grid*dx/lvec}
    if (logy) {dy <- (log10(ymax)-log10(ymin))/grid; vy <- 1+3.32*grid*dy/lvec}
    else {dy <- (ymax-ymin)/grid; vy = grid*dy/lvec}
    
    for (i in seq(1,grid)) {
      if (logx) state[x] <- 10^((i-1)*dx + dx/2 + log10(xmin))
      else state[x] <- (i-1)*dx + dx/2 + xmin
      for (j in seq(1,grid,1)) {
        if (logy) state[y] <- 10^((j-1)*dy + dy/2 + log10(ymin))
        else state[y] <- (j-1)*dy + dy/2 + ymin      
        if (portrait) {
          points(state[x],state[y],pch=pch)
          nsol <- do.call('run',c(list(state=state,parms=parms,odes=odes,timeplot=FALSE,table=TRUE),dots_run))
          lines(cbind(nsol[x+1],nsol[y+1]),col=col)
        }
        if (vector) {
          dt <- sign(odes(time,state,parms)[[1]])
          if (logx) lines(c(state[x],state[x]*vx^dt[x]), c(state[y],state[y]))
          else lines(c(state[x],state[x]+vx*dt[x]), c(state[y],state[y]))
          if (logy) lines(c(state[x],state[x]), c(state[y],state[y]*vy^dt[y]))
          else lines(c(state[x],state[x]), c(state[y],state[y]+vy*dt[y]))
        }
      }
    }
  }
}

dummyEvent <-function(t, state , parms) return(state)

run <- function(tmax=100, tstep=1, state=inistate, parms=pars, odes=feedbackqAOP, ymin=0, ymax=NULL, log="", xlab="Time", ylab="Density", tmin=0, draw=lines, times=NULL, show=NULL, arrest=NULL, events=NULL, after=NULL, tweak=NULL, timeplot=TRUE, traject=FALSE, table=FALSE, add=FALSE, legend=TRUE, solution=FALSE, delay=FALSE, lwd=2, col="black", pch=20, ...) {   
  # run model and make a table, time plot, or trajectory
  if (delay & (solution | !is.null(after))) stop("Don't use solution or after with delay equations")
  if (delay) args_run <- args_run_dde
  dots <- list(...)
  if (!is.null(dots)) {
    unknown <- names(dots[!names(dots) %in% c(args_run,args_plot)])
    if (length(unknown)>0) warning(paste("Unknown argument(s):",unknown,sep=" "))
    dots_run <- dots[names(dots) %in% args_run]
  }else dots_run <- NULL
  nvar <- length(state)
  #if (!is.numeric(x)) x <- index(x,names(state))
  #if (!is.numeric(y)) y <- index(y,names(state))
  #if (!is.null(show)) ishows <- index(show, names(state))
  #else ishows <- seq(1,nvar)
  if (is.null(times)) times <- seq(tmin,tmax,by=tstep)
  else { 
    times <- sort(times)
    tmin <- min(times)
    tmax <- max(times)
  }
  if (!is.null(arrest)) {
    if (!is.null(events)) stop("Don't combine the option arrest with events")
    if (!is.numeric(arrest)) arrest <- sort(as.numeric(parms[arrest]))
    nearby <- nearestEvent(arrest,times)  # Find nearby times
    nearby <- nearby[nearby < arrest]     # Pick those before arrest
    lennear <- length(nearby)
    if (lennear == 1 && nearby[1] == 0) lennear <- 0
    if (lennear > 0) {                    # Add to arrest for safety
      if (nearby[1] == 0) nearby <- nearby[2:lennear]
      arrest <- sort(unique(c(nearby,arrest)))   
    }
    events <- list(func=dummyEvent,time=arrest)
    times <- cleanEventTimes(times,arrest,eps = .Machine$double.eps*10)
    #times <- cleanEventTimes(times,arrest,eps=min(arrest)/10000)
    times <- sort(c(times,arrest))
  }
  if (solution) {                              # Run in one go
    if (!is.null(after)) stop("Don't combine the option after with solution")
    nsol <- sapply(times,odes,state,parms)
    if (is.list(nsol)) {
      nsol <- unlist(nsol)
      if (nvar > 1) dim(nsol) <- c(nvar,length(times))
    }
    if (nvar > 1) nsol <- data.frame(times,t(nsol))
    else nsol <- data.frame(times,nsol)
    names(nsol) <- c("time",names(state))
  }else{
    if (is.null(after)){                       # Run in one go 
      if (!delay) nsol <- as.data.frame(do.call('ode',c(list(times=times,func=odes,y=state,parms=parms,events=events),dots_run)))
      else nsol <- as.data.frame(do.call('dede',c(list(times=times,func=odes,y=state,parms=parms,events=events),dots_run)))
    }else{                                    # After:Run many in steps
      keep <- state
      nsol <- t(sapply(seq(length(times)-1),function(i){
        t <- times[i+1]
        f <- do.call('ode',c(list(times=c(times[i],t),func=odes,y=state,parms=parms),dots_run))
        dim(f) <- c(2,nvar+1)
        state[1:nvar] <- f[2,2:(nvar+1)]
        eval(parse(text=after))
        parms <<- parms
        state <<- state
      }
      ))
      if (nvar > 1) nsol <- as.data.frame(cbind(times,rbind(as.numeric(keep),nsol)))
      else nsol <- as.data.frame(cbind(times,c(keep,nsol)))
      names(nsol) <- c("time",names(state))
      state <- keep
    }
  }
  if (!is.null(tweak)) eval(parse(text=tweak))
  if (timeplot & !traject)
    do.call('timePlot',c(list(data=nsol,tmin=tmin,tmax=tmax,ymin=ymin,ymax=ymax,log=log,add=add,xlab=xlab,ylab=ylab,show=show,draw=draw,lwd=lwd,legend=legend,font.main=font.main,font.sub=font.sub),dots[names(dots) %in% args_plot]))
  if (traject) {
    points(nsol[1,x_plane+1],nsol[1,y_plane+1],pch=pch)
    lines(nsol[,x_plane+1],nsol[,y_plane+1],lwd=lwd,col=col)
  }
  if (table) return(nsol)
  f <- state
  f[1:length(f)] <- as.numeric(nsol[nrow(nsol),2:(nvar+1)])
  return(f)
}

newton <- function(state=inistate, parms=pars, odes=feedbackqAOP, time=1, positive=FALSE, jacobian=FALSE, vector=FALSE, plot=FALSE, silent=FALSE, addone=FALSE, ...) {
  # find a steady state
  # if (!is.numeric(x)) x <- index(x,names(state))
  # if (!is.numeric(y)) y <- index(y,names(state))
  One <- ifelse(addone, 1, 0)
  q <- steady(y=state,func=odes,parms=parms,time=time,positive=positive)#, ...
  if (attr(q,"steady")) {
    equ <- q$y
    equ <- ifelse(abs(equ) < 1e-8, 0, equ)
    jac <- jacobian.full(y=equ,func=odes,parms=parms)
    eig <- eigen(jac)
    dom <- max(Re(eig$values))
    if (!silent) {
      print(equ)
      if (dom < 0) cat("Stable point, ")
      else cat("Unstable point, ")
      cat("eigenvalues:\n")
      print(eig$values)
    }
    if (vector) {cat("Eigenvectors:\n"); print(eig$vectors)}
    if (jacobian) {cat("Jacobian:\n"); print(jac)}
    if (plot) {
      if (dom < 0) points(equ[x_plane]+One,equ[y_plane]+One,pch=19)
      else points(equ[x_plane]+One,equ[y_plane]+One,pch=1)
    }
    if (silent) return(list(state=equ,jacobian=jac,values=eig$values,vectors=eig$vectors))
    return(equ)
  }
  cat("No convergence: start closer to a steady state")
  return(NULL)
}

x_continue <- 1; xmin_continue <- 0; xmax_continue <- 1
y_continue <- 1; ymin_continue <- 0; ymax_continue <- 1
log_continue <- ""; addone_continue <- FALSE

continue <- function(state=inistate, parms=pars, odes=feedbackqAOP, step=0.1, x=1, y=2, time=0, xmin=0, xmax=1,ymin=0, ymax=1.1, xlab="", ylab="", log="", col=c("red","black","blue"), lwd= c(2,1,1), addone=FALSE, positive=FALSE, nvar=FALSE, add=FALSE, ...) {  
  # continue a steady state
  dots <- list(...)
  if (!is.null(dots)) {
    unknown <- names(dots[!names(dots) %in% c(args_steady,args_plot)])
    if (length(unknown)>0) warning(paste("Unknown argument(s):",unknown,sep=" "))
    dots_steady <- dots[names(dots) %in% args_steady]
  }else dots_steady <- NULL
  if (add) {
    x <- x_continue
    y <- y_continue
  } else {
    if (!is.numeric(x)) x <- index(x,names(parms))
    if (!is.numeric(y)) y <- index(y,names(state))
    x_continue <<- x
    y_continue <<- y
  }
  
  p0 <- parms[x]
  q0 <- do.call('steady',c(list(y=state,func=odes,parms=parms,time=time,positive=positive),dots_steady))
  if (!attr(q0,"steady"))
    stop("No convergence: start closer to a steady state")
  cat("Starting at",names(parms[x]),"=",parms[x],"with:\n")
  print(q0$y)
  bary <- q0$y[y]
  if (!add) {
    if (missing(xmax) & parms[x] >= 1) xmax <- 2*parms[x]
    if (missing(xmin) & parms[x] < 0) xmin <- 2*parms[x]
    if (!missing(xmin) & xmin >= parms[x]) stop("xmin should be smaller than parameter")
    if (!missing(xmax) & xmax <= parms[x]) stop("xmax should be larger than parameter")
    if (missing(ymax) & bary >= 1.1) ymax <- 2*bary
    if (missing(ymin) & bary < 0) ymin <- 2*bary
    if (!missing(ymin) & ymin >= bary & !addone) stop("ymin should be smaller than y-variable")
    if (!missing(ymax) & ymax <= bary) stop("ymax should be larger than y-variable")
    if (xlab == "") xlab <- names(p0)
    if (ylab == "") {
      ylab <- names(state)[y]
      if (addone) ylab <- paste(ylab,"+ 1")
    }
    do.call('plot',c(list(1,1,type='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=xlab,ylab=ylab,log=log,font.main=font.main,font.sub=font.sub),dots[names(dots) %in% args_plot]))
    xmin_continue <<- xmin; xmax_continue <<- xmax
    ymin_continue <<- ymin; ymax_continue <<- ymax 
    log_continue <<- log; addone_continue <<- addone
  } else {
    xmin <- xmin_continue; xmax <- xmax_continue
    ymin <- ymin_continue; ymax <- ymax_continue
    log <- log_continue; addone <- addone_continue
    if (ymin >= bary & !addone) stop("Initial point below minimum of current y-axis")
    if (ymax <= bary) stop("Initial point above maximum of current y-axis")
  }
  logx <- ifelse(grepl('x',log), TRUE, FALSE)
  COL <- function(s,i) {
    if (!nvar) return(col[i])
    return(col[length(s[s>1e-9])])
  }
  FUN <- function(lastState,lastDom,step) {
    lastP <- p0
    preLastState <- lastState
    nok <- 0
    while (xmin < lastP & lastP < xmax & ymin < lastState[y]+One & lastState[y] < ymax) {
      if (logx) parms[x] <- lastP*(1+step)
      else parms[x] <- lastP + step
      q <- do.call('steady',c(list(y=lastState,func=odes,parms=parms,time=time,positive=positive),dots_steady))
      newState <- q$y  # should be steady state and closeby
      if (attr(q,"steady") & sum(abs(newState-lastState))/(1e-9+sum(abs(lastState))) < 0.1) {
        jac <- jacobian.full(y=newState,func=odes,parms=parms)
        dom <- sign(max(Re(eigen(jac)$values)))
        if (dom != lastDom) cat("Bifurcation at",names(parms[x]),"=",parms[x],"\n")
        if (logx) lines(c(parms[x]/(1+step),parms[x]),c(lastState[y]+One,newState[y]+One), col=COL(lastState,dom+2),lwd=lwd[dom+2])
        else lines(c(parms[x]-step,parms[x]),c(lastState[y]+One,newState[y]+One), col=COL(lastState,dom+2),lwd=lwd[dom+2])
        preLastState <- lastState
        lastState <- newState
        lastDom <- dom
        lastP <- parms[x]
        if (nok > 10 & abs(step) < actualStep) step <- sign(step)*min(2*abs(step),actualStep)
        nok <- nok + 1
      }else{
        nok <- 0
        if (abs(step) > actualStep/100) step <- step/2
        else{ # Go back one step, overpredict, call steady, and turn
          parms[x] <- lastP
          predState <- lastState + 5*(lastState-preLastState)
          q <- do.call('steady',c(list(y=predState,func=odes,parms=parms,time=time,positive=positive),dots_steady))
          newState <- q$y  # should be steady state and not the same
          if (attr(q,"steady") & sum(abs(newState-lastState))/(1e-9+sum(abs(lastState))) > 0.001) {
            cat("Turning point point at",names(parms[x]),"=",parms[x],"\n")
            jac <- jacobian.full(y=newState,func=odes,parms=parms)
            dom <- sign(max(Re(eigen(jac)$values)))
            middle <- (lastState[y]+newState[y])/2
            lines(c(parms[x],parms[x]),c(lastState[y]+One,middle+One), col=COL(lastState,lastDom+2),lwd=lwd[lastDom+2])
            lines(c(parms[x],parms[x]),c(middle+One,newState[y]+One), col=COL(newState,dom+2),lwd=lwd[dom+2])
            step <- -step
            preLastState <- lastState
            lastState <- newState
            lastDom <- dom
            lastP <- parms[x]
          }else{
            cat("Final point at",names(parms[x]),"=",parms[x],"\n")
            cat("If this looks wrong try changing the step size\n")
            break
          }
        }
      }
    }
  }
  One <- ifelse(addone, 1, 0)
  orgWarn <- getOption("warn")
  options(warn = -1)
  jac <- jacobian.full(y=q0$y,func=odes,parms=parms)
  dom <- sign(max(Re(eigen(jac)$values)))
  if (logx) actualStep <- step
  else actualStep <- step*xmax
  FUN(lastState=q0$y,lastDom=dom,actualStep)
  FUN(lastState=q0$y,lastDom=dom,-actualStep)
  options(warn = orgWarn)
  return(NULL)
}

