## ## for creating a package
.package.Name <- "CARMAgeddon"

setClass( "monocarptr", representation( pointer = "externalptr" ) )

setMethod( "initialize", "monocarptr",
          function(.Object, ...) {
              .Object@pointer <- .Call( "monocarptr__new", ... , PACKAGE=.package.Name)
              .Object
          } )

monocarptr.print <- function(param, dataptr) {
    .Call( "monocarptr__print" , dataptr@pointer , param, PACKAGE=.package.Name)
}
monocarptr.loglik <- function(param, dataptr) {
    .Call( "monocarptr__fun" , dataptr@pointer , param, PACKAGE=.package.Name)
}
monocarptr.grad <- function(param, dataptr) {
    .Call( "monocarptr__grad" , dataptr@pointer , param, PACKAGE=.package.Name)
}
monocarptr.parallelgrad <- function(param, dataptr) {
    .Call( "monocarptr__parallelgrad" , dataptr@pointer , param, PACKAGE=.package.Name)
}
monocarptr.partial <- function(param, dataptr) {
    .Call( "monocarptr__partial" , dataptr@pointer , param, PACKAGE=.package.Name)
}
monocarptr.hess <- function(param, dataptr) {
    .Call( "monocarptr__hess" , dataptr@pointer , param, PACKAGE=.package.Name)
}
monocarptr.hist <- function(param, dataptr) {
    .Call( "monocarptr__hist" , dataptr@pointer , param, PACKAGE=.package.Name)
}
monocarptr.histrand <- function(param, num.iter, dataptr) {
    .Call( "monocarptr__histrand" , dataptr@pointer , param, as.integer(num.iter), PACKAGE=.package.Name)
}
monocarptr.endopt <- function(success) {
    .Call( "monocarptr__endopt" , success, PACKAGE=.package.Name)
}

monocar.cal <- function(init.param, data, exodata, pc.list,
                        verbose = 1,
                        optimizers=c("MMA","BFGS"),
                        maxiters=NULL,
                        tolerances=NULL,
                        epsilon=0,
                        lower.bounds=NULL,
                        upper.bounds=NULL,
                        var.centers=NULL,
                        priors=NULL,
                        newton.raphson=0,
                        time.indices = NULL,
                        compute.hessian=TRUE,
                        compute.partial.scores=TRUE,
                        history.times=NULL,
                        use.parallel=NULL) {

    compute.history <- !is.null(history.times)
    
    mu.const <- pc.list$mu.const
    extra.mu.const <- pc.list$extra.mu.const
    theta.const <- pc.list$theta.const
    alpha.const <- pc.list$alpha.const
    sigma.const <- pc.list$sigma.const
    extra.var.add.const <- pc.list$ev.const
    extra.var.mult.const <- pc.list$evm.const
    extra.var.pow.const <- pc.list$evp.const
    start.const <- pc.list$start.const
    exo.const <- pc.list$exo.const
    sizes <- pc.list$sizes
    
    if (is.null(mu.const))
        mu.const <- matrix(0, sqrt(dim(theta.const)[1]), dim(theta.const)[2])
    if (is.null(extra.mu.const))
        extra.mu.const <- matrix(0, sqrt(dim(theta.const)[1]), dim(theta.const)[2])
    if (is.null(extra.var.add.const))
        extra.var.add.const <- matrix(0, dim(extra.mu.const)[1], dim(extra.mu.const)[2])
    if (is.null(extra.var.mult.const))
        extra.var.mult.const <- cbind(matrix(0, dim(extra.mu.const)[1], dim(extra.mu.const)[2] - 1), 1)
    if (is.null(extra.var.pow.const))
        extra.var.pow.const <- cbind(matrix(0, dim(extra.mu.const)[1], dim(extra.mu.const)[2] - 1), 1)
    
    stopifnot("data.frame" %in% class(data))
    stopifnot(all(c("which.series", "time.start", "time.end", "obs", "var") %in% names(data)))
    
    if (is.null(time.indices))
        time.indices = sort(unique(c(data$time.start, data$time.end)))
    my.order <- order(data$time.end, data$time.start)
    
    stopifnot(length(init.param)+1 == dim(extra.mu.const)[2])
    stopifnot(length(init.param)+1 == dim(theta.const)[2])
    stopifnot(length(init.param)+1 == dim(alpha.const)[2])
    stopifnot(length(init.param)+1 == dim(sigma.const)[2])
    
    if ("house" %in% names(data))
        house.series <- data$house
    else
        house.series <- data$which.series
    
    if (!("which.process" %in% names(data)))
        process.series <- rep(1, data$which.series)
    else if (is.null(data$which.process))
        process.series <- rep(1, data$which.series)
    else
        process.series <- data$which.process
    
    if (!is.null(start.const))
        start.const <- as.matrix(start.const)
    
    if (!is.null(exo.const))
        {
            stopifnot(all(!is.na(exodata)))
            exo.const <- as.matrix(exo.const)
        }
    
    if (is.null(var.centers)) {
        if (sizes[4]==1)
            var.centers <- tapply(data$var,
                                  factor(data$house, levels=0:max(data$house, na.rm=TRUE)),
                                  mean,
                                  na.rm=TRUE)
        else
            var.centers <- tapply(data$var,
                                  factor(data$which.series, levels=0:max(data$which.series, na.rm=TRUE)),
                                  mean,
                                  na.rm=TRUE)

        maxcenters <- max(nrow(extra.var.add.const), nrow(extra.var.mult.const), nrow(extra.var.pow.const))
        if (maxcenters > length(var.centers))
            var.centers <- c(var.centers, rep(1, maxcenters - length(var.centers)))
        var.centers <- ifelse(is.na(var.centers) | (var.centers<= 0.0), 1, var.centers)
    }
    
    screen.width <- getOption("width")
    if (is.null(screen.width))
        screen.width <- -1
    else if (length(screen.width) > 1)
        screen.width <- screen.width[1]

    output.info <- c(compute.partial.scores,
                     compute.hessian,
                     compute.history,
                     verbose)

    {
        num.threads <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", defaultNumThreads())
        if (grepl("^[0-9]+$", num.threads))
            num.threads <- as.numeric(num.threads)
        else
            num.threads <- defaultNumThreads()
        if (is.null(use.parallel))
            use.parallel <- num.threads > 1
        grain.size <- min(16,max(4,2*ceiling(length(c(init.param))/num.threads)))
        mf.raw <- match.call()
        if (is.na(use.parallel)) {
            ptr.data <- list(as.double(init.param),
                             as.integer(output.info),
                             as.numeric(epsilon),
                             start.const,
                             as.matrix(mu.const),
                             as.matrix(extra.mu.const),
                             as.matrix(theta.const),
                             as.matrix(alpha.const),
                             as.matrix(sigma.const),
                             as.matrix(extra.var.add.const),
                             as.matrix(extra.var.mult.const),
                             as.matrix(extra.var.pow.const),
                             exo.const,
                             priors,
                             as.double(time.indices),
                             as.integer(process.series[my.order]),
                             as.integer(data$which.series[my.order]),
                             as.integer(house.series[my.order]),
                             as.double(data$time.start[my.order]),
                             as.double(data$time.end[my.order]),
                             as.matrix(exodata[my.order,]),
                             as.double(data$obs[my.order]),
                             as.double(data$var[my.order]),
                             as.double(var.centers),
                             as.integer(sizes),
                             as.integer(screen.width),
                             as.integer(grain.size))
            names(ptr.data) <- 
                c("param",
                  "output",
                  "epsilon",
                  "start", "mean", "extraMean",
                  "reversion", "simultaneous", "sigma",
                  "extraVarAdd", "extraVarMult", "extraVarPow",
                  "exo",
                  "priors",
                  "timeIndices", "processIndices",
                  "seriesIndices", "houseIndices",
                  "timeStart", "timeEnd",
                  "exoData",
                  "observations", "variances", "varCenters",
                  "sizes", "screenWidth", "grainSize")
            ptr.data$bounds <- list(lower=lower.bounds,
                                    upper=upper.bounds)
            attr(ptr.data, "raw.ptr.data") <- TRUE
            return(ptr.data)
        }
        
        my.point <- new("monocarptr",
                        as.double(init.param),
                        as.integer(output.info),
                        as.numeric(epsilon),
                        start.const, as.matrix(mu.const), as.matrix(extra.mu.const),
                        as.matrix(theta.const), as.matrix(alpha.const), as.matrix(sigma.const),
                        as.matrix(extra.var.add.const),
                        as.matrix(extra.var.mult.const),
                        as.matrix(extra.var.pow.const),
                        exo.const,
                        priors,
                        as.double(time.indices),
                        as.integer(process.series[my.order]),
                        as.integer(data$which.series[my.order]),
                        as.integer(house.series[my.order]),
                        as.double(data$time.start[my.order]),
                        as.double(data$time.end[my.order]),
                        as.matrix(exodata[my.order,]),
                        as.double(data$obs[my.order]),
                        as.double(data$var[my.order]),
                        as.double(var.centers),
                        as.double(history.times),
                        as.integer(sizes),
                        as.integer(screen.width),
                        as.integer(grain.size))

        my.param <- as.double(init.param)
        my.min <- NA

        if (is.null(optimizers))
            optimizers <- c("MMA","BFGS")
        ran.optimizer <- FALSE

        if (length(optimizers) > 0)
            cat("Estimating:\n")
        
        for (i in 1:length(optimizers)) {
            optimizerName <- optimizers[[i]]
            hasOptimizer <- FALSE
            
            if (optimizerName == "DIRECT" ||
                optimizerName == "DIRECTL" ||
                optimizerName == "DIRECT_L")
                optimizerFullName <- "NLOPT_GN_DIRECT_L"
            if (optimizerName == "CRS" || optimizerName == "CRS2")
                {
                    optimizerFullName <- "NLOPT_GN_CRS2_LM"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "MLSL")
                {
                    optimizerFullName <- "NLOPT_G_MLSL_LDS"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "STOGO")
                {
                    optimizerFullName <- "NLOPT_GD_STOGO"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "STOGORAND" || optimizerName == "STOGO_RAND")
                {
                    optimizerFullName <- "NLOPT_GD_STOGO_RAND"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "ISRES")
                {
                    optimizerFullName <- "NLOPT_GN_ISRES"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "ESCH")
                {
                    optimizerFullName <- "NoptimizerNameLOPT_GN_ESCH"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "COBYLA")
                {
                    optimizerFullName <- "NLOPT_LN_COBYLA"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "BOBYQA")
                {
                    optimizerFullName <- "NLOPT_LN_BOBYQA"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "NEWUOA")
                {
                    optimizerFullName <- "NLOPT_LN_NEWUOA"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "PRAXIS")
                {
                    optimizerFullName <- "NLOPT_LN_PRAXIS"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "SUBPLEX" || optimizerName == "SBPLX")
                {
                    optimizerFullName <- "NLOPT_LN_SBPLX"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "NM" || optimizerName == "NELDERMEAD")
                {
                    optimizerFullName <- "NLOPT_LN_NELDERMEAD"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "MMA")
                {
                    optimizerFullName <- "NLOPT_LD_MMA"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "CCSA" || optimizerName == "CCSAQ")
                {
                    optimizerFullName <- "NLOPT_LD_CCSAQ"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "SLSQP")
                {
                    optimizerFullName <- "NLOPT_LD_SLSQP"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "BFGS" || optimizerName == "LBFGS")
                {
                    optimizerFullName <- "NLOPT_LD_LBFGS"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "TNEWTON")
                {
                    optimizerFullName <- "NLOPT_LD_TNEWTON_RESTART"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "VAR2")
                {
                    optimizerFullName <- "NLOPT_LD_VAR2"
                    hasOptimizer <- TRUE
                }
            if (optimizerName == "VAR1")
                {
                    optimizerFullName <- "NLOPT_LD_VAR1"
                    hasOptimizer <- TRUE
                }
            
            if (!hasOptimizer) {
                if (optimizerName != "NONE")
                    cat("\nSkipping unknown optimizer: ", optimizerName, "\n");
            }
            else {
                if (is.null(attr(optimizerName, "opts")))
                    opts <- list("algorithm"=optimizerFullName,
                                 "ftol_abs"=1e-12,
                                 "print_level"=0)
                else
                    opts <- append(list("algorithm"=optimizerFullName),
                                   attr(optimizerName, "opts"))
                if (!("print_level") %in% opts)
                    opts[["print_level"]] <- 0
                if (i <= length(tolerances))
                    opts$xtol_rel = tolerances[i]
                if (i <= length(maxiters))
                    opts$maxeval = maxiters[i]

                if (use.parallel==TRUE)
                    res <- nloptr(x0=my.param,
                                  eval_f=monocarptr.loglik,
                                  eval_grad_f=monocarptr.parallelgrad,
                                  lb=lower.bounds,
                                  ub=upper.bounds,
                                  opts=opts,
                                  dataptr=my.point)
                else
                    res <- nloptr(x0=my.param,
                                  eval_f=monocarptr.loglik,
                                  eval_grad_f=monocarptr.grad,
                                  lb=lower.bounds,
                                  ub=upper.bounds,
                                  opts=opts,
                                  dataptr=my.point)

                ran.optimizer <- TRUE
                
                if (res$status>0) {
                    my.param <- res$solution
                    my.min <- res$objective
                }

                monocarptr.endopt(as.integer(res$status))
                if (verbose > 2)
                    cat(res$message, "\n")
                ## else if (verbose > 0)
                ##     cat(ifelse(res$status>0,"+","-"))
            }
        }

        if (newton.raphson > 0) {
            my.param <- as.double(my.param)
            for (nr.iter in 1:newton.raphson) {
                mygrad <- monocarptr.grad(as.double(my.param), dataptr=my.point)
                myhess <- monocarptr.hess(as.double(my.param), dataptr=my.point)
                myeig <- eigen(myhess)
                my.param <- my.param - as.double(myeig$vectors %*% (t(myeig$vectors) * ifelse(myeig$values <= 0.0, 0.0, 1/myeig$values)) %*% mygrad)
            }
            my.min <- monocarptr.loglik(as.double(my.param), dataptr=my.point)
        }
        
        output <- list(par = my.param,
                       logLik = -my.min)

        if (compute.partial.scores) {
            output$score <- -monocarptr.grad(as.double(my.param), my.point)
            output$grad <- monocarptr.partial(as.double(my.param), my.point)$grad
        }
        
        if (compute.hessian) {
            output$hess <- monocarptr.hess(as.double(my.param), my.point)
        }
        else if (ran.optimizer && (verbose > 0)) {
            cat("\n")
        }

        if (compute.history) {
            my.hist <- monocarptr.hist(as.double(my.param), my.point)
            output$hist.time <- my.hist$hist.time
            output$hist.mean <- my.hist$hist.mean
            output$hist.sd <- my.hist$hist.sd
        }

    }
    
    if (inherits(output$grad, "matrix") && (nrow(output$grad) == length(my.order)))
        output$grad[my.order,] <- output$grad
    return(output)
}

simulate.monocar <- function(object, nsim = 1, seed = NULL,
                                    var=NULL, t1=seq(0,100,1.0), t2=NULL, data=NULL,
                                    series.name = NULL, house.name = NULL, process.name = NULL,
                                    transform = function(x, v) x, transform.var = function(x, v) v,
                                    exo.data=NULL, exovars=NULL, byexo=NULL,
                                    var.centers=NULL, ...) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
        runif(1)
    if (is.null(seed)) 
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    if (is.null(var))
        var <- rep(0, length(t1))
    if (any(is.na(var)))
        stop("Some variances are missing")
    if (any(var<0))
        stop("Some variances are negative")
    if (any(var==Inf))
        stop("Some variances are infinite")
    if (is.null(t2))
        t2 <- t1
    if (nsim==1) {
        new.data <-
            create.ctdata(rep(NA,length(t1)), var, t1, t2, data=data,
                          series.name = series.name, house.name = house.name, process.name = process.name,
                          transform = transform, transform.var = transform.var,
                          time.origin = NULL, time.unit = NULL, exo.data = exo.data)
        sample.data <- simdata.monocar.latent(init=object, data=new.data, exovars=exovars, byexo=byexo, var.centers=var.centers)
        output.data <-
            create.ctdata(sample.data$obs, var, t1, t2, data=data,
                          series.name = series.name, house.name = house.name, process.name = process.name,
                          transform = transform, transform.var = transform.var,
                          time.origin = NULL, time.unit = NULL, exo.data = exo.data)
        attr(output.data, "seed") <- RNGstate
        return(output.data)
    }
    else if (nsim>1) {
        output.list <- list()
        for (i in 1:nsim) {
            output.list[[i]] <-
                simulate.monocar(object, nsim = 1, seed = NULL,
                                 var=var, t1=t1, t2=t2, data=data,
                                 series.name = series.name, house.name = house.name, process.name = process.name,
                                 transform = transform, transform.var = transform.var,
                                 exo.data=exo.data, exovars=exovars, byexo=byexo,
                                 var.centers=var.centers)
            attr(output.list, "seed") <- RNGstate
            return(output.list)
        }
    }
    else if (nsim==0) {
        return(NULL)
    }
    else {
        stop("Bad nsim (number of simulations) argument")
    }
}

simdata.monocar.latent <- function(init,
                                    time.indices=NULL, n.series=NULL, n.processes=NULL,
                                    data=NULL, exovars=NULL, byexo=NULL,
                                    var.centers=NULL, var.by.house=FALSE) {
    if (inherits(init, "monocar"))
        init <- init$estimates

    stopifnot((!is.null(data)) | (!is.null(time.indices)))

    stopifnot((!is.null(data)) | (!is.null(n.series)) | (!is.null(init$theta)) | (!is.null(init$sigma)) | (!is.null(init$mu)))

    restrict <- list()
    
    if (!is.null(init$mu))
        n.series <- length(init$mu)
    else if (!is.null(init$theta))
        n.series <- ncol(init$theta)
    else if (!is.null(init$sigma))
        n.series <- ncol(init$sigma)

    delta.group <- NULL
    mu.by.process <- FALSE
    delta.by.process <- FALSE

    if (is.null(n.processes))
        {
            if (mu.by.process & is.matrix(init$mu))
                n.processes <- nrow(init$mu)
            else if (delta.by.process & is.matrix(init$delta))
                n.processes <- nrow(init$delta)
            else
                n.processes <- 1
        }

    no.data <- is.null(data)
    if (is.null(data))
        data <- data.frame(time.start=time.indices,
                           time.end=time.indices,
                           var=rep(0,length(time.indices)),
                           which.process=rep(n.processes-1,length(time.indices)),
                           which.series=rep(n.series-1,length(time.indices)),
                           house=rep(n.series-1,length(time.indices)))

    if (!("obs" %in% names(data)))
        data$obs <- 0
    
    exodata <- data[,NULL,drop=FALSE]
    if (inherits(exovars, "formula"))
        exodata <- as.data.frame(model.matrix(as.formula(exovars),data))
    else if (inherits(exovars, "character") | inherits(exovars, "integer"))
        exodata <- data[,exovars,drop=FALSE]

    if (is.null(data$house))
        data$house <- data$which.series
    if ("character" %in% class(data$which.series))
        data$which.series <- factor(data$which.series)
    if ("character" %in% class(data$house))
        data$house <- factor(data$house)
    exonames <- names(exodata)
    byexovar <- NULL
    if (!is.null(byexo) & (dim(exodata)[2] > 0))
        if ((byexo != FALSE) & (byexo != "") & (byexo != "none"))
            {
                if ((byexo=="series")|(byexo=="which.series")|(byexo == TRUE))
                    byexovar <- data$which.series
                else if ((byexo=="house")|(byexo=="which.house"))
                    byexovar <- data$house
                else
                    stop('Bad specification of argument "byexo"')
                if ("integer" %in% class(byexovar))
                    byexovar <- factor(byexovar, levels=0:max(byexovar))
                else if (is.character(byexovar))
                    byexovar <- factor(byexovar)
                newexodata <- exodata[,NULL]
                exonames <- NULL
                byexovarnames <- levels(byexovar)
                for (i in 1:length(byexovarnames))
                    {
                        exodatatmp <- as.data.frame(exodata)
                        exonames <- c(exonames, paste(byexovarnames[i], names(exodatatmp), sep=" <- "))
                        names(exodatatmp) <- paste(byexovarnames[i], names(exodatatmp), sep=":")
                        exodatatmp[as.numeric(byexovar)!=i,] <- 0
                        newexodata <- cbind(newexodata, exodatatmp)
                    }
                exodata <- newexodata
                                        #stopifnot(!("byexovar" %in% names(as.data.frame(exodata))))
                                        #exodata <- as.data.frame(model.matrix(~0+byexovar:., as.data.frame(exodata)))
                                        #names(exodata) <- gsub("^byexovar", "", names(exodata))
            }
    n.series <- table(data$which.series)
    n.house <- table(data$house)
    series.house.table <- table(data$which.series, data$house)
    if ("factor" %in% class(data$which.series))
        data$which.series <- as.numeric(data$which.series)-1
    if ("factor" %in% class(data$house))
        data$house <- as.numeric(data$house)-1

    if (is.null(data$which.process))
        data$which.process <- 0
    if ("character" %in% class(data$which.process))
        data$which.process <- factor(data$which.process)  
    if ("factor" %in% class(data$which.process))
        {
            if (any(is.na(data$which.process)))
                data$which.process <- as.numeric(data$which.process)
            else
                data$which.process <- as.numeric(data$which.process)-1
        }
    stopifnot(!(any(is.na(data$which.process)) && any(data$which.process == 0, na.rm=TRUE)))
    if (any(is.na(data$which.process)))
        data$which.process <- ifelse(is.na(data$which.process), 0, as.numeric(data$which.process))

    ##exodata <- exodata[order(data$time.end,data$time.start),,drop=FALSE]

    stopifnot(all(c("time.start","time.end","var","which.process","which.series","house") %in% names(data)))

                                        #stopifnot(min(data$which.series)>=0)
    stopifnot(min(data$house[data$which.series>=0])>=0)
    if (!is.null(init$theta))
        stopifnot(max(data$which.series)<dim(init$theta)[1])    
    if (!is.null(init$sigma)) {
        stopifnot(isSymmetric(init$sigma))
        if (!is.null(init$theta))
            stopifnot(all(dim(init$theta)==dim(init$sigma)))
    }
    if (is.character(restrict$sigma)) {
        if (tolower(restrict$sigma)==substr("diagonal",1,nchar(restrict$sigma)))
            restrict$sigma <- (diag(max(data$which.series)+1)==0)
        else if (tolower(restrict$sigma)==substr("unrestricted",1,nchar(restrict$sigma)))
            restrict$sigma <- matrix(FALSE, max(data$which.series)+1, max(data$which.series)+1)
        else if (tolower(restrict$sigma)==substr("restricted",1,nchar(restrict$sigma)))
            restrict$sigma <- matrix(FALSE, max(data$which.series)+1, max(data$which.series)+1)
        else
            stop('Sigma must be "diagonal", "unrestricted", "restricted", a matrix, or NULL')
    }
    if (is.character(restrict$theta)) {
        if (tolower(restrict$theta)==substr("diagonal",1,nchar(restrict$theta)))
            restrict$theta <- (diag(max(data$which.series)+1)==0)
        else if (tolower(restrict$theta)==substr("unrestricted",1,nchar(restrict$theta)))
            restrict$theta <- matrix(FALSE, max(data$which.series)+1, max(data$which.series)+1)
        else if (tolower(restrict$theta)==substr("restricted",1,nchar(restrict$theta)))
            restrict$theta <- matrix(FALSE, max(data$which.series)+1, max(data$which.series)+1)
        else
            stop('Theta must be "diagonal", "unrestricted", "restricted", a matrix, or NULL')
    }
    if (!is.null(init$delta))
        stopifnot(max(data$house)<length(as.vector(init$delta)))
    if (!is.null(init$var.pow)) {
        if (var.by.house)
            stopifnot(max(data$house)<length(as.vector(init$var.pow)))
        else
            stopifnot(max(data$which.series)<length(as.vector(init$var.pow)))
    }
    if (!is.null(init$var.mult)) {
        if (var.by.house)
            stopifnot(max(data$house)<length(as.vector(init$var.mult)))
        else
            stopifnot(max(data$which.series)<length(as.vector(init$var.mult)))
    }
    if (!is.null(init$var.add)) {
        if (var.by.house)
            stopifnot(max(data$house)<length(as.vector(init$var.add)))
        else
            stopifnot(max(data$which.series)<length(as.vector(init$var.add)))
    }

    endodata <- data[,c("time.start","time.end","obs","var","which.process","which.series","house")]

    paramconst <- getparamconst(endodata, exodata, init, restrict, delta.group, mu.by.process, delta.by.process, var.by.house)

    pc.list <- paramconst$pc.list

    mu.const <- pc.list$mu.const
    extra.mu.const <- pc.list$extra.mu.const
    theta.const <- pc.list$theta.const
    alpha.const <- pc.list$alpha.const
    sigma.const <- pc.list$sigma.const
    extra.var.add.const <- pc.list$ev.const
    extra.var.mult.const <- pc.list$evm.const
    extra.var.pow.const <- pc.list$evp.const
    start.const <- pc.list$start.const
    exo.const <- pc.list$exo.const
    sizes <- pc.list$sizes
    
    if (is.null(mu.const))
        mu.const <- matrix(0, sqrt(dim(theta.const)[1]), dim(theta.const)[2])
    if (is.null(extra.mu.const))
        extra.mu.const <- matrix(0, sqrt(dim(theta.const)[1]), dim(theta.const)[2])
    if (is.null(extra.var.add.const))
        extra.var.add.const <- matrix(0, dim(extra.mu.const)[1], dim(extra.mu.const)[2])
    if (is.null(extra.var.mult.const))
        extra.var.mult.const <- cbind(matrix(0, dim(extra.mu.const)[1], dim(extra.mu.const)[2] - 1), 1)
    if (is.null(extra.var.pow.const))
        extra.var.pow.const <- cbind(matrix(0, dim(extra.mu.const)[1], dim(extra.mu.const)[2] - 1), 1)
    
    stopifnot("data.frame" %in% class(data))
    stopifnot(all(c("which.series", "time.start", "time.end", "var") %in% names(data)))
    
    if (is.null(time.indices))
        time.indices <- sort(unique(c(data$time.start, data$time.end)))
    
    init.param <- paramconst$param
    
    stopifnot(length(init.param)+1 == dim(extra.mu.const)[2])
    stopifnot(length(init.param)+1 == dim(theta.const)[2])
    stopifnot(length(init.param)+1 == dim(alpha.const)[2])
    stopifnot(length(init.param)+1 == dim(sigma.const)[2])
    
    if ("house" %in% names(data))
        house.series <- data$house
    else
        house.series <- data$which.series
    
    if ("which.process" %in% names(data))
        process.series <- data$which.process
    else
        process.series <- data$which.series
    
    if (!is.null(start.const))
        start.const <- as.matrix(start.const)
    
    if (!is.null(exo.const))
        {
            stopifnot(all(!is.na(exodata)))
            exo.const <- as.matrix(exo.const)
        }
    
    if (is.null(var.centers)) {
        var.centers <- tapply(data$var,
                              factor(data$house, levels=0:max(data$house, na.rm=TRUE)),
                              mean,
                              na.rm=TRUE)
        var.centers <- ifelse(is.na(var.centers) | (var.centers<= 0.0), 1, var.centers)
    }

    data$house <- house.series
    data$which.process <- process.series

    if (no.data) {
        data <- data[NULL,,drop=FALSE]
        exodata <- exodata[NULL,,drop=FALSE]
    }

    time.indices.complete <- sort(unique(c(time.indices, data$time.start, data$time.end)))

    screen.width <- getOption("width")
    if (is.null(screen.width))
        screen.width <- -1
    else if (length(screen.width) > 1)
        screen.width <- screen.width[1]
    
    output <- .Call("ousim", as.double(init.param),
                    start.const, as.matrix(mu.const), as.matrix(extra.mu.const),
                    as.matrix(theta.const), as.matrix(alpha.const), as.matrix(sigma.const),
                    as.matrix(extra.var.add.const),
                    as.matrix(extra.var.mult.const),
                    as.matrix(extra.var.pow.const),
                    exo.const,
                    as.double(time.indices.complete),
                    as.integer(data$which.process),
                    as.integer(data$which.series),
                    as.integer(data$house),
                    as.double(data$time.start),
                    as.double(data$time.end),
                    as.matrix(exodata),
                    as.double(data$var),
                    as.double(var.centers),
                    as.integer(sizes),
                    as.integer(screen.width),
                    PACKAGE=.package.Name)

    output$time.indices <- time.indices.complete
    ##output$time.indices.complete <- time.indices.complete
    output$interval.start <- time.indices.complete[-length(time.indices.complete)]
    output$interval.end <- time.indices.complete[-1]
    output$exonames <- exonames
    
    return(output)

}

create.ctdata <- function(x, v, t1, t2 = NULL,
                          data = NULL,
                          series.name = NULL,
                          house.name = NULL,
                          process.name = NULL,
                          transform = NULL,
                          transform.var = NULL,
                          time.origin = NULL,
                          time.unit = NULL,
                          exo.data = NULL,
                          inclusive.end.date = FALSE) {
    if (is.data.frame(data) | is.environment(data)) {
        x <- eval(substitute(x), data)
        v <- eval(substitute(v), data)
        t1 <- eval(substitute(t1), data)
        t2 <- eval(substitute(t2), data)
        series.name <- eval(substitute(series.name), data)
        house.name <- eval(substitute(house.name), data)
        process.name <- eval(substitute(process.name), data)
        exo.data.orig <- eval(substitute(exo.data), data)
        if (inherits(exo.data.orig, "formula"))
            exo.data.orig <- model.frame(exo.data, data)
        time.origin <- eval(substitute(time.origin), data)
        time.unit <- eval(substitute(time.unit), data)
    }
    else {
        exo.data.orig <- exo.data
        if (inherits(exo.data, "formula"))
            exo.data.orig <- model.frame(exo.data, data)
    }
    if (is.null(transform)) {
        transform <- function(x,v) return(x)
        if (is.null(transform.var))
            transform.var <- function(x,v) return(v)
    }
    else if (inherits(transform, "character")) {
        if (!is.null(transform.var))
            cat("Warning: specifying a transform.var when using a named transformation is usually ill-advised\n")
        if (all(tolower(transform) %in% c("arcsin","arc-sin","asin","sin-1","binomial"))) {
            transform <- function(x,v) asin((1-1/5000)*(2*x/v-1))
            if (is.null(transform.var))
                transform.var <- function(x,v) return(1/v)
        }
        else if (all(tolower(transform) %in% c("sqrt","squareroot","square-root","sin-1","poisson"))) {
            transform <- function(x,v) return(sqrt(x+3/8))
            if (is.null(transform.var))
                transform.var <- function(x,v) return(rep(0.25, length(x)))
        }
        else {
            stop("unknown transformation name")
        }
    }
    if (is.null(transform.var)) {
        transform.var <- function(x,v) return(v)
    }
    if (is.matrix(exo.data.orig))
        exo.data.orig <- as.data.frame(exo.data.orig)
    if (is.null(t2))
        t2 <- t1
    if (inherits(t2, "Date") && (inclusive.end.date))
        t2 <- t2+1
    if (is.null(series.name))
        series.name <- deparse(substitute(x))
    if (is.null(house.name))
        house.name <- series.name
    if(!is.null(exo.data.orig)) {
        if (is.data.frame(exo.data.orig)) {
            if (is.null(names(exo.data.orig)))
                names(exo.data.orig) <- paste0("X", 1:nrow(exo.data.orig))
            exo.data <- exo.data.orig
        }
        else if (is.factor(exo.data.orig) | is.numeric(exo.data.orig)) {
            exo.data.frame <- as.data.frame(exo.data.orig)
            names(exo.data.frame) <- deparse(substitute(exo.data))
            exo.data <- exo.data.frame
        }
        else if (is.character(exo.data.orig)) {
            if (length(x) == length(exo.data.orig))
                {
                    exo.data.frame <- as.data.frame(exo.data.orig)
                    names(exo.data.frame) <- deparse(substitute(exo.data))
                }
            else if (all(exo.data.orig %in% names(data)))
                {
                    exo.data.frame <- data[,exo.data.orig]
                }
            exo.data <- exo.data.frame
        }
        else {
            stop("Bad exo.data type")
        }
        stopifnot(dim(as.data.frame(exo.data))[1] == length(x))
    }
    if (is.null(house.name))
        house.name <- series.name
    if (is.null(process.name))
        process.name <- NA
    if (!is.null(time.origin)) {
        if (is.null(time.unit) | !("difftime" %in% class(time.unit))) {
            t1 <- as.numeric(t1-time.origin)
            t2 <- as.numeric(t2-time.origin)
        }
        else {
            t1 <- difftime(t1,time.origin,units=units(time.unit))
            t2 <- difftime(t2,time.origin,units=units(time.unit))
        }
    }
    if (!is.null(time.unit)) {
        if (is.null(time.origin) & (!(is.null(time.unit) | !("difftime" %in% class(time.unit))))) {
            if (inherits(t1, "POSIXct"))
                time.origin <- as.POSIXct("1970-01-01 00:00.00 UTC")
            else if (inherits(t1, "POSIXlt"))
                time.origin <- as.POSIXlt("1970-01-01 00:00.00 UTC")
            else if (inherits(t1, "Date"))
                time.origin <- as.Date("1970-01-01")
            else
                time.origin <- as.POSIXct("1970-01-01 00:00.00 UTC")
            t1 <- as.numeric(difftime(t1,time.origin,units=units(time.unit)))
            t2 <- as.numeric(difftime(t2,time.origin,units=units(time.unit)))
        }
        t1 <- t1/as.numeric(time.unit)
        t2 <- t2/as.numeric(time.unit)
    }
    if(is.null(exo.data))
        return.object <-
            data.frame(time.start=as.numeric(t1), time.end=as.numeric(t2),
                       obs=transform(x,v), var=transform.var(x,v),
                       which.series=series.name, house=house.name, which.process = process.name)
    else
        return.object <-
            data.frame(time.start=as.numeric(t1), time.end=as.numeric(t2),
                       obs=transform(x,v), var=transform.var(x,v),
                       which.series=series.name, house=house.name, which.process = process.name,
                       as.data.frame(exo.data))
    if (!is.null(time.origin))
        attr(return.object, "time.origin") <- time.origin
    if (!is.null(time.unit))
        attr(return.object, "time.unit") <- time.unit
    class(return.object) <- c("ct.data.frame", "data.frame")
    return(return.object)
}


join.ctdata <- function(...) {
    x <- list(...)
    if (length(x)==0)
        return(NULL)
    if ((length(x)==1) && (inherits(x[[1]], "list")))
        x <- x[[1]]
    x <- x[!sapply(x, is.null)]
    stopifnot(all(sapply(x, inherits, what="ct.data.frame")))
    if (length(unique(unlist(sapply(x, attr, which="time.origin")))) > 1)
        stop("Different values of time.origin")
    if (length(unique(unlist(sapply(x, attr, which="time.unit")))) > 1)
        stop("Different values of time.unit")
    if (length(x)==1)
        return(x[[1]])
    req.names <- c("time.start", "time.end", "obs", "var",
                   "which.series", "house", "which.process")
    all.series.numeric <- TRUE
    all.series.factor <- FALSE
    all.series.levels <- NULL
    all.series <- NULL
    all.house.numeric <- TRUE
    all.house.factor <- FALSE
    all.house.levels <- NULL
    all.house <- NULL
    all.process.numeric <- TRUE
    all.process.factor <- FALSE
    all.process.levels <- NULL
    all.process <- NULL
    all.exo <- NULL
    for (y in x) {
        if (all(is.na(y$which.series))) {
            all.series <- c(all.series, rep(NA, length(y$which.series)))
            all.series.levels <- c(all.series.levels, NA)
        }
        else if (is.factor(y$which.series)) {
            all.series <- c(all.series, as.character(y$which.series))
            all.series.levels <- c(all.series.levels, levels(y$which.series))
            all.series.numeric <- FALSE
            all.series.factor <- TRUE
        }
        else if (is.character(y$which.series)) {
            all.series <- c(all.series, as.character(y$which.series))
            all.series.levels <- c(all.series.levels, unique(y$which.series))
            all.series.numeric <- FALSE
        }
        else if (is.numeric(y$which.series)) {
            all.series <- c(all.series, as.character(y$which.series))
            all.series.levels <- c(all.series.levels, unique(y$which.series))
        }
        if (all(is.na(y$house))) {
            all.house <- c(all.house, rep(NA, length(y$house)))
            all.house.levels <- c(all.house.levels, NA)
        }
        else if (is.factor(y$house)) {
            all.house <- c(all.house, as.character(y$house))
            all.house.levels <- c(all.house.levels, levels(y$house))
            all.house.numeric <- FALSE
            all.house.factor <- TRUE
        }
        else if (is.character(y$house)) {
            all.house <- c(all.house, as.character(y$house))
            all.house.levels <- c(all.house.levels, unique(y$house))
            all.house.numeric <- FALSE
        }
        else if (is.numeric(y$house)) {
            all.house <- c(all.house, as.character(y$house))
            all.house.levels <- c(all.house.levels, unique(y$house))
        }
        if (all(is.na(y$which.process))) {
            all.process <- c(all.process, rep(NA, length(y$which.process)))
            all.process.levels <- c(all.process.levels, NA)
        }
        else if (is.factor(y$which.process)) {
            all.process <- c(all.process, as.character(y$which.process))
            all.process.levels <- c(all.process.levels, levels(y$which.process))
            all.process.numeric <- FALSE
            all.process.factor <- TRUE
        }
        else if (is.character(y$which.process)) {
            all.process <- c(all.process, as.character(y$which.process))
            all.process.levels <- c(all.process.levels, unique(y$which.process))
            all.process.numeric <- FALSE
        }
        else if (is.numeric(y$which.process)) {
            all.process <- c(all.process, as.character(y$which.process))
            all.process.levels <- c(all.process.levels, unique(y$which.process))
        }
        all.exo <- c(all.exo, setdiff(names(y), req.names))
    }
    all.series.levels <- unique(all.series.levels)
    all.house.levels <- unique(all.house.levels)
    all.process.levels <- unique(all.process.levels)
    all.exo <- unique(all.exo)
    output <- NULL
    for (y in x) {
        if (!all.series.numeric)
            y$which.series <- as.character(y$which.series)
        if (!all.house.numeric)
            y$house <- as.character(y$house)
        if (!all.process.numeric)
            y$which.process <- as.character(y$which.process)
        for (z in setdiff(all.exo, names(y))) {
            y[[z]] <- NA
        }
        output <- rbind(output, y[,c(req.names, all.exo)])
    }
    if (all.series.factor)
        output$which.series <-
            factor(all.series, levels=setdiff(all.series.levels,NA))
    if (all.house.factor)
        output$house <-
            factor(all.house, levels=setdiff(all.house.levels,NA))
    if (all.process.factor)
        output$which.process <-
            factor(all.process, levels=setdiff(all.process.levels,NA))
    class(output) <- c("ct.data.frame", "data.frame")
    return(output)
}

getparamconst <- function(data, exodata, estimates, restrict,
                          delta.group=TRUE, mu.by.process=FALSE,
                          delta.by.process=FALSE, var.by.house=FALSE) {
    param <- rests <- varnum <- NULL

    sizes <- c(NA,NA,NA,NA)
    if (var.by.house)
        sizes[4] <- 1
    else
        sizes[4] <- 0
    if (is.null(estimates$delta)) {
        sizes[3] <- max(data$house)+1
    }
    else if (is.matrix(estimates$delta)) {
        sizes[3] <- nrow(estimates$delta)
        if (delta.by.process)
            sizes[1] <- ncol(estimates$delta)
    }
    else {
        sizes[3] <- length(estimates$delta)
    }
    
    if (!is.null(estimates$theta))
        sizes[2] <- nrow(estimates$theta)
    else if (!is.null(estimates$sigma))
        sizes[2] <- nrow(estimates$sigma)
    else if (is.null(estimates$mu))
        sizes[2] <- max(data$which.series)+1
    else if (is.matrix(estimates$mu)) {
        sizes[2] <- nrow(estimates$mu)
        if (mu.by.process)
            sizes[1] <- ncol(estimates$mu)
    }
    else {
        sizes[2] <- length(estimates$mu)
    }

    if (is.na(sizes[1]))
        sizes[1] <- max(data$which.process)+1
    if (is.na(sizes[1]))
        sizes[1] <- 1
    
    if (is.logical(delta.group)) {
        if (delta.group & all(!is.na(data$which.series)) & all(!is.na(data$house)) & ((length(unique(data$which.process))<=1) | !(mu.by.process | delta.by.process))) {
            if (all(apply(table(data$which.series, data$house)>0,2,sum)==1)) {
                group.table <- table(factor(data$which.series,levels=0:max(data$which.series)),
                                     factor(data$house,levels=0:max(data$house)))
                delta.group <- lapply(apply(group.table > 0, 1, function(x) list(which(x)-1)), function(x) x[[1]])
            }
            else {
                delta.group <- NULL
            }
        }
        else {
            delta.group <- NULL
        }
    }
    
    muguess <- rep(0, max(data$which.series)+1)
    if (is.null(estimates$mu) & is.null(restrict$mu) &
        !(is.null(delta.group))) {
        muguess <- rep(0, max(data$which.series)+1)
        param <- c(param, muguess)
        rests <- c(rests, rep(T, length(muguess)))
    }
    else if (is.null(estimates$mu) & is.null(restrict$mu) &
             !(is.null(estimates$delta) & is.null(restrict$delta))) {
        muguess <- rep(0, max(data$which.series)+1)
        param <- c(param, muguess)
        rests <- c(rests, rep(T, length(muguess)))
    }
    else if (is.null(estimates$mu)) {
        if (mu.by.process)
            muguess <- tapply(data$obs[data$which.series>=0],
                              list(factor(data$which.series[data$which.series>=0],
                                          levels=0:max(data$which.series)),
                                   factor(data$which.process[data$which.series>=0],
                                          levels=0:max(data$which.process))),
                              FUN=mean, na.rm=T)
        else
            muguess <- tapply(data$obs[data$which.series>=0],
                              factor(data$which.series[data$which.series>=0],
                                     levels=0:max(data$which.series)),
                              FUN=mean, na.rm=T)
        param <- c(param, ifelse(is.na(muguess), 0, muguess))
        if (is.null(restrict$mu))
            rests <- c(rests, is.na(muguess))
        else
            rests <- c(rests, as.vector(restrict$mu))
    }
    else {
        muguess <- estimates$mu
        param <- c(param, as.vector(estimates$mu))
        if (is.null(restrict$mu))
            rests <- c(rests, rep(FALSE, length(as.vector(estimates$mu))))
        else
            rests <- c(rests, as.vector(restrict$mu))
    }
    varnum <- c(varnum, length(rests))

    if (is.null(estimates$delta) & is.null(restrict$delta) & (is.null(delta.group) & !all(rests))) {
        param <- c(param, rep(0,max(data$house)+1))
        rests <- c(rests, rep(TRUE,max(data$house)+1))
    }
    else if (is.null(estimates$delta)) {
        if (delta.by.process)
            deltaguess <- tapply(data$obs[data$which.series>=0 & data$house>=0] -
                                     muguess[1+data$which.series[data$which.series>=0 & data$house>=0]],
                                 list(factor(data$house[data$which.series>=0 & data$house>=0],
                                             levels=0:max(data$house)),
                                      factor(data$which.process[data$which.series>=0 & data$house>=0],
                                             levels=0:max(data$which.process))),
                                 FUN=mean, na.rm=T)
        else
            deltaguess <- tapply(data$obs[data$which.series>=0 & data$house>=0] -
                                     muguess[1+data$which.series[data$which.series>=0 & data$house>=0]],
                                 factor(data$house[data$which.series>=0 & data$house>=0],
                                        levels=0:max(data$house)),
                                 FUN=mean, na.rm=T)
        param <- c(param, ifelse(is.na(deltaguess), 0, deltaguess))
        if (is.null(restrict$delta))
            rests <- c(rests, is.na(deltaguess))
        else
            rests <- c(rests, as.vector(restrict$delta))
    }
    else {
        param <- c(param, as.vector(estimates$delta))
        if (is.null(restrict$delta))
            rests <- c(rests, rep(FALSE, length(as.vector(estimates$delta))))
        else
            rests <- c(rests, as.vector(restrict$delta))
    }
    varnum <- c(varnum, length(rests))

    if (var.by.house)
        var.type.length <- sizes[3]
    else
        var.type.length <- sizes[2]
    
    esttheta <- estimates$theta
    if (is.null(estimates$theta)) {
        if (is.null(estimates$mu))
            esttheta <- diag(max(data$which.series)+1)
        else if (is.matrix(estimates$mu))
            esttheta <- diag(ncol(estimates$mu))
        else
            esttheta <- diag(length(estimates$mu))
    }
    param <- c(param, as.vector(esttheta))
    if (is.null(restrict$theta))
        rests <- c(rests, rep(FALSE, length(as.vector(esttheta))))
    else
        rests <- c(rests, as.vector(restrict$theta))
    varnum <- c(varnum, length(rests))

    alpha.guess <- diag(dim(esttheta)[1])
    if (!is.null(estimates$alpha))
        alpha.guess <- estimates$alpha
    stopifnot(all(dim(esttheta) == dim(alpha.guess)))
    param <- c(param, as.vector(alpha.guess))
    if (is.null(restrict$alpha))
        rests <- c(rests, rep(TRUE, length(as.vector(alpha.guess))))
    else
        rests <- c(rests, as.vector(restrict$alpha))
    varnum <- c(varnum, length(rests))
    
    estsigorig <- estsig <- estimates$sigma
    restsig <- restrict$sigma
    if (is.null(estimates$sigma)) {
        estsigorig <- estsig <- diag(dim(esttheta)[1])
        if (is.null(restsig))
            restsig <- (diag(dim(esttheta)[1]) == 0)
    }
    else if (any(is.na(estimates$sigma)) & !is.null(restrict$sigma))
        {
            esut <- estsig[upper.tri(estsig, TRUE)]
            rut <- restrict$sigma[upper.tri(estsig, TRUE)]
            esut <- ifelse(is.na(esut),
                           esut[!is.na(esut)][cumsum(!is.na(esut) & !rut) + pmin(0, rut)],
                           esut)
            estsig[upper.tri(estsig, TRUE)] <- esut
        }
    eigvals <- eigen(estsig, symmetric=TRUE)$values
    if (any(eigvals < 0)) {
        eigvals <- zapsmall(eigvals)
        eigvecs <- eigen(estsig, symmetric=TRUE)$vectors
        if (any(eigvals < 0))
            stop("Sigma is not positive definite")
        estsig <- tcrossprod(eigvecs %*% diag(eigvals), eigvecs)
    }
    cholpsd <- suppressWarnings(chol(estsig, pivot=TRUE))
    sigma.vec <- cholpsd[order(attr(cholpsd, "pivot")),order(attr(cholpsd, "pivot"))][upper.tri(estsig, TRUE)]
    sigma.vec <- ifelse(is.na(estsigorig[upper.tri(estsig, TRUE)]), NA, sigma.vec)
    param <- c(param, as.vector(sigma.vec))
    if (is.null(restsig) &
        all(sigma.vec[upper.tri(estsig, FALSE)[upper.tri(estsig, TRUE)]]==0))
        rests <- c(rests, (upper.tri(estsig, FALSE)[upper.tri(estsig, TRUE)]))
    else if (is.null(restsig))
        rests <- c(rests, rep(FALSE, length(sigma.vec)))
    else
        rests <- c(rests, (restsig)[upper.tri(restsig, TRUE)])
                                        #rests <- c(rests, (restrict$sigma|t(restrict$sigma))[upper.tri(restrict$sigma, TRUE)])
    varnum <- c(varnum, length(rests))


    if (is.null(estimates$var.add)) {
        if (is.null(restrict$var.add)) {
            param <- c(param, rep(0, var.type.length))
            rests <- c(rests, rep(TRUE, var.type.length))
        }
        else {
            param <- c(param, rep(0, length(as.vector(restrict$var.add))))
            rests <- c(rests, as.vector(restrict$var.add))
        }
    }
    else {
        param <- c(param, sqrt(abs(as.vector(estimates$var.add))))
        if (is.null(restrict$var.add))
            rests <- c(rests, rep(FALSE, length(as.vector(estimates$var.add))))
        else
            rests <- c(rests, as.vector(restrict$var.add))
    }
    varnum <- c(varnum, length(rests))

    if (is.null(estimates$var.mult)) {
        if (is.null(restrict$var.mult)) {
            param <- c(param, rep(1, var.type.length))
            rests <- c(rests, rep(TRUE, var.type.length))
        }
        else {
            param <- c(param, rep(1, length(as.vector(restrict$var.mult))))
            rests <- c(rests, as.vector(restrict$var.mult))
        }
    }
    else {
        param <- c(param, sqrt(abs(as.vector(estimates$var.mult))))
        if (is.null(restrict$var.mult))
            rests <- c(rests, rep(FALSE, length(as.vector(estimates$var.mult))))
        else
            rests <- c(rests, as.vector(restrict$var.mult))
    }
    varnum <- c(varnum, length(rests))

    if (is.null(estimates$var.pow)) {
        if (is.null(restrict$var.pow)) {
            param <- c(param, rep(1, var.type.length))
            rests <- c(rests, rep(TRUE, var.type.length))
        }
        else {
            param <- c(param, rep(1, length(as.vector(restrict$var.pow))))
            rests <- c(rests, as.vector(restrict$var.pow))
        }
    }
    else {
        param <- c(param, as.vector(estimates$var.pow))
        if (is.null(restrict$var.pow))
            rests <- c(rests, rep(FALSE, length(as.vector(estimates$var.pow))))
        else
            rests <- c(rests, as.vector(restrict$var.pow))
    }
    varnum <- c(varnum, length(rests))

    if (is.null(exodata)) {
        stopifnot(is.null(estimates$gamma))
        stopifnot(is.null(restrict$gamma))
    }
    else if (is.null(estimates$gamma)) {
        param <- c(param, rep(0, dim(exodata)[2]))
        rests <- c(rests, rep(FALSE, dim(exodata)[2]))
    }
    else {
        stopifnot(length(estimates$gamma) == dim(exodata)[2])
        param <- c(param, estimates$gamma)
        if (is.null(restrict$gamma)) {
            rests <- c(rests, rep(FALSE, length(as.vector(estimates$gamma))))
        }
        else {
            stopifnot(length(restrict$gamma) == dim(exodata)[2])
            rests <- c(rests, as.vector(restrict$gamma))
        }
    }
    varnum <- c(varnum, length(rests))
    
    if (!is.null(estimates$start)) {
        param <- c(param, as.vector(estimates$start))
        if (is.null(restrict$start))
            rests <- c(rests, rep(FALSE, length(as.vector(estimates$start))))
        else
            rests <- c(rests, as.vector(restrict$start))
    }
    varnum <- c(varnum, length(rests))

    stopifnot(length(param)==length(rests))
    
    rest.mat <- diag(sum(!rests & !is.na(param))+1)[ifelse(rests & !is.na(param), sum(!rests & !is.na(param))+1,
                                                           ifelse(is.na(param),
                                                                  cumsum(!rests & !is.na(param))+pmin(0,rests),
                                                                  cumsum(!rests & !is.na(param)))),]
    rest.mat[,sum(!rests & !is.na(param))+1] <- ifelse((rests==TRUE) & !is.na(param), param, 0)

                                        #stopifnot(all((rests != FALSE) | !is.na(param)))
    
    mu.const <- rest.mat[1:varnum[1],,drop=F]
    extra.mu.const <- rest.mat[(varnum[1]+1):varnum[2],,drop=F]
    theta.const <- rest.mat[(varnum[2]+1):varnum[3],,drop=F]
    alpha.const <- rest.mat[(varnum[3]+1):varnum[4],,drop=F]
    sigma.const <- rest.mat[(varnum[4]+1):varnum[5],,drop=F]
    ev.const <- rest.mat[(varnum[5]+1):varnum[6],,drop=F]
    evm.const <- rest.mat[(varnum[6]+1):varnum[7],,drop=F]
    evp.const <- rest.mat[(varnum[7]+1):varnum[8],,drop=F]
    if (varnum[8] < varnum[9])
        exo.const <- rest.mat[(varnum[8]+1):varnum[9],,drop=F]
    else
        exo.const <- NULL
    start.const <- NULL
    if (!is.null(estimates$start))
        start.const <- rest.mat[(varnum[9]+1):varnum[10],,drop=F]

    pc.list <- list(mu.const=mu.const,
                    extra.mu.const=extra.mu.const,
                    theta.const=theta.const,
                    alpha.const=alpha.const,
                    sigma.const=sigma.const,
                    ev.const=ev.const,
                    evm.const=evm.const,
                    evp.const=evp.const,
                    exo.const=exo.const,
                    start.const=start.const,
                    sizes=sizes)

    return(list(param=param[!rests & !is.na(param)], pc.list=pc.list, delta.group=delta.group))
}


monocar.hist <- function(data, model, history.times=NULL, estimates=NULL, exovars=NULL, byexo="series", delta.group=TRUE, verbose=1, mu.by.process=FALSE, delta.by.process=FALSE, var.by.house=FALSE, ...) {
    if (is.null(history.times))
        history.times <- sort(unique(c(data$time.start, data$time.end)))
    else if (is.integer(history.times) && length(history.times) == 1)
        history.times <- seq(min(data$time.start), max(data$time.end), length.out=history.times)
    stopifnot(!(is.null(model) & is.null(estimates)))
    if (is.null(estimates)) {
        estimates <- model$estimates
    }
    if (is.null(exovars) & !is.null(model)) {
        exovars <- model$exovars
    }
    if (is.null(byexo) & !is.null(model)) {
        byexo <- model$byexo
    }
    if (is.null(delta.group) & !is.null(model)) {
        delta.group <- model$delta.group
    }
    exodata <- data[,NULL,drop=FALSE]
    if (inherits(exovars, "formula"))
        exodata <- as.data.frame(model.matrix(as.formula(exovars),data))
    else if (inherits(exovars, "character") | inherits(exovars, "integer"))
        exodata <- data[,exovars,drop=FALSE]

    if ("monocar" %in% class(estimates)) {
        restrict <- estimates$restrict
        estimates <- estimates$estimates
    }
    else {
        restrict <- list()
    }
    
    if (is.null(data$house))
        data$house <- data$which.series
    if ("character" %in% class(data$which.series))
        data$which.series <- factor(data$which.series)
    if ("character" %in% class(data$house))
        data$house <- factor(data$house)
    exonames <- names(exodata)
    byexovar <- NULL
    if (!is.null(byexo) & (dim(exodata)[2] > 0))
        if ((byexo != FALSE) & (byexo != "") & (byexo != "none"))
            {
                if ((byexo=="series")|(byexo=="which.series")|(byexo == TRUE))
                    byexovar <- data$which.series
                else if ((byexo=="house")|(byexo=="which.house"))
                    byexovar <- data$house
                else
                    stop('Bad specification of argument "byexo"')
                if ("integer" %in% class(byexovar))
                    byexovar <- factor(byexovar, levels=0:max(byexovar))
                else if (is.character(byexovar))
                    byexovar <- factor(byexovar)
                newexodata <- exodata[,NULL]
                exonames <- NULL
                byexovarnames <- levels(byexovar)
                for (i in 1:length(byexovarnames))
                    {
                        exodatatmp <- as.data.frame(exodata)
                        exonames <- c(exonames, paste(byexovarnames[i], names(exodatatmp), sep=" <- "))
                        names(exodatatmp) <- paste(byexovarnames[i], names(exodatatmp), sep=":")
                        exodatatmp[as.numeric(byexovar)!=i,] <- 0
                        newexodata <- cbind(newexodata, exodatatmp)
                    }
                exodata <- newexodata
                                        #stopifnot(!("byexovar" %in% names(as.data.frame(exodata))))
                                        #exodata <- as.data.frame(model.matrix(~0+byexovar:., as.data.frame(exodata)))
                                        #names(exodata) <- gsub("^byexovar", "", names(exodata))
            }
    if ("factor" %in% class(data$which.series))
        data$which.series <- as.numeric(data$which.series)-1
    if ("factor" %in% class(data$house))
        data$house <- as.numeric(data$house)-1
    
    if (is.null(data$which.process))
        data$which.process <- 0
    if ("character" %in% class(data$which.process))
        data$which.process <- factor(data$which.process)  
    if ("factor" %in% class(data$which.process))
        {
            if (any(is.na(data$which.process)))
                data$which.process <- as.numeric(data$which.process)
            else
                data$which.process <- as.numeric(data$which.process)-1
        }
    stopifnot(!(any(is.na(data$which.process)) && any(data$which.process == 0, na.rm=TRUE)))
    if (any(is.na(data$which.process)))
        data$which.process <- ifelse(is.na(data$which.process), 0, as.numeric(data$which.process))
    
    exodata <- exodata[order(data$time.end,data$time.start),,drop=FALSE]

    stopifnot(all(c("time.start","time.end","obs","var","which.process","which.series","house") %in% names(data)))

                                        #stopifnot(min(data$which.series)>=0)    
    stopifnot(min(data$house[data$which.series>=0])>=0)
    stopifnot(isSymmetric(estimates$sigma))
    stopifnot(all(dim(estimates$theta)==dim(estimates$sigma)))
    stopifnot(max(data$which.series)<dim(estimates$theta)[1])
    if (!is.null(estimates$mu))
        stopifnot(max(data$which.series)<length(as.vector(estimates$mu)))
    if (!is.null(estimates$delta))
        stopifnot(max(data$house)<length(as.vector(estimates$delta)))
    if (var.by.house) {
        if (!is.null(estimates$var.pow))
            stopifnot(max(data$house)<length(as.vector(estimates$var.pow)))
        if (!is.null(estimates$var.mult))
            stopifnot(max(data$house)<length(as.vector(estimates$var.mult)))
        if (!is.null(estimates$var.add))
            stopifnot(max(data$house)<length(as.vector(estimates$var.add)))
    }
    else {
        if (!is.null(estimates$var.pow))
            stopifnot(max(data$which.series)<length(as.vector(estimates$var.pow)))
        if (!is.null(estimates$var.mult))
            stopifnot(max(data$which.series)<length(as.vector(estimates$var.mult)))
        if (!is.null(estimates$var.add))
            stopifnot(max(data$which.series)<length(as.vector(estimates$var.add)))
    }
    
    endodata <- data[order(data$time.end,data$time.start),
                     c("time.start","time.end","obs","var","which.process","which.series","house")]

    paramconst <- getparamconst(endodata, exodata, estimates, restrict, delta.group, mu.by.process, delta.by.process, var.by.house)

    model.est <- monocar.cal(paramconst$param,
                             endodata, exodata,
                             paramconst$pc.list,
                             optimizers='NONE',
                             verbose=verbose,
                             compute.hessian=FALSE,
                             compute.partial.scores=FALSE,
                             history.times=history.times,
                             ...)
    history.return <- list(param = estimates,
                           time = model.est$hist.time,
                           mean = model.est$hist.mean,
                           sd = model.est$hist.sd)
    class(history.return) <- "monocar.hist"
    return(history.return)
}

monocar.estimate <- function(data, init=list(), restrict=list(),
                             subset=NULL, exovars=NULL, byexo="series",
                             enforce.bounds=TRUE, lower.bounds=NULL, upper.bounds=NULL,
                             delta.group=TRUE, verbose=1,
                             mu.by.process=FALSE, delta.by.process=FALSE, var.by.house=FALSE,
                             remove.unused=TRUE, ...) {
    subset.obs <- eval(substitute(subset), data)
    if (!is.null(subset.obs))
        data <- data[subset.obs,]
    exodata <- data[,NULL,drop=FALSE]
    if (inherits(exovars, "formula"))
        exodata <- as.data.frame(model.matrix(as.formula(exovars),data))
    else if (inherits(exovars, "character") | inherits(exovars, "integer"))
        exodata <- data[,exovars,drop=FALSE]

    if (is.null(data$house))
        data$house <- data$which.series
    if (("character" %in% class(data$which.series)) |
        (is.factor(data$which.series) & remove.unused))
        data$which.series <- factor(data$which.series)
    if (("character" %in% class(data$house)) |
        (is.factor(data$house) & remove.unused))
        data$house <- factor(data$house)
    byexovar <- NULL
    exonames <- names(exodata)
    exonamesorig <- names(exodata)
    byexosimp <- "none"
    if (!is.null(byexo) & (dim(exodata)[2] > 0))
        if ((byexo != FALSE) & (byexo != "") & (byexo != "none"))
            {
                if ((byexo=="series")|(byexo=="which.series")|(byexo == TRUE))
                    {
                        byexovar <- data$which.series
                        byexosimp <- "series"
                    }
                else if ((byexo=="house")|(byexo=="which.house"))
                    {
                        byexovar <- data$house
                        byexosimp <- "house"
                    }
                else
                    {
                        stop('Bad specification of argument "byexo"')
                    }
                if ("integer" %in% class(byexovar))
                    byexovar <- factor(byexovar, levels=0:max(byexovar))
                else if (is.character(byexovar))
                    byexovar <- factor(byexovar)
                newexodata <- exodata[,NULL]
                exonames <- NULL
                byexovarnames <- levels(byexovar)
                for (i in 1:length(byexovarnames))
                    {
                        exodatatmp <- as.data.frame(exodata)
                        exonames <- c(exonames, paste(byexovarnames[i], names(exodatatmp), sep=" <- "))
                        names(exodatatmp) <- paste(byexovarnames[i], names(exodatatmp), sep=":")
                        exodatatmp[as.numeric(byexovar)!=i,] <- 0
                        newexodata <- cbind(newexodata, exodatatmp)
                    }
                exodata <- newexodata
                                        #stopifnot(!("byexovar" %in% names(as.data.frame(exodata))))
                                        #exodata <- as.data.frame(model.matrix(~0+byexovar:., as.data.frame(exodata)))
                                        #names(exodata) <- gsub("^byexovar", "", names(exodata))
            }
    n.series <- table(data$which.series)
    n.house <- table(data$house)
    series.house.table <- table(data$which.series, data$house)
    if ("factor" %in% class(data$which.series))
        data$which.series <- as.numeric(data$which.series)-1
    if ("factor" %in% class(data$house))
        data$house <- as.numeric(data$house)-1
    if (any(!is.na(data$which.process)))
        process.series.table <- table(data$which.process, data$which.series)
    else
        process.series.table <- NULL

    if (is.null(data$which.process))
        data$which.process <- 0
    if ("character" %in% class(data$which.process))
        data$which.process <- factor(data$which.process)  
    if ("factor" %in% class(data$which.process))
        {
            if (any(is.na(data$which.process)))
                data$which.process <- as.numeric(data$which.process)
            else
                data$which.process <- as.numeric(data$which.process)-1
        }
    stopifnot(!(any(is.na(data$which.process)) && any(data$which.process == 0, na.rm=TRUE)))
    if (any(is.na(data$which.process)))
        data$which.process <- ifelse(is.na(data$which.process), 0, as.numeric(data$which.process))

    exodata <- exodata[order(data$time.end,data$time.start),,drop=FALSE]

    stopifnot(all(c("time.start","time.end","obs","var","which.process","which.series","house") %in% names(data)))

                                        #stopifnot(min(data$which.series)>=0)
    stopifnot(min(data$house[data$which.series>=0])>=0)
    if (!is.null(init$theta))
        stopifnot(max(data$which.series)<dim(init$theta)[1])
    if (!is.null(init$sigma)) {
        stopifnot(isSymmetric(init$sigma))
        if (!is.null(init$theta))
            stopifnot(all(dim(init$theta)==dim(init$sigma)))
    }
    if (is.character(restrict$sigma)) {
        if (tolower(restrict$sigma)==substr("diagonal",1,nchar(restrict$sigma)))
            restrict$sigma <- (diag(max(data$which.series)+1)==0)
        else if (tolower(restrict$sigma)==substr("unrestricted",1,nchar(restrict$sigma)))
            restrict$sigma <- matrix(FALSE, max(data$which.series)+1, max(data$which.series)+1)
        else if (tolower(restrict$sigma)==substr("restricted",1,nchar(restrict$sigma)))
            restrict$sigma <- matrix(FALSE, max(data$which.series)+1, max(data$which.series)+1)
        else
            stop('Sigma must be "diagonal", "unrestricted", "restricted", a matrix, or NULL (=diagonal)')
    }
    if (is.character(restrict$theta)) {
        if (tolower(restrict$theta)==substr("diagonal",1,nchar(restrict$theta)))
            restrict$theta <- (diag(max(data$which.series)+1)==0)
        else if (tolower(restrict$theta)==substr("unrestricted",1,nchar(restrict$theta)))
            restrict$theta <- matrix(FALSE, max(data$which.series)+1, max(data$which.series)+1)
        else if (tolower(restrict$theta)==substr("restricted",1,nchar(restrict$theta)))
            restrict$theta <- matrix(FALSE, max(data$which.series)+1, max(data$which.series)+1)
        else
            stop('Theta must be "diagonal", "unrestricted", "restricted", a matrix, or NULL (=unrestricted)')
    }
    if (!is.null(init$delta))
        stopifnot(max(data$house)<length(as.vector(init$delta)))
    if (var.by.house) {
        if (!is.null(init$var.pow))
            stopifnot(max(data$house)<length(as.vector(init$var.pow)))
        if (!is.null(init$var.mult))
            stopifnot(max(data$house)<length(as.vector(init$var.mult)))
        if (!is.null(init$var.add))
            stopifnot(max(data$house)<length(as.vector(init$var.add)))
    }
    else  {
        if (!is.null(init$var.pow))
            stopifnot(max(data$which.series)<length(as.vector(init$var.pow)))
        if (!is.null(init$var.mult))
            stopifnot(max(data$which.series)<length(as.vector(init$var.mult)))
        if (!is.null(init$var.add))
            stopifnot(max(data$which.series)<length(as.vector(init$var.add)))
    }

    endodata <- data[order(data$time.end,data$time.start),
                     c("time.start","time.end","obs","var","which.process","which.series","house")]

    paramconst <- getparamconst(endodata, exodata, init, restrict, delta.group, mu.by.process, delta.by.process, var.by.house)
    
    if (!is.null(attr(init,"force.start"))) {
        if (length(attr(init,"force.start")) != length(paramconst$param))
            stop(paste("Error: starting values have",length(attr(init,"force.start")),"elements, but need",length(paramconst$param),"elements"))
        paramconst$param <- as.double(attr(init,"force.start"))
    }
    
    mu.const <- paramconst$pc.list$mu.const
    extra.mu.const <- paramconst$pc.list$extra.mu.const
    theta.const <- paramconst$pc.list$theta.const
    alpha.const <- paramconst$pc.list$alpha.const
    sigma.const <- paramconst$pc.list$sigma.const
    ev.const <- paramconst$pc.list$ev.const
    evm.const <- paramconst$pc.list$evm.const
    evp.const <- paramconst$pc.list$evp.const
    start.const <- paramconst$pc.list$start.const
    exo.const <- paramconst$pc.list$exo.const
    delta.group <- paramconst$delta.group
    
    if (is.null(lower.bounds))
        lower.bounds <- rep(-Inf, dim(ev.const)[2]-1)
    if (enforce.bounds>1) {
        if (dim(evp.const)[1] > 0)
            lower.bounds <-
                pmax(lower.bounds,
                     ifelse(apply(rbind(evp.const)[,-dim(evp.const)[2],drop=FALSE]>0,2, sum)>0,
                            0, -Inf))
    }
    if (enforce.bounds>0) {
        if (dim(rbind(ev.const,evm.const))[1] > 0)
            lower.bounds <-
                pmax(lower.bounds,
                     ifelse(apply(rbind(ev.const,
                                        evm.const)[,-dim(ev.const)[2],drop=FALSE]>0,2, sum)>0,
                            0, -Inf))
    }
    if (all(lower.bounds == -Inf))
        lower.bounds <- NULL

    if (!is.null(delta.group)) {
        stopifnot(inherits(delta.group, "list"))
        stopifnot(length(delta.group)<=dim(mu.const)[1])
        stopifnot(max(sapply(delta.group, function(x) if (is.null(x)) 0 else max(x)))<dim(extra.mu.const)[1])
    }
    
    model.est <- monocar.cal(paramconst$param,
                             endodata, exodata,
                             paramconst$pc.list,
                             lower.bounds = lower.bounds,
                             upper.bounds = upper.bounds,
                             verbose = verbose,
                             history.times = NULL,
                             ...)

    if (!is.null(attr(model.est, "raw.ptr.data"))) {
        attr(model.est, "raw.ptr.data") <- NULL
        return(model.est)
    }
    
    if (!is.null(delta.group)) {
        for (i in 1:length(delta.group)) {
            whichgroups <- unique(delta.group[[i]])+1
            if (length(whichgroups)>0) {
                meanconst <- apply(extra.mu.const[whichgroups,,drop=FALSE],2,mean)
                mu.const[i,] <- mu.const[i,]+meanconst
                for (j in whichgroups)
                    extra.mu.const[j,] <- extra.mu.const[j,]-meanconst
            }
        }
    }

    initorig <- init
    if (is.null(init$theta))
        init$theta <- diag(sqrt(nrow(theta.const)))
    if (is.null(init$sigma))
        init$sigma <- diag(nrow(init$theta))
    log.likelihood <- model.est$logLik
    par.orig <- as.vector(model.est$par)
    par <- c(par.orig, 1)
    sigma.chol.est <- matrix(0, dim(init$sigma)[1], dim(init$sigma)[2])
    sigma.chol.est[upper.tri(sigma.chol.est, TRUE)] <- as.vector(sigma.const %*% par)
    estimates <- list(mu = as.vector(mu.const %*% par),
                      delta = as.vector(extra.mu.const %*% par),
                      theta = matrix(as.vector(theta.const %*% par), dim(init$theta)),
                      alpha = matrix(as.vector(alpha.const %*% par), dim(init$theta)),
                      sigma = crossprod(sigma.chol.est),
                      var.add = as.vector(ev.const %*% par)^2,
                      var.mult = as.vector(evm.const %*% par)^2,
                      var.pow = as.vector(evp.const %*% par))
    if (!is.null(start.const)) {
        estimates$start <- as.vector(start.const %*% par)
    }
    if (!is.null(exo.const)) {
        estimates$gamma <- as.vector(exo.const %*% par)
        names(estimates$gamma) <- names(exodata)
        if ((!is.null(byexo)) & (length(estimates$gamma) > 0))
            {
                if ((byexosimp != "none") & (length(levels(byexovar))>0) &
                    ((length(estimates$gamma) %% length(levels(byexovar)))==0))
                    {
                        estimates$gamma <- matrix(estimates$gamma, nrow=length(levels(byexovar)), byrow=TRUE)
                        colnames(estimates$gamma) <- exonamesorig
                        rownames(estimates$gamma) <- levels(byexovar)
                    }
            }
    }
    if (is.null(model.est$hess)) {
        vcov.orig <- NA
        hess.orig <- NA
        standard.errors <- NA
        p.values <- NA
    }
    else {
        tol.solve = max(.Machine$double.eps^4, sqrt(sqrt(.Machine$double.xmin)))
        if (any(!is.finite(model.est$hess))) {
            warning("Hessian contains missing values; can't compute standard errors")
            vcov.orig <- NA
            standard.errors <- NA
            p.values <- NA
        }
        if (any(abs(diag(model.est$hess)) < tol.solve)) {
            pseudo.solve <- function(hess) {
                nm <- abs(diag(hess)) >= tol.solve
                vcov.orig <- diag(Inf, nrow(hess), ncol(hess))
                vcov.orig[nm,nm] <- solve(hess[nm,nm], tol=tol.solve)
                return(vcov.orig)
            }
            vcov.orig <- try(pseudo.solve(model.est$hess), silent=TRUE)
            if (!inherits(vcov.orig, "try-error"))
                warning("Hessian has zeros on diagonal\nSome standard errors will be infinite.")
        } else {
            vcov.orig <- try(solve(model.est$hess, tol=tol.solve), silent=TRUE)
        }
        hess.orig <- model.est$hess
        if (inherits(vcov.orig, "try-error")) {
            warning(paste0(as.character(attr(vcov.orig, "condition")), "No standard errors computed"))
            vcov.orig <- NA
            standard.errors <- NA
            p.values <- NA
        }
        else {
            vcov <- cbind(rbind(vcov.orig,0),0)
            `%vcov%` <- function(x, y) {
                return(x %*% vcov %*% y)
            }
            if (any(is.infinite(diag(vcov)))) {
                vcov <- ifelse(is.infinite(vcov), sign(vcov)*(0+1i), vcov)
                `%vcov%` <- function(x, y) {
                    return(ifelse(Im(x %*% vcov %*% y) != 0.0, Inf*sign(Im(x %*% vcov %*% y)), Re(x %*% vcov %*% y)))
                }
            }
            sigma.se <- estimates$sigma*NA
            for (i in 1:dim(estimates$sigma)[1]) {
                for (j in i:dim(estimates$sigma)[1]) {
                    a.const <- sigma.const[(i*(i-1)/2+1):(i*(i+1)/2),,drop=F]
                    b.const <- sigma.const[(j*(j-1)/2+1):(j*(j-1)/2+i),,drop=F]
                    a.param <- a.const %*% par
                    b.param <- b.const %*% par
                    if (all((apply(a.const[,-dim(a.const)[2],drop=F]==0, 1, all) & apply(b.const[,-dim(b.const)[2],drop=F]==0, 1, all)) |
                                apply(a.const==0, 1, all) |
                                    apply(b.const==0, 1, all))) {
                        sigma.se[j,i] <- sigma.se[i,j] <- NA
                    }
                    else {
                        sigma.se[j,i] <- sigma.se[i,j] <- sqrt(as.vector((t(rbind(b.param,a.param)) %*% rbind(a.const,b.const)) %vcov% (t(rbind(a.const,b.const)) %*% rbind(b.param,a.param))))
                    }
                }
            }
            standard.errors <- list(mu = ifelse(apply(mu.const[,-dim(mu.const)[2],drop=F]==0,1,all),
                                        NA,
                                        sqrt(diag(mu.const %vcov% t(mu.const)))),
                                    delta = ifelse(apply(extra.mu.const[,-dim(extra.mu.const)[2],drop=F]==0,1,all),
                                        NA,
                                        sqrt(diag(extra.mu.const %vcov% t(extra.mu.const)))),
                                    theta = matrix(ifelse(apply(theta.const[,-dim(theta.const)[2],drop=F]==0,1,all),
                                        NA,
                                        sqrt(diag(theta.const %vcov% t(theta.const)))), dim(init$theta)),
                                    alpha = matrix(ifelse(apply(alpha.const[,-dim(alpha.const)[2],drop=F]==0,1,all),
                                        NA,
                                        sqrt(diag(alpha.const %vcov% t(alpha.const)))), dim(init$theta)),
                                    sigma = sigma.se,
                                    var.add = ifelse(apply(ev.const[,-dim(ev.const)[2],drop=F]==0,1,all),
                                                     NA,
                                                     2 * abs(as.vector(ev.const %*% par)) *
                                                     sqrt(diag(ev.const %vcov% t(ev.const)))),
                                    var.mult = ifelse(apply(evm.const[,-dim(evm.const)[2],drop=F]==0,1,all),
                                        NA,
                                        2 * abs(as.vector(evm.const %*% par)) *
                                        sqrt(diag(evm.const %vcov% t(evm.const)))),
                                    var.pow = ifelse(apply(evp.const[,-dim(evp.const)[2],drop=F]==0,1,all),
                                        NA,
                                        sqrt(diag(evp.const %vcov% t(evp.const)))))
            p.values <- list(mu=2*pnorm(abs(estimates$mu/standard.errors$mu), lower.tail=F),
                             delta=2*pnorm(abs(estimates$delta/standard.errors$delta), lower.tail=F),
                             theta=2*pnorm(abs(estimates$theta/standard.errors$theta), lower.tail=F),
                             alpha=2*pnorm(abs(estimates$alpha/standard.errors$alpha), lower.tail=F),
                             sigma=2*pnorm(abs(estimates$sigma/standard.errors$sigma), lower.tail=F),
                             var.add=2*pnorm(abs(estimates$var.add/standard.errors$var.add), lower.tail=F),
                             var.mult=2*pnorm(abs(estimates$var.mult/standard.errors$var.mult), lower.tail=F),
                             var.pow=2*pnorm(abs(estimates$var.pow/standard.errors$var.pow), lower.tail=F))
            if (!is.null(start.const)) {
                standard.errors$start <- ifelse(apply(start.const[,-dim(start.const)[2],drop=F]==0,1,all),
                                                NA,
                                                sqrt(diag(start.const %vcov% t(start.const))))
                p.values$start <- 2*pnorm(abs(estimates$start/standard.errors$start), lower.tail=F)
            }
            if (!is.null(exo.const)) {
                standard.errors$gamma <- ifelse(apply(exo.const[,-dim(exo.const)[2],drop=F]==0,1,all),
                                                NA,
                                                sqrt(diag(exo.const %vcov% t(exo.const))))
                names(standard.errors$gamma) <- names(exodata)
                if ((!is.null(byexo)) & (length(standard.errors$gamma) > 0))
                    {
                        if ((byexo != FALSE) & (length(levels(byexovar))>0) &
                            ((length(standard.errors$gamma) %% length(levels(byexovar)))==0))
                            {
                                standard.errors$gamma <- matrix(standard.errors$gamma, nrow=length(levels(byexovar)), byrow=TRUE)
                                colnames(standard.errors$gamma) <- exonamesorig
                            }
                    }
                p.values$gamma <- 2*pnorm(abs(estimates$gamma/standard.errors$gamma), lower.tail=F)
            }
            if (nrow(estimates$theta)==nrow(series.house.table)) {
                colnames(estimates$theta) <- rownames(estimates$theta) <- rownames(series.house.table)
                colnames(standard.errors$theta) <- rownames(standard.errors$theta) <- rownames(series.house.table)
                colnames(p.values$theta) <- rownames(p.values$theta) <- rownames(series.house.table)
            }
            if (nrow(estimates$sigma)==nrow(series.house.table)) {
                colnames(estimates$sigma) <- rownames(estimates$sigma) <- rownames(series.house.table)
                colnames(standard.errors$sigma) <- rownames(standard.errors$sigma) <- rownames(series.house.table)
                colnames(p.values$sigma) <- rownames(p.values$sigma) <- rownames(series.house.table)
            }
            if (!is.null(estimates$alpha))
                if (nrow(estimates$alpha)==nrow(series.house.table)) {
                    colnames(estimates$alpha) <- rownames(estimates$alpha) <- rownames(series.house.table)
                    colnames(standard.errors$alpha) <- rownames(standard.errors$alpha) <- rownames(series.house.table)
                    colnames(p.values$alpha) <- rownames(p.values$alpha) <- rownames(series.house.table)
                }
            if (is.matrix(estimates$mu)) {
                if (nrow(estimates$mu)==nrow(series.house.table)) {
                    colnames(estimates$mu) <- rownames(estimates$mu) <- rownames(series.house.table)
                    colnames(standard.errors$mu) <- rownames(standard.errors$mu) <- rownames(series.house.table)
                    colnames(p.values$mu) <- rownames(p.values$mu) <- rownames(series.house.table)
                }
            }
            else {
                if (length(estimates$mu)==nrow(series.house.table)) {
                    names(estimates$mu) <- rownames(series.house.table)
                    names(standard.errors$mu) <- rownames(series.house.table)
                    names(p.values$mu) <- rownames(series.house.table)
                }
            }
            if (is.matrix(estimates$delta)) {
                if (nrow(estimates$delta)==nrow(series.house.table)) {
                    colnames(estimates$delta) <- rownames(estimates$delta) <- rownames(series.house.table)
                    colnames(standard.errors$delta) <- rownames(standard.errors$delta) <- rownames(series.house.table)
                    colnames(p.values$delta) <- rownames(p.values$delta) <- rownames(series.house.table)
                }
            }
            else {
                if (length(estimates$delta)==nrow(series.house.table)) {
                    names(estimates$delta) <- rownames(series.house.table)
                    names(standard.errors$delta) <- rownames(series.house.table)
                    names(p.values$delta) <- rownames(series.house.table)
                }
            }
            if (length(estimates$var.add)==ncol(series.house.table)) {
                names(estimates$var.add) <- colnames(series.house.table)
                names(standard.errors$var.add) <- colnames(series.house.table)
                names(p.values$var.add) <- colnames(series.house.table)
            }
            if (length(estimates$var.mult)==ncol(series.house.table)) {
                names(estimates$var.mult) <- colnames(series.house.table)
                names(standard.errors$var.mult) <- colnames(series.house.table)
                names(p.values$var.mult) <- colnames(series.house.table)
            }
            if (length(estimates$var.pow)==ncol(series.house.table)) {
                names(estimates$var.pow) <- colnames(series.house.table)
                names(standard.errors$var.pow) <- colnames(series.house.table)
                names(p.values$var.pow) <- colnames(series.house.table)
            }
        }
    }
    
    returnobject <- list(coef=par.orig, vcov=vcov.orig,
                         grad=model.est$grad, score=model.est$score, hess=hess.orig,
                         constraints=paramconst$pc.list,
                         exovars=exovars, byexo=byexosimp, exonames=exonames,
                         delta.group=delta.group,
                         logLik=log.likelihood,
                         n.series=n.series,
                         n.house=n.house,
                         series.house.table=series.house.table,
                         process.series.table=process.series.table,
                         estimates=estimates,
                         standard.errors=standard.errors,
                         p.values=p.values,
                         restrict=restrict)
    class(returnobject) <- c("monocar")
    return(returnobject)
}

logLik.monocar <- function(object, ...) {
    if (!missing(...)) 
        warning("extra arguments discarded")
    p <- length(object$coef)
    val <- object$logLik
    attr(val, "df") <- p
    class(val) <- "logLik"
    val
}

summary.monocar <- function(object, reverse.offdiag = TRUE, ...) {
    exoseries <- ((object$byexo == "series") & (is.matrix(object$estimates$gamma)))
    estimates <- object$estimates
    if (reverse.offdiag) {
        estimates$theta <- estimates$theta *
            (diag(2, nrow(object$estimates$theta))-1)
    }
    else {
        estimates$theta <- estimates$theta *
            matrix(1, nrow(object$estimates$theta),
                   ncol(object$estimates$theta))
    }
    if (exoseries)
        exoseries <- (nrow(as.matrix(object$estimates$gamma)) == nrow(object$estimates$theta))
    mulab <- "(Mean)"
    if (length(estimates$mu) > nrow(estimates$theta)) {
        nproc <- (length(estimates$mu) %/% nrow(estimates$theta))
        if (is.null(object$process.series.table) ||
            length(rownames(object$process.series.table)) != nproc)
            mulab <- paste0("(Mean, process ",1:nproc,")")
        else
            mulab <- paste0("(Mean, ",rownames(object$process.series.table),")")
    }
    if (exoseries)
    {
            labfunc <- function(x) c(t(cbind(x$theta,x$gamma,matrix(x$mu, nrow(x$theta)))))
            cmat <- 
                cbind("Estimate"=labfunc(estimates),
                      "Std.Err"=labfunc(object$standard.errors),
                      "Z value"=labfunc(object$estimates)/labfunc(object$standard.errors),
                      "Pr(>z)"=labfunc(object$p.values))
            rownames(cmat) <- 
                c(paste(rep(rownames(object$series.house.table),
                            each=nrow(object$series.house.table)+ncol(object$estimates$gamma)+length(mulab)),
                        rep(c(paste("<-",c(rownames(object$series.house.table),
                                           colnames(object$estimates$gamma))), mulab),
                            times=nrow(object$series.house.table))))
        }
    else
        {
            labfunc <- function(x) c(t(cbind(x$theta,matrix(x$mu, nrow(x$theta)),as.vector(x$gamma))))
            cmat <- 
                cbind("Estimate"=labfunc(estimates),
                      "Std.Err"=labfunc(object$standard.errors),
                      "Z value"=labfunc(object$estimates)/labfunc(object$standard.errors),
                      "Pr(>z)"=labfunc(object$p.values))
            rownames(cmat) <- 
                c(paste(rep(rownames(object$series.house.table), each=nrow(object$series.house.table)+length(mulab)),
                        rep(c(paste("<-",rownames(object$series.house.table)), mulab),
                            times=nrow(object$series.house.table))),
                  object$exonames)
        }
    output <- list(estimates=cmat,
                   logLik=object$logLik,
                   N=apply(object$series.house.table,1,sum),
                   reverse.offdiag=reverse.offdiag)
    class(output) <- "summary.monocar"
    return(output)
}


print.monocar <- function(x, reverse.offdiag = TRUE, digits = max(3L, getOption("digits") - 3L), ...) {
    estimates <- x$estimates
    if (reverse.offdiag) {
        cat("", "Estimates (off-diagonal elements of Theta reversed in sign):",sep="\n")
        estimates$theta <- estimates$theta *
            (diag(2, nrow(x$estimates$theta))-1)
    }
    else {
        cat("", "Estimates:",sep="\n")
        estimates$theta <- estimates$theta *
            matrix(1, nrow(x$estimates$theta),
                   ncol(x$estimates$theta))
    }
    exoseries <- ((x$byexo == "series") & (is.matrix(x$estimates$gamma)))
    if (exoseries)
        exoseries <- (nrow(as.matrix(x$estimates$gamma)) == nrow(x$estimates$theta))
    if (exoseries)
        {
            labfunc <- function(z)
                c(t(cbind(z$theta,z$gamma,z$mu)))
            output <- labfunc(estimates)
            names(output) <- 
                c(paste(rep(rownames(x$series.house.table),
                            each=nrow(x$series.house.table)+ncol(x$estimates$gamma)+1),
                        rep(c(paste("<-",c(rownames(x$series.house.table),
                                           colnames(x$estimates$gamma))), "(Mean)"),
                            times=nrow(x$series.house.table))))
        }
    else
        {
            labfunc <- function(z)
                c(t(cbind(z$theta,z$mu,as.vector(z$gamma))))
            output <- labfunc(estimates)
            names(output) <- 
                c(paste(rep(rownames(x$series.house.table), each=nrow(x$series.house.table)+1),
                        rep(c(paste("<-",rownames(x$series.house.table)), "(Mean)"),
                            times=nrow(x$series.house.table))),
                  x$exonames)
        }
    print(output, digits=digits, ...)
    cat("\n")
}

print.summary.monocar <- function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"), ...) {
    cat("", "Estimates:",sep="\n")
    printCoefmat(x$estimates,
                 digits = digits, signif.stars = signif.stars, 
                 na.print = "NA", ...)
    if (x$reverse.offdiag) {
        if (getOption("width")>=70)
            cat("Off-diagonal elements of Theta reversed in sign to aid interpretation.", sep="\n")
        else
            cat("Off-diagonal elements of Theta reversed in sign.", sep="\n")
    }
    cat("",
        paste("log-likelihood:", formatC(x$logLik, 1, format="f")),
        paste("N:", paste(names(x$N), x$N, collapse=", ")),
        "",
        sep="\n")
}

monocar.ptr <- function(...) {
    raw.ptr.data <- monocar.estimate(..., use.parallel=NA, verbose=0)
    Data <- list(param.init = raw.ptr.data$param)
    Data$mon.names <- Data$parm.names <- paste0("P[",1:length(raw.ptr.data$param),"]")
    Data$y <- raw.ptr.data$observations
    Data$n <- length(Data$y)
    Data$bounds <- raw.ptr.data$bounds
    Data$ptr <- with(raw.ptr.data, new("monocarptr",
                                 param, output, epsilon,
                                 start, mean, extraMean,
                                 reversion, simultaneous, sigma,
                                 extraVarAdd, extraVarMult, extraVarPow,
                                 exo,
                                 priors,
                                 timeIndices, processIndices,
                                 seriesIndices, houseIndices,
                                 timeStart, timeEnd,
                                 exoData,
                                 observations, variances, varCenters,
                                 double(0L),
                                 sizes, screenWidth, grainSize))
    return(Data)
}
