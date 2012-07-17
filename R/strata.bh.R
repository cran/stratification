`strata.bh` <-
function(x,bh,n=NULL,CV=NULL,Ls=3,certain=NULL,alloc=list(q1=0.5,q2=0,q3=0.5),takenone=0,bias.penalty=1,takeall=0,takeall.adjust=FALSE,rh=rep(1,Ls),model=c("none","loglinear","linear","random"),model.control=list())
{
    # Validation des arguments et initialisation de variables locales
    xgiven <- N <- N1 <- findn <-  L <- q1 <- q2 <- q3 <- beta <- sig2 <-  ph <- pcertain <- gamma <- epsilon <- A <- B <- C <- NULL
    out.check <- checkargs(x=x,bh=bh,n=n,CV=CV,Ls=Ls,certain=certain,alloc=alloc,takenone=takenone,bias.penalty=bias.penalty,takeall=takeall,
          takeall.adjust=takeall.adjust,rh=rh,model=model,model.control=model.control)
    for(i in 1:length(out.check)) assign(names(out.check)[i],out.check[[i]])
    
    bhfull <- c(min(x),bh,max(x)+1)

    # stratification
    Nh <- Nc <- EYh <- EYc <- VYh <- nh <- nh.nonint <- opti.nh <- opti.nhnonint <- TY <- se <- bias <- prop <- stratumID <- NULL
    out.get <- bh2nh(xgiven,N,bhfull,findn,n,CV,L,certain,q1,q2,q3,A,B,C,bias.penalty,takeall.adjust,rh,model,beta,sig2,ph,pcertain,gamma,epsilon)
    for(i in 1:length(out.get)) assign(names(out.get)[i],out.get[[i]])
    
    # Sortie des résultats
    out <- list(Nh=Nh, nh=nh, n=sum(nh)+Nc, nh.nonint=nh.nonint, certain.info=c(Nc=Nc,meanc=EYc),
                opti.nh=opti.nh, opti.nhnonint=opti.nhnonint, meanh=ifelse(Nh==0,NA,EYh),
                varh=ifelse(Nh==0,NA,VYh), mean=TY/N, RMSE=se/N, RRMSE=se/TY, relativebias=bias,
                propbiasMSE=prop, stratumID=stratumID, takeall=length(B), call=match.call(),
                date=date(), args=list(x=xgiven, bh=as.vector(bh), n=n, CV=CV, Ls=Ls,
                    certain=certain, alloc=alloc, takenone=takenone, bias.penalty=bias.penalty,
                    takeall=takeall, rh=rh, model=model[1], model.control=model.control))
    class(out)<-"strata"
    return(out)  
}


`strata.internal` <-
function(x,N,bh,findn,n,CV,Ls,Nc,EYc,alloc,takenone,bias.penalty,takeall,rh,model,model.control)
{ # Équivalent de strata.bh, programmé en C, appelé par strata.LH (Kozak)
    L <- Ls + takenone
    Nh <- nhnonint <- nh <- EYh <- VYh <- rep(0,L)
    opti.nh <- opti.nhnonint <- 0
    if ("none"==model) {
    resuC <- .C("strataCnone",as.double(x),as.integer(N),as.double(bh),as.integer(findn),as.integer(n),as.double(CV),
                as.integer(L),as.integer(Nc),as.double(EYc),as.double(alloc$q1),as.double(alloc$q2),as.double(alloc$q3),
                as.integer(takenone),as.double(bias.penalty),as.integer(takeall),as.double(rh),
                Nh=as.integer(Nh),EYh=as.double(EYh),VYh=as.double(VYh),optinh=as.double(opti.nh),
                optinhnonint=as.double(opti.nhnonint),nhnonint=as.double(nhnonint),nh=as.integer(nh),
                PACKAGE="stratification")
    } else if ("loglinear"==model) {
    resuC <- .C("strataCloglinear",as.double(x),as.integer(N),as.double(bh),as.integer(findn),as.integer(n),as.double(CV),
                as.integer(L),as.integer(Nc),as.double(EYc),as.double(alloc$q1),as.double(alloc$q2),as.double(alloc$q3),
                as.integer(takenone),as.double(bias.penalty),as.integer(takeall),
                as.double(rh),as.double(model.control$beta),as.double(model.control$sig2),as.double(model.control$ph),
                Nh=as.integer(Nh),EYh=as.double(EYh),VYh=as.double(VYh),optinh=as.double(opti.nh),
                optinhnonint=as.double(opti.nhnonint),nhnonint=as.double(nhnonint),nh=as.integer(nh),
                PACKAGE="stratification")
    } else if ("linear"==model) {
    resuC <- .C("strataClinear",as.double(x),as.integer(N),as.double(bh),as.integer(findn),as.integer(n),as.double(CV),
                as.integer(L),as.integer(Nc),as.double(EYc),as.double(alloc$q1),as.double(alloc$q2),as.double(alloc$q3),
                as.integer(takenone),as.double(bias.penalty),as.integer(takeall),
                as.double(rh),as.double(model.control$beta),as.double(model.control$sig2),as.double(model.control$gamma),
                Nh=as.integer(Nh),EYh=as.double(EYh),VYh=as.double(VYh),optinh=as.double(opti.nh),
                optinhnonint=as.double(opti.nhnonint),nhnonint=as.double(nhnonint),nh=as.integer(nh),
                PACKAGE="stratification")
    } else if ("random"==model) {
    resuC <- .C("strataCrandom",as.double(x),as.integer(N),as.double(bh),as.integer(findn),as.integer(n),as.double(CV),
                as.integer(L),as.integer(Nc),as.double(EYc),as.double(alloc$q1),as.double(alloc$q2),as.double(alloc$q3),
                as.integer(takenone),as.double(bias.penalty),as.integer(takeall),as.double(rh),as.double(model.control$epsilon),
                Nh=as.integer(Nh),EYh=as.double(EYh),VYh=as.double(VYh),optinh=as.double(opti.nh),
                optinhnonint=as.double(opti.nhnonint),nhnonint=as.double(nhnonint),nh=as.integer(nh),
                PACKAGE="stratification")
    }
    out <- list(Nh=resuC$Nh,nh=resuC$nh,nhnonint=resuC$nhnonint,opti.nh=resuC$optinh,
                opti.nhnonint=resuC$optinhnonint,meanh=resuC$EYh,varh=resuC$VYh)
    return(out)
}
