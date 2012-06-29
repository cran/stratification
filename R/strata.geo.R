`strata.geo` <-
function(x,n=NULL,CV=NULL,Ls=3,certain=NULL,alloc=list(q1=0.5,q2=0,q3=0.5),rh=rep(1,Ls),model=c("none","loglinear","linear","random"),model.control=list())
{
    # Validation des arguments et initialisation de variables locales
    xgiven <- N <- findn <- L <- q1 <- q2 <- q3 <- beta <- sig2 <- ph <- pcertain <- gamma <- epsilon <- NULL
    takenone <- bias.penalty <- takeall.adjust <- A <- B <- C <- NULL
    out.check <- checkargs(x=x,n=n,CV=CV,Ls=Ls,certain=certain,alloc=alloc,rh=rh,model=model,model.control=model.control)
        for(i in 1:length(out.check)) assign(names(out.check)[i],out.check[[i]])
    if (any(x==0)) stop("the geometric method accepts only positive 'x' values") 
    
    # Détermination des bornes
    a<-x[1]
    r<-(x[length(x)]/x[1])^(1/L)
    bhfull<-a*r^(0:L)
    bhfull[L+1]<-bhfull[L+1]+1
    
    # stratification
    Nh <- Nc <- EYh <- EYc <- VYh <- nh <- nh.nonint <- opti.nh <- opti.nhnonint  <- TY <- se <- bias <- prop <- stratumID <- NULL
    out.get <- bh2nh(xgiven,N,bhfull,findn,n,CV,L,certain,q1,q2,q3,A,B,C,bias.penalty,takeall.adjust,rh,model,beta,sig2,ph,pcertain,gamma,epsilon)
    for(i in 1:length(out.get)) assign(names(out.get)[i],out.get[[i]])
    # Note : se est la racine carrée de la variance de TY
    # on donne en sortie la moyenne = TY/N
    # donc RMSE est le root mean squared error de la moyenne, il vaut donc se/N
    # CV ou RRMSE est le relative root mean squared error de la moyenne ou de la somme : se/TY = (se/N)/(TY/N)
    
    # Sortie des résultats
    out<-list(Nh=Nh,nh=nh,n=sum(nh)+Nc,nh.nonint=nh.nonint,certain.info=c(Nc=Nc,meanc=EYc),opti.nh=opti.nh,opti.nhnonint=opti.nhnonint,bh=as.vector(bhfull[-c(1,L+1)]),
              meanh=ifelse(Nh==0,NA,EYh),varh=ifelse(Nh==0,NA,VYh),mean=TY/N,stderr=se/N,CV=se/TY,stratumID=stratumID,takeall=length(C),
              call=match.call(),date=date(),args=list(x=xgiven,n=n,CV=CV,Ls=Ls,certain=certain,alloc=alloc,rh=rh,model=model[1],model.control=model.control))
    class(out)<-"strata"
    return(out)
}
