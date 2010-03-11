`strata.cumrootf` <-
function(x,nclass,n=NULL,CV=NULL,Ls=3,certain=NULL,alloc=list(q1=0.5,q2=0,q3=0.5),rh=rep(1,Ls),model=c("none","loglinear","linear","random"),model.control=list())
{
    # Validation des arguments et initialisation de variables locales
    xgiven <- N <- findn <- L <- q1 <- q2 <- q3 <- beta <- sig2 <- ph <- pcertain <- gamma <- epsilon <- NULL
    takenone <- bias.penalty <- takeall.adjust <- A <- B <- C <- NULL
    out.check <- checkargs(x=x,n=n,CV=CV,Ls=Ls,certain=certain,alloc=alloc,rh=rh,model=model,model.control=model.control)
          for(i in 1:length(out.check)) assign(names(out.check)[i],out.check[[i]])
    if (!missing(nclass)) {   
          if (!((length(nclass)==1)&&isTRUE((nclass%%1)==0)&&isTRUE(nclass>=Ls))) stop("'nclass' must be an integer greater or equal to Ls")    
    } else {
          nclass<-min(L*10,length(x))
          on.exit(warning("'nclass' value has been chosen arbitrarily"))
    }  
    
    # Détermination des bornes
        # Découpage des classes
        cla<-x[1]+(x[length(x)]-x[1])*((0:nclass)/nclass) # L'étandue des données (min(x) à max(x)) est découpée en nclass classes de même longueur
        cla[nclass+1]<-cla[nclass+1]+1
        factor.c<-vector(length=length(x)) # Chaque donnée est associée à une classe
        for (i in 1:nclass) { factor.c[x>=cla[i]&x<cla[i+1]]<-i }
        # Calcul des cum sqrt(f)
        freq.c<-rep(0,nclass)
        pres.c<-tapply(x,factor.c,length) # Calcul du nombre de données par classe
        freq.c[as.numeric(names(pres.c))]<-as.vector(pres.c) # Les fréquences nulles sont conservées
        csfreq.c<-cumsum(sqrt(freq.c))
        # Calcul des sum sqrt(f) pour tous les regroupements potentiels (matrice L par 2^(L-1))
            # Regroupement potentiel = fréquences cumulées dans les strates juste inférieures ou supérieures à la fréquence cumulée cible.
            # Plus petite strate : 2 choix de bornes sup = la 1ière borne de classe pour laquelle sum sqrt(F) < 'but' 
            #                                              ou 1ière borne de classe pour laquelle sum sqrt(F) > 'but'.
            # 2e strate : connaissant la borne de la première strate, on a encore 2 choix de bornes.
            # Ainsi de suite jusqu'à la (L-1)e strate (la borne sup de la Le strate est connue = max(x)) : 2^(L-1) regroupements potentiels.
        but<-csfreq.c[nclass]/L # Fréquence cumulée à atteindre par classe
        nclass.temp<-vector(length=2^(L-1)-1) # Vecteur servant au calcul du nombre de classes par strate pour chaque regroupement potentiel
        nclassh<-ssfreqh<-matrix(0,L,2^(L-1))
        sous<-0 # Fréquences cumulées à soustraire (la longueur de 'sous' est 2^(h-1) pour la strate h)
        k<-1
        for (i in 1:(L-1)) # Les matrices 'nclassh' et 'ssfreqh' sont construites ligne par ligne, donc strate par strate.
        {
            k1<-k
            for(j in 1:length(sous)) # Pour la ligne h, on a besoin de calculer seulement 2^(h-1) valeurs, les autres valeurs peuvent être déduites.
            {
                a<-csfreq.c-sous[j] # On ramène à zéro les fréquences cumulées à partir du début de la strate pour laquelle le calcul est fait.
                b<-a[a>0&a<but] # On met de côté les classes qui forment la strate quand sum sqrt(f) < 'but'.
                nclass.temp[k]<-length(b) # Nombre de données dans la strate
                k<-k+1
            }
            nclassh[i,]<-rep(c(t(cbind(nclass.temp[k1:(k-1)],nclass.temp[k1:(k-1)]+1))),each=(2^(L-2))/(k-k1))
            cumnclass<-apply(nclassh[1:i,,drop=FALSE],2,sum)
            if (isTRUE(i==(L-1)))
            { # Certains regroupements impossibles peuvent avoir été considérés (sum(nclassh[,i])>nclass), ils sont ici retirés des calculs.
                nclassh<-nclassh[,cumnclass<=nclass]
                ssfreqh<-ssfreqh[,cumnclass<=nclass]
                cumnclass<-cumnclass[cumnclass<=nclass]
            }
            cumss<-if (isTRUE(i==1)) 0 else apply(ssfreqh[1:(i-1),,drop=FALSE],2,sum)
            ssfreqh[i,]<-ifelse(cumnclass==0,0,csfreq.c[cumnclass])-cumss
            pos<-seq(1,2^(L-1),(2^(L-2))/(k-k1)) # Détermination des colonnes pour lesquelles un calcul devra être fait à la boucle suivante.
            sous<-apply(ssfreqh,2,sum)[pos]
        }
        nclassh[L,]<-nclass-cumnclass # Déduction de la dernière ligne de nclassh (la somme de chaque colonne doit être nclass)
        ssfreqh[L,]<-csfreq.c[nclass]-apply(ssfreqh,2,sum) # Déduction de la dernière ligne de ssfreqh
        # Identification du meilleur regroupement et des bornes associées
            # Le regroupement choisi est celui avec la plus petite norme du vecteur des différences entre les sum sqrt(f) par strate et 'but'.
        nrgr<-order(apply((ssfreqh-but)^2,2,sum))[1]
        bhfull<-cla[cumsum(c(1,nclassh[,nrgr]))]
    
    # stratification
    Nh <- Nc <- EYh <- EYc <- VYh <- nh <- nh.nonint <- opti <- TY <- se <- bias <- prop <- stratumID <- NULL
    out.get <- bh2nh(xgiven,N,bhfull,findn,n,CV,L,certain,q1,q2,q3,A,B,C,bias.penalty,takeall.adjust,rh,model,beta,sig2,ph,pcertain,gamma,epsilon)
    for(i in 1:length(out.get)) assign(names(out.get)[i],out.get[[i]])
    # Note : se est la racine carrée de la variance de TY
    # on donne en sortie la moyenne = TY/N
    # donc RMSE est le root mean squared error de la moyenne, il vaut donc se/N
    # CV ou RRMSE est le relative root mean squared error de la moyenne ou de la somme : se/TY = (se/N)/(TY/N)

    # Sortie des résultats
    out<-list(Nh=Nh,nh=nh,n=sum(nh)+Nc,nh.nonint=nh.nonint,certain.info=c(Nc=Nc,meanc=EYc),opti.criteria=opti,bh=as.vector(bhfull[-c(1,L+1)]),
              meanh=ifelse(Nh==0,NA,EYh),varh=ifelse(Nh==0,NA,VYh),mean=TY/N,stderr=se/N,CV=se/TY,stratumID=stratumID,nclassh=nclassh[,nrgr],takeall=length(C),
              call=match.call(),date=date(),args=list(x=xgiven,nclass=nclass,n=n,CV=CV,Ls=Ls,certain=certain,alloc=alloc,rh=rh,model=model,model.control=model.control))
    class(out)<-"strata"
    return(out)
}
