`strata.cumrootf` <-
function(x,nclass,n=NULL,CV=NULL,Ls=3,certain=NULL,alloc=list(q1=0.5,q2=0,q3=0.5),rh=rep(1,Ls),model=c("none","loglinear","linear","random"),model.control=list())
{
    # Validation des arguments et initialisation de variables locales
    xgiven <- N <- N1 <- findn <- L <- q1 <- q2 <- q3 <- beta <- sig2 <- ph <- pcertain <- gamma <- epsilon <- NULL
    takenone <- bias.penalty <- takeall.adjust <- A <- B <- C <- NULL
    out.check <- checkargs(x=x,n=n,CV=CV,Ls=Ls,certain=certain,alloc=alloc,rh=rh,model=model,model.control=model.control)
          for(i in 1:length(out.check)) assign(names(out.check)[i],out.check[[i]])
    if (!missing(nclass)) {   
          if (!((length(nclass)==1)&&isTRUE((nclass%%1)==0)&&isTRUE(nclass>=Ls))) stop("'nclass' must be an integer greater or equal to Ls")    
    } else {
          nclass <- min(L*15, N1)
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
     but<-csfreq.c[nclass]/L # Fréquence cumulée cible à atteindre par classe
     nclass.temp<-vector(length=2^(L-1)-1) # Vecteur servant au calcul du nombre de classes par strate pour chaque regroupement potentiel
     csfreqh.temp<-rep(0,2^(L-1)) # Vecteur servant au calcul de la somme de sqrt(f) par strate pour chaque regroupement potentiel
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
       cumnclass<-colSums(nclassh[1:i,,drop=FALSE])
       cumss<-csfreqh.temp # ou cumss<-if (isTRUE(i==1)) 0 else colSums(ssfreqh[1:(i-1),,drop=FALSE])
       csfreqh.temp[cumnclass!=0]<-csfreq.c[cumnclass] # retourne NA si cumnclass>nclass
       csfreqh.temp[cumnclass==0]<-0 # cas particulier  
       ssfreqh[i,]<-csfreqh.temp-cumss
       pos<-seq(1,2^(L-1),(2^(L-2))/(k-k1)) # Détermination des colonnes pour lesquelles un calcul devra être fait à la boucle suivante.
       sous<-colSums(ssfreqh)[pos]
     }
     nclassh[L,]<-nclass-cumnclass # Déduction de la dernière ligne de nclassh (la somme de chaque colonne doit être nclass)
     ssfreqh[L,]<-csfreq.c[nclass]-sous # Déduction de la dernière ligne de ssfreqh (sous est la même chose que colSums(ssfreqh) car pos<-1:(2^(L-1))
     # Identification du meilleur regroupement et des bornes associées
       # Le regroupement choisi est celui avec la plus petite norme du vecteur des différences entre les sum sqrt(f) par strate et 'but'.        
       # On retire les regroupements contenant un nombre de classes nul ou négatif (cas particulier dernière strate) pour au moins une strate
     out<-apply(nclassh<=0,2,any)
     nrgr<-order(colSums((ssfreqh[,!out]-but)^2))[1]
     bhfull<-cla[cumsum(c(1,nclassh[,!out][,nrgr]))]
     
     ### Note 13 mai 2010 :
     ### J'ai corrigé un fonctionnement incorrect lorsque il n'y a qu'une seule borne possible pour la première borne (csfreq[1]>but).
     ### J'ai ensuite effectué plusieurs tests. Le programme fonctionne bien même si des regroupements potentiels contiennent des nombres
     ### de classes nuls pour des strates ou un nombre cumulatif de classes supérieur à nclass. Le programme fonctionne aussi bien même
     ### s'il doit couper à un endroit où des classes avec féquences nuls sont présentes. 
        
    
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
              meanh=ifelse(Nh==0,NA,EYh),varh=ifelse(Nh==0,NA,VYh),mean=TY/N,stderr=se/N,CV=se/TY,stratumID=stratumID,nclassh=nclassh[,nrgr],takeall=length(C),
              call=match.call(),date=date(),args=list(x=xgiven,nclass=nclass,n=n,CV=CV,Ls=Ls,certain=certain,alloc=alloc,rh=rh,model=model,model.control=model.control))
    class(out)<-"strata"
    return(out)
}
