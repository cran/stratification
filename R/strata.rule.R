strata.geo <- function(x, n=NULL, CV=NULL, Ls=3, certain=NULL, alloc=list(q1=0.5,q2=0,q3=0.5), 
                           rh=rep(1,Ls), model=c("none","loglinear","linear","random"), model.control=list())
{
  ### Fonction externe : voir fiche d'aide pour la documentation
  
  # Validation des arguments et initialisation de variables :
  call.ext <- match.call()
  out <- valid_args(obj_fct = as.list(environment()), call.ext)
  # Variables generales
  N <- out$N; findn <- out$findn; L <- out$L; rhL <- out$rhL;
  # Variables relatives a la strate certain
  certain <- out$certain; xnoc <- out$xnoc; Nc <- out$Nc; Nnoc <- out$Nnoc;
  # Variables relatives a l'allocation
  q1 <- out$q1; q2 <- out$q2; q3 <- out$q3;
  # Variables relatives au model
  nmodel <- out$nmodel; beta <- out$beta; sig2 <- out$sig2; ph <- out$ph; pcertain <- out$pcertain; 
  gamma <- out$gamma; epsilon <- out$epsilon;
  # Variable pour la sortie : liste des arguments
  args <- out$args;
  # Variable initialise car takenone n'est pas un argument en entree ici
  takenone  <- out$takenone;
  
  # Verification unique a la methode geometrique
  if (any(x==0)) stop("the geometric method accepts only positive 'x' values") 
  
  # Initialisation de quelques variables supplementaires 
  out <- init_stat(obj_fct = as.list(environment()))
  EX <- out$EX;  EX2 <- out$EX2; EYc <- out$EYc;
  bias.penalty <- 1  ## a initialiser car pas arg et pas cree par valid_args car inutile aux autres validations
  
  # Determination des bornes pleines
  bhfull <- strata.geo.internal(obj_fct = as.list(environment()))
    
  # Stratification selon ces bornes
  out <- strata.bh.internal(bhfull = bhfull, takeallin = 0, takeall.adjust = TRUE, 
                            obj_fct = as.list(environment()))
  
  # Pour la sortie, je dois modifier un peu la liste. 
  out_rule(out, bhfull)
}


strata.geo.internal <- function(obj_fct)
{
  # Determination des bornes selon la regle geometrique de Gunning et Horgan (2004)
  # La sortie est les bornes des strates, incluant b0 et bL (vecteur de longueur L + 1).

  # Pour tirer de obj_fct les variables dont on a besoin ici
  xnoc <- obj_fct$xnoc
  L <- obj_fct$L

  # Calculs
  a <- min(xnoc)
  r <- (max(xnoc)/a)^(1/L)
  bhfull <- a*r^(0:L)
  bhfull[L+1] <- bhfull[L+1] + 1
  bhfull             
}


######################################################################################################

strata.cumrootf <- function(x, n=NULL, CV=NULL, Ls=3, certain=NULL, alloc=list(q1=0.5,q2=0,q3=0.5), 
                   rh=rep(1,Ls), model=c("none","loglinear","linear","random"), model.control=list(), nclass=NULL)
{
  ### Fonction externe : voir fiche d'aide pour la documentation
  
  # Validation des arguments et initialisation de variables :
  call.ext <- match.call()
  out <- valid_args(obj_fct = as.list(environment()), call.ext)
  # Variables generales
  N <- out$N; findn <- out$findn; L <- out$L; rhL <- out$rhL;
  # Variables relatives a la strate certain
  certain <- out$certain; xnoc <- out$xnoc; Nc <- out$Nc; Nnoc <- out$Nnoc;
  # Variables relatives a l'allocation
  q1 <- out$q1; q2 <- out$q2; q3 <- out$q3;
  # Variables relatives au model
  nmodel <- out$nmodel; beta <- out$beta; sig2 <- out$sig2; ph <- out$ph; pcertain <- out$pcertain; 
  gamma <- out$gamma; epsilon <- out$epsilon;
  # Variable pour la sortie : liste des arguments
  args <- out$args;
  # Variable initialise car takenone n'est pas un argument en entree ici
  takenone  <- out$takenone; 
  # Variable nclass a laquelle on a peut-etre du attribuer une valeur par defaut
  nclass  <- out$nclass; 
  
  # Initialisation de quelques variables supplementaires
  out <- init_stat(obj_fct = as.list(environment()))
  EX <- out$EX;  EX2 <- out$EX2; EYc <- out$EYc;
  bias.penalty <- 1  ## a initialiser car pas arg et pas cree par valid_args car inutile aux autres validations
  
  # Determination des bornes pleines
  out <- strata.cumrootf.internal(obj_fct = as.list(environment()))
  bhfull <- out$bhfull; nclassh <- out$nclassh;  

  # Stratification selon ces bornes
  out <- strata.bh.internal(bhfull = bhfull, takeallin = 0, takeall.adjust = TRUE, 
                            obj_fct = as.list(environment()))
  
  # Pour la sortie, je dois modifier un peu la liste. 
  out_rule(out, bhfull, nclassh)
}


strata.cumrootf.internal <- function(obj_fct)
{
  # Determination des bornes selon la regle du cumulative foot frequency de Dalenius et Hodges (1959)
  # La sortie est les bornes des strates, incluant b0 et bL (vecteur de longueur L + 1).

  # Pour tirer de obj_fct les variables dont on a besoin ici
  xnoc <- obj_fct$xnoc
  nclass <- obj_fct$nclass
  L <- obj_fct$L
  
  # Decoupage des classes
  cla<-min(xnoc)+(max(xnoc)-min(xnoc))*((0:nclass)/nclass) 
    ## L'etandue des donnees (min(x) a max(x)) est decoupee en nclass classes de meme longueur
  cla[nclass+1]<-cla[nclass+1]+1
  factor.c<-vector(length=length(xnoc)) # Chaque donnee est associee a une classe
  for (i in 1:nclass) { factor.c[xnoc>=cla[i]&xnoc<cla[i+1]]<-i }
  # Calcul des cum sqrt(f)
  freq.c<-rep(0,nclass)
  pres.c<-tapply(xnoc,factor.c,length) # Calcul du nombre de donnees par classe
  freq.c[as.numeric(names(pres.c))]<-as.vector(pres.c) # Les frequences nulles sont conservees
  csfreq.c<-cumsum(sqrt(freq.c))
  # Calcul des sum sqrt(f) pour tous les regroupements potentiels (matrice L par 2^(L-1))
  # Regroupement potentiel = frequences cumulees dans les strates juste inferieures ou superieures a la frequence cumulee cible.
  # Plus petite strate : 2 choix de bornes sup = la 1iere borne de classe pour laquelle sum sqrt(F) < 'but' 
  #                                              ou 1iere borne de classe pour laquelle sum sqrt(F) > 'but'.
  # 2e strate : connaissant la borne de la premiere strate, on a encore 2 choix de bornes.
  # Ainsi de suite jusqu'a la (L-1)e strate (la borne sup de la Le strate est connue = max(x)) : 2^(L-1) regroupements potentiels.
  but<-csfreq.c[nclass]/L # Frequence cumulee cible a atteindre par classe
  nclass.temp<-vector(length=2^(L-1)-1) # Vecteur servant au calcul du nombre de classes par strate pour chaque regroupement potentiel
  csfreqh.temp<-rep(0,2^(L-1)) # Vecteur servant au calcul de la somme de sqrt(f) par strate pour chaque regroupement potentiel
  nclassh<-ssfreqh<-matrix(0,L,2^(L-1))
  sous<-0 # Frequences cumulees a soustraire (la longueur de 'sous' est 2^(h-1) pour la strate h)
  k<-1
  for (i in 1:(L-1)) # Les matrices 'nclassh' et 'ssfreqh' sont construites ligne par ligne, donc strate par strate.
  {
    k1<-k
    for(j in 1:length(sous)) # Pour la ligne h, on a besoin de calculer seulement 2^(h-1) valeurs, les autres valeurs peuvent etre deduites.
    {
      a<-csfreq.c-sous[j] # On ramene a zero les frequences cumulees a partir du debut de la strate pour laquelle le calcul est fait.
      b<-a[a>0&a<but] # On met de cote les classes qui forment la strate quand sum sqrt(f) < 'but'.
      nclass.temp[k]<-length(b) # Nombre de donnees dans la strate
      k<-k+1
    }
    nclassh[i,]<-rep(c(t(cbind(nclass.temp[k1:(k-1)],nclass.temp[k1:(k-1)]+1))),each=(2^(L-2))/(k-k1))
    cumnclass<-colSums(nclassh[1:i,,drop=FALSE])
    cumss<-csfreqh.temp # ou cumss<-if (isTRUE(i==1)) 0 else colSums(ssfreqh[1:(i-1),,drop=FALSE])
    csfreqh.temp[cumnclass!=0]<-csfreq.c[cumnclass] # retourne NA si cumnclass>nclass
    csfreqh.temp[cumnclass==0]<-0 # cas particulier  
    ssfreqh[i,]<-csfreqh.temp-cumss
    pos<-seq(1,2^(L-1),(2^(L-2))/(k-k1)) # Determination des colonnes pour lesquelles un calcul devra etre fait a la boucle suivante.
    sous<-colSums(ssfreqh)[pos]
  }
  nclassh[L,]<-nclass-cumnclass # Deduction de la derniere ligne de nclassh (la somme de chaque colonne doit etre nclass)
  ssfreqh[L,]<-csfreq.c[nclass]-sous # Deduction de la derniere ligne de ssfreqh (sous est la meme chose que colSums(ssfreqh) car pos<-1:(2^(L-1))
  # Identification du meilleur regroupement et des bornes associees
  # Le regroupement choisi est celui avec la plus petite norme du vecteur des differences entre les sum sqrt(f) par strate et 'but'.        
  # On retire les regroupements contenant un nombre de classes nul ou negatif (cas particulier derniere strate) pour au moins une strate
  out<-apply(nclassh<=0,2,any)
  nrgr<-order(colSums((ssfreqh[, !out, drop=FALSE]-but)^2))[1]
  bhfull<-cla[cumsum(c(1,nclassh[, !out, drop=FALSE][,nrgr]))]
  
  list(bhfull=bhfull, nclassh=nclassh[, !out, drop=FALSE][,nrgr])

  ### Note 13 mai 2010 :
  ### J'ai corrige un fonctionnement incorrect lorsque il n'y a qu'une seule borne possible pour la premiere borne (csfreq[1]>but).
  ### J'ai ensuite effectue plusieurs tests. Le programme fonctionne bien meme si des regroupements potentiels contiennent des nombres
  ### de classes nuls pour des strates ou un nombre cumulatif de classes superieur a nclass. Le programme fonctionne aussi bien meme
  ### s'il doit couper a un endroit ou des classes avec fequences nuls sont presentes. 
}

