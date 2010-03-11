`pbh2bh` <-
function(pbh,x1)
{ # makes the conversion from stratum boundaries expressed in terms of data rank (pbh), 
  # to stratum boundaries expressed on the scale of the data (bh).
    bh<-vector(length=length(pbh))
    for (i in 1:length(pbh))
        if (pbh[i]<=1) {
            bh[i]<-x1[1]
        } else if (pbh[i]>=length(x1)) {
            bh[i]<-(x1[length(x1)]+x1[length(x1)-1])/2
        } else bh[i]<-(x1[pbh[i]]+x1[pbh[i]-1])/2
    return(bh)
}   

testNh <- function(Nh,A,B,C,minNh) { # teste les conditions sur les Nh
     all(Nh[c(B,C)]>=minNh) && all(Nh[A]>=0)
     # retourne TRUE si tous les Nh des strates échantillonnées sont supérieurs ou égaux à minNh
     # et si tous les Nh des strates à tirage nul sont positifs ou nuls; retourne FALSE sinon
}

testnh <- function(nh,B,C) { # teste les conditions sur les nh
     all(nh[c(B,C)]>0)
     # retourne TRUE si tous les nh des strates échantillonnées sont positifs; retourne FALSE sinon
}

`checkargs` <- 
function(x,n,CV,Ls,certain,alloc,takenone,bias.penalty,takeall,takeall.adjust,rh,model,
     model.control,initbh,bh,algo,algo.control,y,rh.postcorr,testcertain=TRUE)
{
# Fonction qui permet de valider les arguments donnés en entrée à toutes les fonctions publiques du package
# Cette fonction retourne un objet nommé out qui contient certains arguments mis en forme et d'autres variables
# déduites des arguments.

out <- list()

if (!missing(Ls)) {
    if (!((length(Ls)==1)&&isTRUE((Ls%%1)==0)&&isTRUE(Ls>=2))) stop("'Ls' must be an interger greater or equal to 2")
    out$L <- Ls
}

if (!missing(x)) {
    out$xgiven <- x
    out$x <- sort(x)
    out$N <- length(x)
    if (!(is.vector(x)&&is.numeric(x))) stop("'x' must be a numeric vector")
    if (any(x<0)) stop("'x' must take non-negative values only")
    if (length(unique(x))<Ls) stop("it is impossible to form Ls strata containing at least one unit with the given 'x'")
}

if (!missing(n)) {
    if (!is.null(n)) { if (!((length(n)==1)&&isTRUE((n%%1)==0)&&isTRUE(n>0))) stop("'n' must be an integer greater than 0") }
} else out$n <- NULL

if (!missing(CV)) {
    if (!is.null(CV)) { if (!((length(CV)==1)&&is.numeric(CV))) stop("'CV' must be a numeric") }
} else out$CV <- NULL

if(!missing(n)&!missing(CV)) {
    if (is.null(n)&&is.null(CV)) stop("The argument 'n' or the argument 'CV' must be inputed")
    if (!is.null(n)&&!is.null(CV)) stop("Only one of the arguments 'n' and 'CV' can be inputed")
    out$findn <- if (is.null(n)&&!is.null(CV)) TRUE else FALSE
} else out$findn <- TRUE

if (!missing(certain)) {
    if(testcertain) {
        if (!is.null(certain)) {
            # changement de format si on n'a pas un vecteur
            out$certain <- c(certain)
            if (is.list(certain)) out$certain=unlist(out$certain)
            # validation principale
            if (!(all((out$certain%%1)==0)&&all(out$certain>0)&&all(out$certain<=out$N)))
                stop("'certain' must be a vector of integers between 1 and N, the length of 'x'")
            # exclusion des valeurs dans le vecteur x
            out$x <- sort(x[-certain]) 
            # validations pour éviter bugs
            if (!is.null(n)) if(length(unique(out$certain))>n-Ls)
                stop("'certain' must contain at most 'n'-'Ls' unique values") 
            if (length(unique(out$x))<Ls)
                stop("'certain' contains too much values: after removing units sampled with certainty, the x vector does not contain enough unique units to form 'Ls' strata") 
        } else out$certain <- NULL
    } else out$certain <- certain 
}

if (!missing(alloc)) {
    alloc<-as.list(alloc)
    if(!length(alloc)==3) stop("'alloc' must be a list of length 3")
    if(is.null(names(alloc))) names(alloc) <- c("q1","q2","q3")
    out$q1 <- if(is.null(alloc$q1)) 0.5 else alloc$q1
    out$q2 <- if(is.null(alloc$q2)) 0 else alloc$q2
    out$q3 <- if(is.null(alloc$q3)) 0.5 else alloc$q3
    if (!((length(out$q1)==1)&&is.numeric(out$q1)&&isTRUE(out$q1>=0))) stop("'alloc$q1' must be a numeric greater or equal to 0")
    if (!((length(out$q2)==1)&&is.numeric(out$q2)&&isTRUE(out$q2>=0))) stop("'alloc$q2' must be a numeric greater or equal to 0")
    if (!((length(out$q3)==1)&&is.numeric(out$q3)&&isTRUE(out$q3>=0))) stop("'alloc$q3' must be a numeric greater or equal to 0")
    out$alloc <- list(q1=out$q1,q2=out$q2,q3=out$q3)
}

if(!missing(takenone)) { ### dépend de Ls
    if (!((length(takenone)==1)&&(isTRUE(takenone==0)||isTRUE(takenone==1)))) stop("'takenone' must be 0 or 1")
    out$L <- Ls + takenone    
} else {
     takenone <- out$takenone <- 0
     out$bias.penalty <- 1
}

if(!missing(bias.penalty)) {
    if (!((length(bias.penalty)==1)&&is.numeric(bias.penalty)&&isTRUE(bias.penalty>=0)&&isTRUE(bias.penalty<=1)))
        stop("'bias.penalty' must be a numeric between 0 and 1 inclusively")
}
    
if(!missing(takeall)) { ### dépend de Ls
    if (!((length(takeall)==1)&&isTRUE((takeall%%1)==0)&&(isTRUE(takeall>=0)||isTRUE(takeall<=Ls-1))))
        stop("'takeall' must be an integer between 0 and 'Ls'-1 inclusively")
} else takeall <- 0 

if(!missing(takeall.adjust)) { ### dépend de Ls
    if (!is.logical(takeall.adjust)) stop("'takeall.adjust' must be a logical")
} else out$takeall.adjust <- TRUE

if (!missing(model)) {
    model.accepted=c("none","loglinear","linear","random")
    if (!(model[1]%in%model.accepted))  stop(paste("'model' must be one of the following character strings:",
          paste(model.accepted[-length(model.accepted)],collapse=", "),"or",model.accepted[length(model.accepted)]))
    out$model <- model[1]
}

if (!missing(model.control)) { ### dépend de certain, takenone et model
    if(!is.list(model.control)) stop("'model.control' must be a list")
    if( length(model.control)>0 && is.null(names(model.control)) ) stop("'The elements of the list 'model.control' must be named")
    out$beta <- if(is.null(model.control$beta)||out$model=="none") 1 else model.control$beta
    out$sig2 <- if(is.null(model.control$sig2)||out$model=="none") 0 else model.control$sig2
    out$ph <- if(is.null(model.control$ph)||out$model=="none") 1 else model.control$ph
    if (missing(takenone)) ptakenone <- NULL 
    else if (takenone==0) ptakenone <- NULL 
    else { ptakenone <- if(is.null(model.control$ptakenone)||out$model=="none") 1 else model.control$ptakenone } 
    if (missing(certain)) out$pcertain <- NULL 
    else if (is.null(out$certain)) out$pcertain <- NULL 
    else { out$pcertain <- if(is.null(model.control$pcertain)||out$model=="none") 1 else model.control$pcertain }
    out$gamma <- if(is.null(model.control$gamma)) 0 else model.control$gamma
    out$epsilon <- if(is.null(model.control$epsilon)) 0 else model.control$epsilon 
    if (!((length(out$beta)==1)&&is.numeric(out$beta))) stop("'model.control$beta' must be a numeric")
    if (!((length(out$sig2)==1)&&is.numeric(out$sig2)&&isTRUE(out$sig2>=0))) stop("'model.control$sig2' must be a numeric greater or equal to 0")
    if (!(((length(out$ph)==1)||(length(out$ph)==Ls))&&is.numeric(out$ph)&&all(out$ph>=0)&&all(out$ph<=1)))
        stop("'model.control$ph' must be a single numeric or a numeric vector of length Ls. Each element of 'model.control$ph' must be between 0 and 1 inclusively")
    if (!(is.null(ptakenone)||((length(ptakenone)==1)&&is.numeric(ptakenone)&&isTRUE(ptakenone>=0)&&isTRUE(ptakenone<=1))))
        stop("'model.control$ptakenone' must be a numeric between 0 and 1 inclusively")
    if (!(is.null(out$pcertain)||((length(out$pcertain)==1)&&is.numeric(out$pcertain)&&isTRUE(out$pcertain>=0)&&isTRUE(out$pcertain<=1))))
        stop("'model.control$pcertain' must be a numeric between 0 and 1 inclusively")
    if (!((length(out$gamma)==1)&&is.numeric(out$gamma)&&isTRUE(out$gamma>=0))) stop("'model.control$gamma' must be a numeric greater or equal to 0")
    if (!((length(out$epsilon)==1)&&is.numeric(out$epsilon)&&isTRUE(out$epsilon>=0))) stop("'model.control$epsilon' must be a numeric between 0 and 1 inclusively")
    if (isTRUE(length(out$ph)==1)) out$ph<-rep(out$ph,Ls)
    out$model.control <- if (identical(model[1],"none")) list() else
                         if (identical(model[1],"loglinear")) list(beta=out$beta,sig2=out$sig2,ph=out$ph) else
                         if (identical(model[1],"linear")) list(beta=out$beta,sig2=out$sig2,gamma=out$gamma) else
                         if (identical(model[1],"random")) list(epsilon=out$epsilon)
    if (identical(model[1],"loglinear")&&!is.null(ptakenone)) out$model.control <- c(out$model.control,ptakenone=ptakenone)  
    if (identical(model[1],"loglinear")&&!is.null(out$pcertain)) out$model.control <- c(out$model.control,pcertain=out$pcertain)  
    out$ph<-c(ptakenone,out$ph) # ph est de longueur L, ce qui est requis pour les calculs dans l'algo de Sethi
}

if(!missing(rh)) { ### dépend de ptakenone
    if (!(((length(rh)==1)||(length(rh)==Ls))&&is.numeric(rh)&&all(rh>=0)&&all(rh<=1)))
        stop("'rh' must be a single numeric or a numeric vector of length Ls. Each element of 'rh' must be between 0 and 1 inclusively")
    if (isTRUE(length(rh)==1)) out$rh<-rep(rh,Ls) else out$rh<-rh
    if (!is.null(ptakenone)) out$rh<-c(NA,out$rh)
}

if (!missing(initbh)) {
# En commentaires : méthode envisagée mais mise de côté
    out$initbh <-  if(is.null(initbh)) quantile(out$x,probs=(1:(Ls-1))/Ls) else initbh
#    out$initbh <- initbh
#    if(!is.null(out$initbh)) {
         if ((length(out$initbh)==Ls-1)&&(takenone==1)) {
               b1<-quantile(out$x,probs=0.01)
               if (b1==min(out$x)) b1<-unique(out$x)[2]
               out$initbh <- c(b1,out$initbh)
#               out$initbh <- c(min(out$x),out$initbh)
         }
         if (!((length(out$initbh)==Ls+takenone-1)&&is.numeric(out$initbh))) stop("'initbh' must be a numeric vector of length Ls+takenone-1 or Ls-1")
#    }
}

if (!missing(bh)) {
    if (!((length(bh)==Ls+takenone-1)&&is.numeric(bh))) stop("'bh' must be a numeric vector of length Ls+takenone-1")
}

if (!missing(algo)) {
    if (!(algo[1]%in%c("Kozak","Sethi"))) stop("'algo' must be the character string 'Kozak' or 'Sethi'")
    if("Sethi"==algo[1]) {
          arret1 <- if(missing(CV)) { TRUE } else { if(is.null(CV)) TRUE else FALSE }
          if(arret1) stop("To perform stratification minimizing RRMSE, please use Kozak's algorithm. In this package, Sethi's algorithm can only be used to perform stratification minimizing sample size.")
          if(out$model%in%c("linear","random")) stop("To take into account a 'linear' or 'random' model between X and Y, please use Kozak's algorithm. In this package, Sethi's algorithm can only be used with no model or with the loglinear model.")
    } 
    out$algo <- algo[1]
}

if (!missing(algo.control)) {
    if(!is.list(algo.control)) stop("'algo.control' must be a list")
    if(length(algo.control)>0&&is.null(names(algo.control))) stop("'The elements of the list 'algo.control' must be named")
    
    out$method <- if(is.null(algo.control$method)) "original" else algo.control$method
    if (!(out$method%in%c("original","modified"))) stop("'method' must be the character string 'original' or 'modified'")

    if(is.null(algo.control$maxiter)) {
        out$maxiter <- if( out$method=="modified") 3000 else if(out$method=="original") 10000 else 500
    } else out$maxiter <- algo.control$maxiter
    out$minNh <- if(is.null(algo.control$minNh)) 2 else algo.control$minNh
    out$maxstep <- if(is.null(algo.control$maxstep)) 3 else algo.control$maxstep
    out$maxstill <- if(is.null(algo.control$maxstill)) 100 else algo.control$maxstill
    if(is.null(algo.control$rep)) {
        out$rep <- if( out$method=="modified") 1 else 3
    } else out$rep <- algo.control$rep
    # utile seulement si je veux arrêter les répétitions après un certain nombre de répétitions identiques
    # out$stoprep <- if(is.null(algo.control$rep)) TRUE else FALSE
    out$minsol <- if(is.null(algo.control$minsol)) 1000 else algo.control$minsol
    
    if (!((length(out$minNh)==1)&&isTRUE((out$minNh%%1)==0)&&isTRUE(out$minNh>=2))) stop("'minNh' must be an integer greater or equal to 2")
    if (!((length(out$maxstep)==1)&&isTRUE((out$maxstep%%1)==0)&&(isTRUE(out$maxstep>=2)||isTRUE(out$maxstep<=10))))
        stop("'maxstep' must be an interger between 2 and 10 inclusively")
    if (!((length(out$maxiter)==1)&&isTRUE((out$maxiter%%1)==0)&&isTRUE(out$maxiter>0))) stop("'maxiter' must be a positive integer")
    if (!((length(out$maxstill)==1)&&isTRUE((out$maxstill%%1)==0)&&isTRUE(out$maxstill>0))) stop("'maxstill' must be a positive integer")
    if (is.character(out$rep)) {
          if ( !( length(out$rep)==1 && out$rep=="change" ) ) stop("'rep' must be a positive integer or 'change'") 
    } else { 
          if ( !( length(out$rep)==1 && isTRUE((out$rep%%1)==0) && isTRUE(out$rep>0) ) )
               stop("'rep' must be a positive integer or 'change'")
    }
    if (!((length(out$minsol)==1)&&isTRUE((out$minsol%%1)==0)&&isTRUE(out$minsol>=100)&&isTRUE(out$minsol<=1000000)))
          stop("'minsol' must be an integer between 100 and 1 000 000 inclusively")
}

if (!missing(y)) {
    if (!is.null(y)) { if (!(is.vector(y)&&is.numeric(y))) stop("'y' must be a numeric vector") }
}

if (!missing(rh.postcorr)) {
    if (!is.logical(rh.postcorr)) stop("'rh.postcorr' must be a logical")
}

### Détermination des ensembles de strates
out$A <- if(takenone==0) NULL else 1:takenone
out$B <- (takenone+1):(out$L-takeall)
out$C <- if(takeall==0) NULL else (out$L-takeall+1):out$L

return(out)

}


initbh.robust <- function(x1,N1,Ls,takenone,takeall,minNh,wtx1) {
     initbh <- vector(length=Ls-1)
     wtx1_copy <- wtx1
   ### strates takeall s'il y en a ##
     # On les veux les plus petites possibles, tout en respectant minNh, car elles peuvent causer un n négatif.
     # On ne soucie pas de pouvoir calculer une variance des ces strates (N1h>1 ou au moins 2 unités différentes)
     # car cette variance peut être nulle sans causer de problèmes dans les calculs subséquents.
     if (takeall!=0) {
          for (i in 1:takeall) {
               pos <- sum(cumsum(rev(wtx1_copy))<minNh)
               initbh[Ls-i] <- x1[N1-pos]
               N1 <- N1 - pos - 1
               wtx1_copy <- wtx1_copy[1:N1]
          }
     }
   ### strates takesome ###
     # On veut une répartition la plus égale possible en terme des N1h (nombres d'unités différentes dans les strates),
     # ce qui donne des strates avec des plus petites variances que de prendre une répartition la plus égale possible
     # en terme des Nh. Si les N1h ne peuvent être égaux, on ajoute une unité unique dans les premières strates
     # jusqu'à ce que sum(Nh)=N. Il faut aussi s'assurer que minNh soit respecté. S'il ne l'est pas, les strates avec N1h
     # environ égaux sont modifiées de la façon suivante :
     # On veut augmenter Nh dans la ou les strates avec un Nh<checkNh. Pour ce faire, on essaye tous les déplacements possibles
     # d'une valeur (représentant possiblement plus d'une unité) à partir de n'importe quelle strate vers la strate avec le plus 
     # petit Nh. Si une ou des solutions vérifient checkNh, on arrête là et on conserve celle qui donne la plus petite variance des Nh
     # car on les veut aussi les plus égaux possible. Si aucune solution ne vérifie checkNh, on reprend la procédure à partir de la 
     # solution qui donne la plus petite variance des Nh. En cas d'égalité des variances, on favorise des petites strates pour les 
     # grandes unités.      
     if (Ls-takeall>1) {
          N1hobj <- N1/(Ls-takeall)
          N1h <- rep(floor(N1hobj),Ls-takeall)
          plus <- N1-sum(N1h)
          if (plus>0) N1h[1:plus] <- N1h[1:plus]+1
          pos <- cumsum(N1h[-(Ls-takeall)])+1 
          # Il faut maintenant s'assurer que minNh est respecté
          Nh <- vector(length=Ls-takeall)
          for (i in 1:(Ls-takeall)) Nh[i] <- sum(wtx1[ifelse(i>1,pos[i-1],1):ifelse(i<Ls-takeall,pos[i]-1,N1)])
          checkNh <- all(Nh>=minNh)
          iter <- 0
          while (!checkNh && iter<2*minNh*Ls) {
               change  <- matrix(rep(c(pos,NA,NA),Ls-takeall-1),byrow=TRUE,nrow=Ls-takeall-1,ncol=Ls-takeall+1)
               # id colonnes de change :  pos (longueur Ls-takeall-1) + checkNh + varNh 
               posmin <- which.min(Nh)
               idrow <- 1
               for(j in rev((1:(Ls-takeall))[-posmin])) { # ici rev permet de favoriser des petites strates pour les grandes unités
                    if (posmin<j) {
                         change[idrow,posmin:(j-1)] <- pos[posmin:(j-1)] + 1
                    } else {
                         change[idrow,j:(posmin-1)] <- pos[j:(posmin-1)] - 1
                    }
                    for (i in 1:(Ls-takeall)) Nh[i] <- sum(wtx1[ifelse(i>1,change[idrow,i-1],1):ifelse(i<Ls-takeall,change[idrow,i]-1,N1)])
                    change[idrow,Ls-takeall] <- all(Nh>=minNh)
                    change[idrow,Ls-takeall+1] <- var(Nh)
                    idrow <- idrow + 1
               }
               checkNh <- sum(change[,Ls-takeall])>0
               if (checkNh) change <- change[as.logical(change[,Ls-takeall]),,drop=FALSE]
               keeprow <- which.min(change[,Ls-takeall+1])
               pos <- change[keeprow,1:(Ls-takeall-1)]
               iter <- iter + 1 # pour éviter les boucles infinies imprévues
          }
          initbh[1:(Ls-takeall-1)] <- x1[pos]
     }
   ### strate takenone s'il y en a ###
     # strate takenone initiale nulle car elle est la principale cause d'un n négatif.
     if (1==takenone) initbh <- c(min(x1),initbh)
     return(initbh)
}


anticip <- 
function(x,stratumID,model,beta,sig2,ph,pcertain,gamma,epsilon,A,L)
{ # Calcule les moments anticipés par strate en fonction des données, des bornes et du modèle
    Nh <- as.vector(table(stratumID))
    if (identical(model[1],"none")) {
        EYh <- tapply(x,stratumID,mean); EYh[is.na(EYh)] <- 0
        VYh <- tapply(x,stratumID,var)*(Nh-1)/Nh; VYh[is.na(VYh)] <- 0
    } else if (identical(model[1],"loglinear")) {
        ph <- c(ph,pcertain)
        phih <- tapply(x^beta,stratumID,sum)
        psih <- tapply(x^(2*beta),stratumID,sum)        
        EYh <- ifelse(Nh==0,0,ph*phih/Nh)
        VYh <- ifelse(Nh==0,0,ph*(exp(sig2)*psih/Nh-ph*(phih/Nh)^2))
    } else if (identical(model[1],"linear")) {
        EXh <- tapply(x,stratumID,mean)
        EX2h <- tapply(x^2,stratumID,mean)
        VXh <- EX2h-EXh^2
        EXgammah <- tapply(x^gamma,stratumID,mean)        
        EYh <- ifelse(Nh==0,0,beta*EXh)
        VYh <- ifelse(Nh==0,0,(beta^2)*VXh+sig2*EXgammah)
    } else if (identical(model[1],"random")) {
        EX <- mean(x)
        EX2 <- mean(x^2)
        EXh <- tapply(x,stratumID,mean)
        EX2h <- tapply(x^2,stratumID,mean)
        EYh <- ifelse(Nh==0,0,(1-epsilon)*EXh+epsilon*EX)
        VYh <- ifelse(Nh==0,0,(1-epsilon)*EX2h+epsilon*EX2-EYh^2)
    }
    VYh <- ifelse(VYh<0,0,VYh)
    TY <- sum(Nh*EYh)
    TAY <- sum(Nh[A]*EYh[A])
    if(length(Nh)==L) { Nc <- 0 ; EYc <- NA
    } else { Nc <- Nh[L+1] ; EYc <- EYh[L+1]; names(EYc) <- NULL } 
    # Pour enlever l'élément associé à la strate certain
    Nh <- Nh[1:L] ; EYh <- EYh[1:L] ; VYh <- VYh[1:L]
    if (model=="loglinear") { 
        phih <- phih[1:L] ; psih <- psih[1:L] 
        out<-list(Nh=Nh,Nc=Nc,phih=phih,psih=psih,EYh=EYh,EYc=EYc,VYh=VYh,TY=TY,TAY=TAY)
    } else { out<-list(Nh=Nh,Nc=Nc,EYh=EYh,EYc=EYc,VYh=VYh,TY=TY,TAY=TAY) }
    return(out)
}

RMSE <- 
function(bias.penalty,TAY,Nh,VYh,nh,rh,B,C,TY)
{ # Calcule le RMSE (se), relativebias (bias) et propbiasMSE (bias)
     se<- if (all(is.na(nh[B]))||any(nh[B]==0)) NA else sqrt((bias.penalty*TAY)^2+sum(((Nh^2)*VYh*(1/(nh*rh)-1/Nh))[c(B,C)]))
     prop<-(bias.penalty*TAY)^2/(se^2)
     bias<-(bias.penalty*TAY)/TY
     out <- list(se=se,prop=prop,bias=bias)
     return(out)
# Note : dans l'article et dans le programme survey estimator = somme de Y (Ty),
# alors que dans le package survey estimator = moyenne de Y (Ey).
# Lien entre les deux : Ey = Ty/N (voir définition de mean dans la sortie d'une fonction)
# RMSE de Ey = (RMSE de Ty)/N car Var(Ey) = Var(Ty/N) = Var(Ty)/N^2 et
#                                 biais(Ey) = Ty/N pop - Ty/N ech = Tay/N = biais(Ty)/N
# Définition RMSE : sqrt(biais^2+variance) donc RMSE de Ey = sqrt((biais(Ty)^2+Var(Ty))/N^2) =  sqrt((biais(Ty)^2+Var(Ty)))/N = (RMSE de Ty)/N
# Le facteur bias.penalty devant la biais est le même pour les deux RMSE
# RRMSE de Ey = RRMSE de Ty car RRMSE de Ey = RMSE de Ey / Ey = ((RMSE de Ty)/N)/(Ty/N) = RMSE de Ty / Ty = RRMSE de Ty
}

getstratum <-
function(xgiven,N,L,certain,bhfull)
{ # Pour chaque observation, identifie sa strate.
     stratumID <- vector(length=N)
     for (i in 1:L) stratumID[xgiven>=bhfull[i]&xgiven<bhfull[i+1]] <- i
     stratumID[certain] <- "certain"
     lev <- if(is.null(certain)) 1:L else c(1:L,"certain")
     stratumID <- factor(stratumID,levels=lev)
     return(stratumID)
}

getnh <-
function(L,A,B,C,n,findn,T1,ah,Nh)
{ # Calcule les nh selon la formule de l'article (en utilisant n (napprox ou ngiven))
     nh.nonint <- nh <- vector(length=L)
     nh.nonint[A] <- nh[A] <- 0
     nh.nonint[B] <- (n-T1)*ah
     nh.nonint[C] <- nh[C] <- Nh[C]
     # Arrondissement des nh dans les strates takesome
     if(all(is.na(nh.nonint[B]))) {
          nh[B] <- nh.nonint[B]
     } else {
          if (findn) { nh[B] <- ceiling(nh.nonint[B])
          } else { 
               # Si des nh sont < 1 mais supérieurs à 0, je les force à être ramenés à 1
               nh[B] <- ifelse(nh.nonint[B]>0&nh.nonint[B]<1,1,NA)
               Bnonint <- (1:L)[is.na(nh)] # position des chiffres à arrondir
               # nombre de chiffres à arrondir vers le haut et vers le bas
               nhaut <- n - sum(nh[A]) - sum(nh.nonint[B]>0&nh.nonint[B]<1) - sum(floor(nh.nonint[Bnonint])) - T1
               nbas <- length(Bnonint) - nhaut
               reste <- nh.nonint[Bnonint]%%1
               rankreste <- rank(reste,ties.method="first")
               idbas <- rankreste%in%(1:nbas)
               # Arrondissement vers le bas pour les nbas nombre avec les plus petites parties décimales 
               if(nbas>0) nh[Bnonint][idbas] <- floor(nh.nonint[Bnonint][idbas])
               # Arrondissement vers le haut pour les nhaut nombre avec les plus grandes parties décimales 
               if(nhaut>0) nh[Bnonint][!idbas] <- ceiling(nh.nonint[Bnonint][!idbas])
               # Dimminution d'unités si nhaut est un chiffre négatif
               if (nhaut < 0) nh[Bnonint][rankreste[-1*nhaut]] <- nh[Bnonint][rankreste[-1*nhaut]] - 1
#               nhaut <- n - sum(nh[A]) - sum(floor(nh.nonint[B])) - T1
#               nbas <- length(B) - nhaut
#               reste <- ifelse(nh.nonint[B]>0&nh.nonint[B]<1,1,nh.nonint[B]%%1)
#               # Si des nh sont < 1 mais supérieur à 0, je les force à être ramené à 1
#               rankreste <- rank(reste,ties.method="first")
#               idbas <- rankreste%in%(1:nbas)
#               # Arrondissement vers le bas pour les nbas nombre avec les plus petites parties décimales 
#               if(nbas>0) nh[B][idbas] <- floor(nh.nonint[B][idbas])
#               # Arrondissement vers le haut pour les nhaut nombre avec les plus grandes parties décimales 
#               if(nhaut>0) nh[B][!idbas] <- ceiling(nh.nonint[B][!idbas])
          }                                        
     }     
     out <- list(nh=nh,nh.nonint=nh.nonint)
     return(out)
}

getAlloc <-
function(q1,q2,q3,B,C,Nh,Nc,EYh,VYh) 
{ # Calcule ah et T1, les deux stat nécessaires pour faire l'allocation de n dans les strates
     gammah <- (Nh^2)^q1*(EYh^2)^q2*VYh^q3     
     ah <- if(sum(gammah[B])==0) NA else gammah[B]/sum(gammah[B])
     T1 <- sum(Nh[C])+ Nc
     out <- list(gammah=gammah,ah=ah,T1=T1)
     return(out)
}

optiCriteria <-
function(findn,n,CV,q1,q2,q3,bias.penalty,B,C,rh,Nh,VYh,TY,TAY,T1,ah,gammah)
{ # Calcule n (si CV donné) ou RRMSE (si n donné) selon les formules de l'article     
     U <- sum(((Nh^2)*VYh)[B]/(ah*rh[B]),na.rm=TRUE)
          U1 <- sum(((Nh^2)*VYh)[B]/(gammah[B]*rh[B]),na.rm=TRUE)
          U2 <- sum(gammah[B])
     V2 <- (bias.penalty*TAY)^2
     V3 <- sum((Nh*VYh)[B])
     V4 <- sum((Nh*VYh*(1-1/rh))[C])
     if (findn) {
        V1 <- (CV*TY)^2
        V <- V1-V2+V3+V4
        opti <- T1 + U/V      # n  
     } else {
        V <- NULL
        opti = sqrt(V2 + (U / (n - T1)) - V3 - V4) / TY     # RRMSE
     }
     out <- list(U=U,U1=U1,U2=U2,V=V,opti=opti)
     # U, U1, U2 et V sont seulement utilisés par l'algo de Sethi dans le calcul des dérivées
     return(out)
}

takeallauto <-
function(nh,Nh,L,A,B,C)
{ # Ajustement automatique pour les strates recensement : on change B et C
#     newcensus<-sum(ceiling(nh[B])>=Nh[B]|nh[B]<0)  # pourquoi le test nh[B]<0?
     newcensus <- if(all(is.na(nh[B]))) 0 else sum(nh[B]>Nh[B])
     if (newcensus>0&&length(C)<L-1-length(A)) {  
          B<-B[-length(B)]
          C<-c(L-length(C),C)
          valid <- FALSE
     } else valid<-TRUE
     out <- list(B=B,C=C,valid=valid)
     return(out)
# J'ai changé la condition nh[B])>=Nh[B] pour nh[B])>Nh[B] car on ne veux pas de dépassement. Si nh[B])=Nh[B]
# on atteint le CV cible ou le n cible, c'est correct. Ce changement règle certains cas problématiques que j'avais
# rencontré avec Kozak qui restait pris dans un minimum global à cause d'un ajustement pour strate takeall
# qui n'avait pas lieu d'être fait (dans premières strates par exemple).
}

bh2nh <- 
function(xgiven,N,bhfull,findn,n,CV,L,certain,q1,q2,q3,A,B,C,bias.penalty,takeall.adjust,rh,model,
         beta,sig2,ph,pcertain,gamma,epsilon)
{ # Calcule le plan de stratification, i.e. nh et même plus, à partir des bornes données 
  # et de tous les arguments donnés par l'utilisateur

    # Délimitation des strates
    stratumID <- getstratum(xgiven,N,L,certain,bhfull)

    # Calcul des moments anticipés
    moments<-anticip(xgiven,stratumID,model,beta,sig2,ph,pcertain,gamma,epsilon,A,L)
          Nh <- moments$Nh ; Nc <-  moments$Nc; EYh <- moments$EYh ; EYc <- moments$EYc 
          VYh <- moments$VYh ; TY <- moments$TY ; TAY <- moments$TAY 
     
    valid<-FALSE
    while(!valid) {
        # Calcul de ah, T1, n (si nécessaire) 
        out.Alloc <- getAlloc(q1,q2,q3,B,C,Nh,Nc,EYh,VYh)
            gammah <- out.Alloc$gammah ; ah <- out.Alloc$ah ; T1 <- out.Alloc$T1
        if (findn) {
            out.opti <- optiCriteria(findn,n,CV,q1,q2,q3,bias.penalty,B,C,rh,Nh,VYh,TY,TAY,T1,ah,gammah)
            n <- out.opti$opti
        }       
    
        # Détermination des nh
        out.nh <- getnh(L,A,B,C,n,findn,T1,ah,Nh)
            nh.nonint <- out.nh$nh.nonint; nh <- out.nh$nh
                  
        # Vérification nh<Nh sinon -> strates fixées strates-recensement une à une
        if (takeall.adjust) {
            out.adjust<-takeallauto(nh,Nh,L,A,B,C)
                    B <- out.adjust$B ;  C <- out.adjust$C ; valid <- out.adjust$valid
       } else valid<-TRUE
    }
    
    # Calcul du critère à optimiser (napprox ou RRMSEapprox)
    if (findn) {
        opti <- n
    } else {
        out.opti <- optiCriteria(findn,n,CV,q1,q2,q3,bias.penalty,B,C,rh,Nh,VYh,TY,TAY,T1,ah,gammah)
        opti <- out.opti$opti
    }  

    # Pour le calcul du RRMSE avec les nh entiers
    nh[B] <- pmax(pmin(nh[B],Nh[B]),0) # dernier ajustement , avant c'était pmax(pmin(nh[B],Nh[B]),1)
    if (all(is.na(nh[B]))) {
             warning("nh values for the takesome strata cannot be calculated because all takesome strata variances are null and q3 is nonzero")
    } else {
          if(all(nh[B]==0)) {
               if(all(VYh[B]==0)) { # Si CV cible + toutes les variances takesome sont nulles -> tous les nh takesome seront nuls
                    warning("nh values for the takesome strata are all null because all takesome strata variances are null, therefore the RRMSE cannot be calculated")
               } else if(all(nh.nonint[B]<0)) {
                    warning("nh.nonint values for the takesome strata are all negative and the corresponding nh have been rounded to zero, therefore the RRMSE cannot be calculated")
               } else {
                    warning("nh values for the takesome strata are all null, therefore the RRMSE cannot be calculated")
               }
          } else if (any(nh[B]==0)) {
               warning("some nh are null, possibly because of zero stratum variances and nonzero q3, therefore the RRMSE cannot be calculated")
          }   
    }
    out.RMSE <- RMSE(bias.penalty,TAY,Nh,VYh,nh,rh,B,C,TY)
          se <- out.RMSE$se ; prop <- out.RMSE$prop ; bias <- out.RMSE$bias
    
    # Sortie des résultats
    out<-list(Nh=Nh,Nc=Nc,EYh=EYh,EYc=EYc,VYh=VYh,nh=nh,nh.nonint=nh.nonint,TY=TY,se=se,bias=bias,prop=prop,opti=opti,stratumID=stratumID,C=C)
    return(out)
}


########################################################################################################################
## 15 février 2010
## Nouvelle approche : plus de checkNh et checknh car il n'est pas facile de créer
## des fonctions rapides qui modifient bien les bornes. Ça demanderait plus de temps
## et de réflexion pour bien faire ces fonctions et ce n'est vraiment pas une priorité.
## À partir de maintenant, si les bornes initiales ne vérifient pas les conditions, on les
## modifient pour les bornes par défaut qui sont les plus robustes. Si ces bornes ne rencontrent
## pas les conditions, aucunes bornes ne le peuvent. Si on a laissé tomber un argument initbh
## fournit par l'utilisateur afin d'utiliser les bornes par défaut, un message d'avis est affiché.

## Si jamais je reprend mon travail sur ces fonctions:
## Voici où en étaient mes réflexions :
## Les deux algo checkNh et checknh devraient avoir le même objectif : se repprocher le plus possible
## des bornes par défaut robustes. Ils seraient donc basés sur une même fonction, qui serait simplement
## roulée plus longtemps si checkNh roule, devient ok, mais checknh n'est toujours pas ok.
## À chaque étape je déplace un x1 de la strate la plus en surplus vers la strate la plus en déficit
## en favorisant la diminutions des strates takenone, puis takeall. Si plusieurs strates ont le même déficit,
## je prendrais le plus petit i car on préfère que les premières strates takesome soient plus grandes
## que les dernières. On prendrait la ou les unités (en cas d'égalités) à déplacer dans la strate 1- takenone,
## 2 - takeall, 3- la plus en surplus par rapport à son N1h cible. 

## idée : enregistrer l'historique + ajouter argument check.trace pour savoir si l'utilisateur
## veut voir ça dans la sortie. Pour moi, ce serait utile en cours de travail afin de mieux
## comprendre mon programme.

## idée : faire en C car ça risque d'être long à rouler.

########### inutilisé à partir du 15 fev 2010 #################
## Cette fonction était seulement utilisée pas checkNh et je n'ai plus besoin de cette
## fonction. En conséquence, getNh n'est plus utilisé.
`getNh` <-
function(pbh,L,wtx1,N1)
{ # simply calculates Nh, the number of units in each stratum, 
  # using the stratum boundaries expressed in terms of data rank (pbh)  
    resuC <- .C("getNhC",as.integer(pbh),as.integer(L),as.integer(wtx1),as.integer(N1),Nh=as.integer(rep(0,L)),PACKAGE="stratification")
    return(resuC$Nh)
  # code possible en R : beaucoup plus long!
    #fac <- c(rep(1:(L-1),pbh-c(1,pbh[-(L-1)])),rep(L,N1-pbh[L-1]+1))
    #fac <- factor(fac,levels=1:L)
    #Nh <- tapply(wtx1,list(stratum=fac),sum)

}

########### inutilisé à partir du 15 fev 2010 #################
## problème : les bornes cibles ne sont pas infaillibles ##
`checkNh` <-
function(x1,wtx1,bh,Ls,minNh,takenone=0)
{ # Verifies that every strata, except the take-none stratum if present, contains at least minNh units.
  # If not, the boundaries are modified in order to respect this condition (that's the real purpose of the function).
  # This function is used only on the initial strata boundaries.
  # Inside the C code for Kozak's algorithm, the verification Nh>=minNh is done for each set of strata boundaries.
  # In the algorithm, the boundaries are not changed if they don't respect the condition. They are simlpy discarded.  
    L <- Ls + takenone
    N1 <- length(x1)
    pbh <- vector(length=L-1)
    for (i in 1:(L-1)) {
        dif <- x1-bh[i]
        pbh[i] <- length(dif[dif<0])+1
    }
    A<- if(takenone==0) NULL else 1:takenone
    BC<-(takenone+1):L
    Nh<-getNh(pbh=pbh,L=L,wtx1=wtx1,N1=N1)
    detail<-NULL
    # Modification des Nh au besoin
    if (any(Nh[A]<0)|any(Nh[BC]<minNh)) {
        warn<-TRUE
        # Bornes cibles : répartition la plus possible égale 
#        bha<-quantile(x,probs=(1:(L-1))/L)
        bha<-vector(length=L-1)
        Nhcible <- sum(wtx1)/L
        revwtx1 <- rev(wtx1) # commencer par les dernière strates car on laisse ce qu'on veut dans la strate takenone
        nremove <- 0
        for (i in (L-1):1) {
            cumul<-cumsum(revwtx1)
            revpos<-sum(cumul<=Nhcible)
            if(0==revpos || sum(cumul[1:revpos])<minNh) revpos <- revpos + 1
            revwtx1 <- revwtx1[-(1:revpos)]
            nremove <- nremove + revpos
            pos <- length(x1) - nremove
            bha[i] <- (x1[pos]+x1[pos+1])/2
        }
        pbha<-vector(length=L-1)
        for (i in 1:(L-1)) {
            dif<-x1-bha[i]
            pbha[i]<-length(dif[dif<0])+1
        }
        Nha<-getNh(pbh=pbha,L=L,wtx1=wtx1,N1=N1)
        detail<-matrix(c(bha,Nha,bh,Nh),nrow=2,byrow=TRUE)
        while (any(Nh[A]<0)|any(Nh[BC]<minNh)){
            # Identification des bornes à modifier
            jbh<-0
            # strates takenone
            if (isTRUE(takenone>0)) if (Nh[1]<0) jbh<-c(jbh,1) 
            if (isTRUE(takenone>1)) for (j in 2:A[length(A)]) if (Nh[j]<0) jbh<-c(jbh,j-1,j)
            # strates takesome ou takeall
            if (isTRUE(Nh[BC[1]]<minNh)) if (isTRUE(takenone>0)) jbh<-c(jbh,BC[1]-1,BC[1]) else jbh<-c(jbh,1)
            if (isTRUE(L>(1+takenone))) if (Nh[L]<minNh) jbh<-c(jbh,L-1) 
            if (isTRUE(L>(2+takenone))) for (j in BC[2]:(L-1)) if (Nh[j]<minNh) jbh<-c(jbh,j-1,j)
            tab <- table(jbh[-1])
            jbh <- as.numeric(names(tab))
            # Sélection de la borne qui sera modifiée, i.e. celle qui est le plus loin de sa borne cible
            # et modification pour la rapprocher de la borne cible.
            jmodif <- jbh[order((pbh-pbha)[jbh]^2,decreasing=TRUE)][1]
            pbh[jmodif] <- if (pbh[jmodif]<pbha[jmodif]) pbh[jmodif]+1 else pbh[jmodif]-1
            Nh<-getNh(pbh=pbh,L=L,wtx1=wtx1,N1=N1)
            bh<-pbh2bh(pbh,x1)
            detail<-rbind(detail,c(bh,Nh))
        }
        colnames(detail) <- c(paste("b",1:(L-1),sep=""),paste("N",1:L,sep=""))
        rownames(detail) <- c("target",0:(nrow(detail)-2))
    } else warn<-FALSE
    out<-list(bh=bh,pbh=pbh,warn=warn,detail=detail)
}


########### inutilisé à partir du 15 fev 2010 #################
## problème : peut proposer des strates avec nh=0 ##
## raison : on vide takenone et on rend takeall les plus petites possibles
## mais on arrête là et ce n'est parfois pas suffisant pour s'assurer d'avoir
## des nh tous >0. Il faudrait aussi modifier les bornes vers des bornes cibles
## robustes, un peu comme le fait checkNh
`checknh` <-
function(x,x1,N,pbh,findn,n,CV,Ls,Nc,EYc,alloc,A,B,C,bias.penalty,rh,model,nmodel,model.control,minNh)
{ # Verifies that every sampled strata sample size is positive. 
  # If not, the boundaries are modified in order to respect this condition (that's the real purpose of the function).
  # This function is used only on the initial strata boundaries.
  # Inside the C code for Kozak's algorithms, the verification nh>0 is done for each set of strata boundaries.
  # In the algorithm, the boundaries are not changed if they don't respect the condition. They are simlpy discarded.  
    bh<-pbh2bh(pbh,x1)
    calculn<-strata.internal(x=x,N=N,bh=bh,findn=findn,n=n,CV=CV,Ls=Ls,Nc=Nc,EYc=EYc,alloc=alloc,takenone=length(A),
                             bias.penalty=bias.penalty,takeall=length(C),rh=rh,model=model,model.control=model.control) 
    L <- Ls + length(A) 
    # Modification des strates initiales afin d'avoir nh>0 pour toute strate échantillonnée
    if (any(calculn$nh[c(B,C)]<=0)) {
        warn<-TRUE
        # Programmé en C pour que ce soit plus rapide
        resuC<-.C("checknhC",as.double(x),as.double(x1),as.integer(N),as.integer(findn),as.integer(n),as.double(CV),
                  as.integer(L),as.integer(Nc),as.double(EYc),as.double(alloc$q1),as.double(alloc$q2),as.double(alloc$q3),
                  as.integer(length(A)),as.double(bias.penalty),as.integer(length(C)),as.double(rh),
                  as.integer(nmodel),as.double(model.control$beta),as.double(model.control$sig2),
                  as.double(model.control$ph),as.double(model.control$gamma),as.double(model.control$epsilon),
                  as.integer(minNh),as.integer(calculn$Nh),pbh=as.integer(pbh),PACKAGE="stratification")
        if (all(pbh==resuC$pbh)) stop("some initial nh for sampled strata are <=0 but the initial boundaries cannot be modified automatically: please change the 'initbh' argument")
        pbh<-resuC$pbh
    } else warn<-FALSE
    out<-list(pbh=pbh,warn=warn)
}

########################################################################################################################
