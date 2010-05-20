`strata.LH` <-
function(x,initbh=NULL,n=NULL,CV=NULL,Ls=3,certain=NULL,alloc=list(q1=0.5,q2=0,q3=0.5),takenone=0,bias.penalty=1,takeall=0,rh=rep(1,Ls),
         model=c("none","loglinear","linear","random"),model.control=list(),algo=c("Kozak","Sethi"),algo.control=list())
{
    bhi_type <- if (is.null(initbh)) "default" else "initbh"  
    # Validation des arguments et initialisation de variables locales
    xgiven <- N <- findn <-  L <- q1 <- q2 <- q3 <- beta <- sig2 <-  ph <- pcertain <- gamma <- epsilon <- NULL
    maxiter <- method <- minNh <- maxstep <- maxstill <- rep <- minsol <- A <- B <- C <- NULL
    out.check <- checkargs(x=x,initbh=initbh,n=n,CV=CV,Ls=Ls,certain=certain,alloc=alloc,takenone=takenone,bias.penalty=bias.penalty,
         takeall=takeall,rh=rh,model=model,model.control=model.control,algo=algo,algo.control=algo.control)
    for(i in 1:length(out.check)) assign(names(out.check)[i],out.check[[i]])
    
    Ai <- A;  Bi <- B;  Ci <- C; # regroupements initiaux à conserver car je les reprends à chaque nouvelle répétition de l'algo
    desc.iter<-NULL
    tab <- table(x)
    x1 <- as.numeric(names(tab))
    wtx1 <- as.numeric(tab)
    N1 <- length(wtx1)

    # test : est-il possible de former Ls strates ayant au moins minNh unités?
#    wtx1_copy <- if(0==takenone) wtx1 else wtx1[-1] # on veut au moins une unité dans la strate recensement au départ
    wtx1_copy <- wtx1
    for(i in 1:Ls) {
        nx1h <- sum(cumsum(wtx1_copy)<minNh)
        if (length(wtx1_copy)==nx1h) {
            stop("it is impossible to form Ls strata containing at least 'minNh'=",minNh," units with the given 'x' and 'certain'")
        } else wtx1_copy <- wtx1_copy[-(1:(nx1h+1))]
    }

    # test : si q3 est différent de 0, est-il possible de calculer une variance dans chaque strate à tirage partiel?
    if (q3!=0 && N1<Ls*2-takeall) stop("with the given 'x' and 'certain', some takesome strata variances are sure to be null, therefore it is impossible to form Ls strata having all positive nh with a nonzero 'q3'")
    
#############################################################################################

## Kozak
if (algo=="Kozak")
{
    rhC <- ifelse(is.na(rh),1,rh) # car on ne peut pas passer un NA à une fonction C, même si la variable n'est pas utilisé dans les calculs

    # Calculs de Nc (strate certain) et EYc
    # Note : TY varie en fonction des strates pour le modèle loglinéaire avec 
    # des taux de mortalité différents entre les strates.
    stratumID <- getstratum(xgiven,N,L=1,certain,bhfull=c(min(x),max(x)+1))
    moments<-anticip(xgiven,stratumID,model,beta,sig2,ph=1,pcertain,gamma,epsilon,A=NULL,L=1)
        Nc <-  moments$Nc; EYc <-  moments$EYc
    if(0==Nc) EYc <- 0; # car on ne peut pas passer NA à une fonction C

    # Essayer toutes les solutions possibles si elles sont en nombre inférieur à minsol
    nsol <- if(0==takenone) choose(N1-1,L-1) else choose(N1,L-1)
    if(nsol<minsol) {
        warnsol <- TRUE
        sol<-combn((2-takenone):N1,L-1)
        sol.detail <- matrix(NA,nrow=nsol,ncol=3*L+3)
        colnames(sol.detail)<-c(paste("b",1:(L-1),sep=""),paste("N",1:L,sep=""),paste("n",1:L,sep=""),"opti","takeall","Nhok","nhok")
        for (i in 1:nsol) {
             sol.detail[i,1:(L-1)] <- x1[sol[,i]]
             # J'ai pensé à d'abord calculer Nh et N1h et ne pas faire les calculs subséquents si Nhok ou nhok est faux.
             # Par contre, même en utilisant getNhC, ça ralentit la fonction, probablement parce que certains calculs se trouvent 
             # à être faits en double. Étant donné que le but était d'accélérer la fonction, je laisse tomber cette idée.
             A <- Ai; B <- Bi; C <- Ci
             valid<-FALSE
             while(!valid) {
                 calculn<-strata.internal(x=x,N=length(x),bh=x1[sol[,i]],findn=findn,n=n,CV=CV,Ls=Ls,
                                          Nc=Nc,EYc=EYc,alloc=list(q1=q1,q2=q2,q3=q3),takenone=length(A),
                                          bias.penalty=bias.penalty,takeall=length(C),rh=rhC,model=model,
                                          model.control=list(beta=beta,sig2=sig2,ph=ph,gamma=gamma,epsilon=epsilon))
                 # Vérification nh<Nh sinon -> strates fixées strates-recensement une à une
                 out.adjust<-takeallauto(nh=calculn$nh,Nh=calculn$Nh,L,A,B,C)
                      B <- out.adjust$B ;  C <- out.adjust$C ; valid <- out.adjust$valid
             }
             sol.detail[i,L:(2*L-1)] <- calculn$Nh
             sol.detail[i,(2*L):(3*L-1)] <- ifelse(-2147483647==calculn$nh,NA,calculn$nh)
             sol.detail[i,"opti"] <- calculn$opti
             sol.detail[i,"takeall"] <- length(C)
             sol.detail[i,"Nhok"] <- testNh(calculn$Nh,A,B,C,minNh)
             sol.detail[i,"nhok"] <- testnh(calculn$nh,B,C)
        }
        sol.detail.check <- sol.detail[as.logical(sol.detail[,"Nhok"])&as.logical(sol.detail[,"nhok"]),,drop=FALSE]
        if (nrow(sol.detail.check)==0) stop("it is impossible to form Ls sampled strata containing at least 'minNh' units and having all positive nh with the given 'x', 'certain' and 'alloc'")
        solmin <- which.min(sol.detail.check[,"opti"])
    } else {

    # Sinon faire rouler l'aglo de Kozak
        warnsol <- FALSE
        nmodel <- if (identical(model,"none")) 0 else if (identical(model,"loglinear")) 1 else 
                  if (identical(model,"linear")) 2 else if (identical(model,"random")) 3
        # desc.rep (data.frame) va contenir des informations pour chaque répétition de l'algorithme
        desc.rep <- data.frame(matrix(NA,nrow=ifelse("change"==rep,27,rep),ncol=3*(L-1)+10))
        dimnames(desc.rep)<-list(1:nrow(desc.rep),c("opti","takeall",paste("b",1:(L-1),sep=""),"niter","maxstep","maxstill","initbh.type",
                                 paste("initb",1:(L-1),sep=""),"Nhok.initbh","nhok.initbh",paste("b",1:(L-1),"i",sep=""),"Nhok.robust","nhok.robust"))
        # bhi.robust (matrice numérique) va contenir les bornes initiales par défaut pour chaque valeur de takeall utilisée 
        # + opti dont on a besoin pour l'algo de Kozak + résultats des tests sur Nh et nh.
        # J'enregistre ces résultats afin de ne pas avoir besoin de les obtenir inutilement plus d'une fois.
        bhi.robust <- matrix(NA,nrow=Ls-takeall,ncol=L+2,dimnames=list(takeall:(Ls-1),c(paste("b",1:(L-1),"i",sep=""),"Nhok","nhok","opti")))
        if ("change"==rep) {
            bhi_type <- c("default","geo","cumrootf")
            maxstep_val <- c(3,min(10,N1),min(ceiling(N1/10),100))
            repid <- 3 # nombre de répétitions identiques (mêmes paramètres)
        } else {
            bhi_type <- bhi_type
            maxstep_val <- maxstep
            repid <- rep
        }
        # On fait rouler l'algo de Kozak un certain nombre de fois selon la valeur de 'rep'
        rowrep <- 0
        for (idbhi in bhi_type)
        {           
            # calcul de bh, Nhok, nhok et opti pour takeall donné en entrée + sauvegarde
            if ("geo"==idbhi) {
                 bh <- suppressWarnings(strata.geo(x=x,n=n,CV=CV,Ls=Ls,certain=certain,alloc=list(q1=q1,q2=q2,q3=q3),rh=rh,model=model,
                        model.control=list(beta=beta,sig2=sig2,ph=ph,gamma=gamma,epsilon=epsilon))$bh)
            } else if("cumrootf"==idbhi) {
                 bh <- suppressWarnings(strata.cumrootf(x=x,nclass=min(L*10,length(x)),n=n,CV=CV,Ls=Ls,certain=certain,
                       alloc=list(q1=q1,q2=q2,q3=q3),rh=rh,model=model,model.control=list(beta=beta,sig2=sig2,ph=ph,gamma=gamma,epsilon=epsilon))$bh)
            } else bh <- initbh
            # bhi (matrice numérique) va contenir la valeur de bornes demandées + opti (pour algo Kozak)
            # + résultats des tests sur Nh et nh pour chaque valeur de takeall utilisée.
            bhi <- rep(NA,L+2*(Ls-takeall))
            names(bhi) <- c(paste("initb",1:(L-1),sep=""),"Nhok",paste("nhok",takeall:(Ls-1),sep=""),paste("opti",takeall:(Ls-1),sep=""))
            bhi[1:(L-1)] <- bh
            calculn<-strata.internal(x=x,N=length(x),bh=bh,findn=findn,n=n,CV=CV,Ls=Ls,Nc=Nc,EYc=EYc,
                                     alloc=list(q1=q1,q2=q2,q3=q3),takenone=length(Ai),
                                     bias.penalty=bias.penalty,takeall=length(Ci),rh=rhC,model=model,
                                     model.control=list(beta=beta,sig2=sig2,ph=ph,gamma=gamma,epsilon=epsilon))
            bhi["Nhok"] <- testNh(calculn$Nh,Ai,Bi,Ci,minNh)
            bhi[paste("nhok",length(Ci),sep="")] <- testnh(calculn$nh,Bi,Ci)
            bhi[paste("opti",length(Ci),sep="")] <- calculn$opti           
            for (idmaxstep in maxstep_val)
            {
                idmaxstill <- if ("change"!=rep) maxstill else floor(idmaxstep*100/3)
                
                for (idrep in 1:repid)
                {
                    A <- Ai; B <- Bi; C <- Ci
                    valid<-FALSE
                    run <- 0
                    rowrep <- rowrep+1
                    desc.rep[rowrep,"initbh.type"] <- idbhi
                    desc.rep[rowrep,"maxstep"] <- idmaxstep
                    desc.rep[rowrep,"maxstill"] <- idmaxstill
                    desc.rep[rowrep,paste("initb",1:(L-1),sep="")] <- bh <- bhi[1:(L-1)]
                    desc.rep[rowrep,"Nhok.initbh"] <- as.logical(bhi["Nhok"])
                    # Algorithme pour déterminer les bornes optimales
                    while(!valid)
                    {
                         run <- run+1
                         # Vérification du respect des conditions sur les nh pour les bornes initiales demandées et takeall
                         nhok <- as.logical(bhi[paste("nhok",length(C),sep="")])
                         if (is.na(nhok)) { # si nhok pas encore été évalué, l'obtenir et enregistrer le résultat dans bhi
                              calculn<-strata.internal(x=x,N=length(x),bh=bh,findn=findn,n=n,CV=CV,Ls=Ls,Nc=Nc,EYc=EYc,
                                                       alloc=list(q1=q1,q2=q2,q3=q3),takenone=length(A),
                                                       bias.penalty=bias.penalty,takeall=length(C),rh=rhC,model=model,
                                                       model.control=list(beta=beta,sig2=sig2,ph=ph,gamma=gamma,epsilon=epsilon))
                              nhok <- testnh(calculn$nh,B,C)
                              bhi[paste("nhok",length(C),sep="")] <- nhok
                              bhi[paste("opti",length(C),sep="")] <- calculn$opti
                         }
                         desc.rep[rowrep,"nhok.initbh"] <- nhok
                         opti <- bhi[paste("opti",length(C),sep="")]
                         # Si non respect des conditions -> bornes robustes
                         if (!desc.rep[rowrep,"Nhok.initbh"]||!desc.rep[rowrep,"nhok.initbh"]) {
                              bh <- bhi.robust[length(C)==rownames(bhi.robust),1:(L-1)] # ligne pour le bon takeall
                              if (any(is.na(bh))) {
                                   bh <- initbh.robust(x1=x1,N1=N1,Ls=Ls,takenone=length(A),takeall=length(C),minNh=minNh,wtx1=wtx1)
                                   bhi.robust[length(C)==rownames(bhi.robust),1:(L-1)] <- bh
                                   calculn<-strata.internal(x=x,N=length(x),bh=bh,findn=findn,n=n,CV=CV,Ls=Ls,Nc=Nc,EYc=EYc,
                                                            alloc=list(q1=q1,q2=q2,q3=q3),takenone=length(A),
                                                            bias.penalty=bias.penalty,takeall=length(C),rh=rhC,model=model,
                                                            model.control=list(beta=beta,sig2=sig2,ph=ph,gamma=gamma,epsilon=epsilon))
                                   bhi.robust[length(C)==rownames(bhi.robust),"Nhok"] <- testNh(calculn$Nh,A,B,C,minNh)
                                   bhi.robust[length(C)==rownames(bhi.robust),"nhok"] <- testnh(calculn$nh,B,C)
                                   bhi.robust[length(C)==rownames(bhi.robust),"opti"] <- calculn$opti
                              }
                              desc.rep[rowrep,"Nhok.robust"] <- as.logical(bhi.robust[length(C)==rownames(bhi.robust),"Nhok"])
                              desc.rep[rowrep,"nhok.robust"] <- as.logical(bhi.robust[length(C)==rownames(bhi.robust),"nhok"])
                              opti <- bhi.robust[length(C)==rownames(bhi.robust),"opti"]
                         }
                         desc.rep[rowrep,paste("b",1:(L-1),"i",sep="")] <- bh
                         pbh <- vector(length=L-1)
                         for (i in 1:(L-1)) {
                             dif <- x1-bh[i]
                             pbh[i] <- sum(dif<0)+1
                         }
                         iter <- 0
                         # Gestion des erreurs, on saute à la prochaine itération de la boucle si problème avec les bornes par défaut
                         if (!is.na(desc.rep[rowrep,"Nhok.robust"])&&(!desc.rep[rowrep,"Nhok.robust"]||!desc.rep[rowrep,"nhok.robust"])) {
                              valid <- TRUE
                              next                                                                                                        }
                         # Boucle principale : optimisation des bornes
                         desciter<-rep(0,(maxiter+1)*(L+2)) # pour la sortie du programme C
                         if ("modified"==method) {
                             resuC<-.C("KozakModif",as.double(x),as.double(x1),as.integer(wtx1),as.integer(length(x)),as.integer(N1),
                                       as.integer(findn),as.integer(n),as.double(CV),as.integer(L),as.integer(Nc),as.double(EYc),
                                       as.double(q1),as.double(q2),as.double(q3),as.integer(length(A)),as.double(bias.penalty),
                                       as.integer(length(C)),as.double(rhC),as.integer(nmodel),
                                       as.double(beta),as.double(sig2),as.double(ph),as.double(gamma),as.double(epsilon),
                                       as.integer(minNh),as.integer(idmaxstep),as.integer(maxiter),
                                       opti=as.double(opti),pbh=as.integer(pbh),desciter=as.double(desciter),
                                       iter=as.integer(0),Nh=as.integer(rep(0,L)),nhnonint=as.double(rep(0,L)),
                                       nh=as.integer(rep(0,L)),PACKAGE="stratification")
                         } else if ("original"==method) {
                             resuC<-.C("KozakOrig",as.double(x),as.double(x1),as.integer(wtx1),as.integer(length(x)),as.integer(N1),
                                       as.integer(findn),as.integer(n),as.double(CV),as.integer(L),as.integer(Nc),as.double(EYc),
                                       as.double(q1),as.double(q2),as.double(q3),as.integer(length(A)),as.double(bias.penalty),
                                       as.integer(length(C)),as.double(rhC),as.integer(nmodel),
                                       as.double(beta),as.double(sig2),as.double(ph),as.double(gamma),as.double(epsilon),
                                       as.integer(minNh),as.integer(idmaxstep),as.integer(maxiter),as.integer(idmaxstill),
                                       opti=as.double(opti),pbh=as.integer(pbh),desciter=as.double(desciter),
                                       iter=as.integer(0),Nh=as.integer(rep(0,L)),nhnonint=as.double(rep(0,L)),
                                       nh=as.integer(rep(0,L)),PACKAGE="stratification")
                         }
                         pbh<-resuC$pbh
                         opti<-resuC$opti
                         iter<-resuC$iter
                         desc.iter<-rbind(desc.iter,cbind(matrix(resuC$desciter[1:(iter*(L+2))],ncol=L+2,byrow=TRUE),run,rowrep))
                         if (identical(method,"original")) { desc.iter<-desc.iter[!apply(desc.iter[,1:(L+2)]==0,1,all),] }

                         # Vérification nh<Nh sinon -> strates fixées strates-recensement une à une
                         out.adjust<-takeallauto(nh=resuC$nh,Nh=resuC$Nh,L,A,B,C)
                                 B <- out.adjust$B ;  C <- out.adjust$C ; valid <- out.adjust$valid
                    }
        
### Je ne sais plus pourquoi j'avais mis ce bout de code. Ça ressemble à une patch pour corriger un bug.
### Je doute que ce soit encore utile, si le besoin survient, je remettrai le code actif.
###                    # Calcul de opti pour les bornes obtenues
###                    napprox <- if (1==findn) opti else n
###                    if(napprox>N & takenone>0) pbh[1]<-1 # Si napprox est plus grand que N et que la strate takenone (s'il y a lieu)
###                                                         # n'est pas vide, je vais la forcer à être vide
###                    bh<-pbh2bh(pbh,x1)
###                    calculn<-strata.internal(x=x,N=length(x),bh=bh,findn=findn,n=n,CV=CV,Ls=Ls,Nc=Nc,EYc=EYc,
###                                             alloc=list(q1=q1,q2=q2,q3=q3),takenone=length(A),
###                                             bias.penalty=bias.penalty,takeall=length(C),rh=rhC,model=model,
###                                             model.control=list(beta=beta,sig2=sig2,ph=ph,gamma=gamma,epsilon=epsilon))
###                    desc.rep[rowrep,1:(L+2)] <- c(calculn$opti,length(C),bh,iter)
#############################################################################

                    bh<-pbh2bh(pbh,x1)
                    desc.rep[rowrep,1:(L+2)] <- c(opti,length(C),bh,iter) # si j'ai fait next précédemment on a bh=bhi et iter=0,

                    # Vérification d'arrêt de la boucle
                    # utile seulement si je veux arrêter les répétitions après un certain nombre de répétitions identiques
                    # if(stoprep) { if (dim(desc.rep)[1]==10 && length(unique(desc.rep[,1]))==1 ) break }
                }
            }
        } # fin des répétitions

        repmin <- which.min(desc.rep[,1])

        # Préparation pour la sortie des résultats
        iter <- desc.rep[repmin,"niter"]
        bhi <- as.numeric(desc.rep[repmin,paste("b",1:(L-1),"i",sep="")])
        converge <- if (iter>=maxiter) FALSE else TRUE
        if (!converge) warning("the algorithm did not converge: the maximum number of iterations was reached")
        if (rep!="change") {
          if (!is.na(desc.rep[repmin,"Nhok.initbh"])&&!desc.rep[repmin,"Nhok.initbh"])
               warning("some initial sampled strata contain less than minNh units : robust initial boundaries have been used")
          if (!is.na(desc.rep[repmin,"nhok.initbh"])&&!desc.rep[repmin,"nhok.initbh"])
               warning("some initial sampled strata have non-positive nh : robust initial boundaries have been used")
          if (sum(desc.rep[,"niter"])==0)
               stop("please give an appropriate 'initbh' argument or a bigger 'minsol' argument: the robust initial boundaries give sampled strata with less than 'minNh' units and/or with non-positive nh")
        } else {
          if (sum(desc.rep[,"niter"])==0)
               stop("please give an appropriate 'initbh' argument or a bigger 'minsol' argument: the rep=='change' parameter gives no solution here because the geometric, cumrootf, default and robust initial boundaries give sampled strata with less than 'minNh' units and/or with non-positive nh")
        }
        dimnames(desc.iter)<-list(c(1:dim(desc.iter)[1]),c(paste("b",1:(L-1),sep=""),"opti","step","iter","run","rep"))
        if(1==rep) desc.iter <- desc.iter[,-(L+4),drop=FALSE]
    }
    
    # Calcul de plusieurs stat pour les bornes optimales. En fait, ces stat sont calculées à chaque
    # itération de l'algo, mais elles ne sont pas enregistrées. Ça demanderait plus de mémoire et ça alourdirait le programme inutilement.
    bh <- if(nsol<minsol) sol.detail.check[solmin,1:(L-1)] else as.numeric(desc.rep[repmin,paste("b",1:(L-1),sep="")])
    takeall.final <- if(nsol<minsol) sol.detail.check[solmin,"takeall"] else desc.rep[repmin,"takeall"]
    calculn<-strata.internal(x=x,N=length(x),bh=bh,findn=findn,n=n,CV=CV,Ls=Ls,Nc=Nc,EYc=EYc,
                             alloc=list(q1=q1,q2=q2,q3=q3),takenone=length(A),
                             bias.penalty=bias.penalty,takeall=takeall.final,rh=rhC,model=model,
                             model.control=list(beta=beta,sig2=sig2,ph=ph,gamma=gamma,epsilon=epsilon))         
    opti<-calculn$opti
    Nh<-calculn$Nh
    nh.nonint <- calculn$nhnonint
    nh <- calculn$nh
        
    # Calcul des moments anticipés
    EYh<-calculn$meanh
    VYh<-calculn$varh
    TAY<-sum(Nh[A]*EYh[A])
    TY<-sum(Nh*EYh)+Nc*EYc

    # strates pour la sortie
    bhfull<-c(min(x),bh,max(x)+1)
    stratumID <- getstratum(xgiven,N,L,certain,bhfull)

    # Préparation pour la sortie des résultats
    if (warnsol) warning("the number of possible solutions was smaller than 'minsol', therefore Kozak's algorithm was not run, instead every possible strata boundaries were tried")
    if(0==Nc) EYc <- NA;
    algo.control <- list(method=method,minNh=minNh,maxiter=maxiter,maxstep=maxstep)
    if (method=="original")  {
          if(rep!="change") {
               algo.control <- c(algo.control,list(maxstill=maxstill,rep=rep))
          } else {
               algo.control <- c(algo.control[-length(algo.control)],list(rep=rep))
          }
    }
}

#############################################################################################

## Sethi
if (algo=="Sethi")
{
    converge<-TRUE
    tol0<-min(x)*1e-8
    bhi <- initbh
   
    # Algorithme pour déterminer les bornes optimales
    valid<-FALSE
    run<-0
    while(!valid)
    {
        bhfull <- c(min(x),bhi,max(x)+1)
        iter<-0
        run<-run+1
        diff<-1
        epsilon<-0.00001
        while((iter<maxiter)&&(max(diff)>=epsilon))
        { 
            # Calcul de la taille d'echantillon n courante
            stratumID <- getstratum(xgiven,N,L,certain,bhfull)
            moments<-anticip(xgiven,stratumID,model="loglinear",beta,sig2,ph,pcertain,gamma,epsilon,A,L)
               Nh <- moments$Nh ; Nc <-  moments$Nc; phih <- moments$phih ; psih <- moments$psih ;
               EYh <- moments$EYh ; EYc <- moments$EYc ; VYh <- moments$VYh ;
               TY <- moments$TY ; TAY <- moments$TAY
            out.Alloc <- getAlloc(q1,q2,q3,B,C,Nh,Nc,EYh,VYh)
                gammah <- out.Alloc$gammah ; ah <- out.Alloc$ah ; T1 <- out.Alloc$T1
            out.opti <- optiCriteria(findn,n,CV,q1,q2,q3,bias.penalty,B,C,rh,Nh,VYh,TY,TAY,T1,ah,gammah)
               U <- out.opti$U ; U1 <- out.opti$U1 ; U2 <- out.opti$U2 ; V <- out.opti$V
               opti <- out.opti$opti
            # cat(U,V,"\n")
            desc.iter<-rbind(desc.iter,c(bhfull[-c(1,L+1)],opti,iter,run))
            # Arrêt de l'algorithme si des strates sont vides ou des variances nulles causent des divisions par zéro
            if (any(Nh==0)) {
                warning("The algorithm did not converge: division by zero caused by an empty stratum. Other intial boundaries could solve the problem.")
                nbh<-NA
            } else if (any(VYh==0)&&isTRUE(q3!=0&&q3!=1)) {
                warning("The algorithm did not converge: division by zero caused by a 0 stratum variance. Other intial boundaries could solve the problem.")
                nbh<-NA
            } else {
                # Calcul des dérivées
                dNNh<-rep(1,L); dNphih<-dNpsih<-rep(0,L)
                dENh<--ph*phih/(Nh^2)
                dEphih<-ph/Nh
                dEpsih<-rep(0,L)
                dVNh<-ph*(-exp(sig2)*psih/(Nh^2)+2*ph*(phih^2)/(Nh^3))
                dVphih<--2*ph^2*phih/(Nh^2)
                dVpsih<-ph*exp(sig2)/Nh
                dTYNh<-rep(0,L); dTYphih<-ph; dTYpsih<-rep(0,L)
                dTAYNh<-rep(0,L); dTAYphih<-c(ph[A],rep(0,length(B)+length(C))); dTAYpsih<-rep(0,L)
                dT1<-function(dN){c(rep(0,length(A)+length(B)),dN[C])[order(c(A,B,C))]}
                dT1Nh<-dT1(dNNh); dT1phih<-dT1(dNphih); dT1psih<-dT1(dNpsih)
                dU1 <- function(dN,dE,dV) {
                    dU11<-(2-2*q1)*Nh^(1-2*q1)*dN*EYh^(-2*q2)*VYh^(1-q3)
                    dU12<-Nh^(2-2*q1)*(-2*q2)*EYh^(-2*q2-1)*dE*VYh^(1-q3)
                    dU13<- if (isTRUE(q3==1)) 0 else Nh^(2-2*q1)*EYh^(-2*q2)*(1-q3)*VYh^(-q3)*dV
                    c(rep(0,length(A)),((dU11+dU12+dU13)/rh)[B],rep(0,length(C)))[order(c(A,B,C))]
                }
                dU1Nh<-dU1(dNNh,dENh,dVNh); dU1phih<-dU1(dNphih,dEphih,dVphih); dU1psih<-dU1(dNpsih,dEpsih,dVpsih); 
                dU2 <- function(dN,dE,dV) {
                    dU21<-2*q1*Nh^(2*q1-1)*dN*EYh^(2*q2)*VYh^q3
                    dU22<-Nh^(2*q1)*2*q2*EYh^(2*q2-1)*dE*VYh^q3
                    dU23<- if (isTRUE(q3==0)) 0 else Nh^(2*q1)*EYh^(2*q2)*q3*VYh^(q3-1)*dV
                    c(rep(0,length(A)),(dU21+dU22+dU23)[B],rep(0,length(C)))[order(c(A,B,C))]
                }
                dU2Nh<-dU2(dNNh,dENh,dVNh); dU2phih<-dU2(dNphih,dEphih,dVphih); dU2psih<-dU2(dNpsih,dEpsih,dVpsih); 
                dV1Nh<-CV^2*2*TY*dTYNh; dV1phih<-CV^2*2*TY*dTYphih; dV1psih<-CV^2*2*TY*dTYpsih;
                #cat(dV1Nh,dV1phih,dV1psih,"\n")
                dV2 <- function(dTAY){c(bias.penalty^2*2*TAY*dTAY[A],rep(0,length(B)+length(C)))[order(c(A,B,C))]}
                dV2Nh<-dV2(dTAYNh); dV2phih<-dV2(dTAYphih); dV2psih<-dV2(dTAYpsih);
                #cat(dV2Nh,dV2phih,dV2psih,"\n")
                dV3 <- function(dN,dV){c(rep(0,length(A)),(dN*VYh+Nh*dV)[B],rep(0,length(C)))[order(c(A,B,C))]}
                dV3Nh<-dV3(dNNh,dVNh); dV3phih<-dV3(dNphih,dVphih); dV3psih<-dV3(dNpsih,dVpsih);
                #cat(dV3Nh,dV3phih,dV3psih,"\n")
                dV4 <- function(dN,dV){c(rep(0,length(A)+length(B)),((dN*VYh+Nh*dV)*(1-1/rh))[C])[order(c(A,B,C))]}
                dV4Nh<-dV4(dNNh,dVNh); dV4phih<-dV4(dNphih,dVphih); dV4psih<-dV4(dNpsih,dVpsih);
                #cat(dV4Nh,dV4phih,dV4psih,"\n")
                dNh<-dT1Nh+((dU1Nh*U2+U1*dU2Nh)*V-U*(dV1Nh-dV2Nh+dV3Nh+dV4Nh))/(V^2)
                dphih<-dT1phih+((dU1phih*U2+U1*dU2phih)*V-U*(dV1phih-dV2phih+dV3phih+dV4phih))/(V^2)
                dpsih<-dT1psih+((dU1psih*U2+U1*dU2psih)*V-U*(dV1psih-dV2psih+dV3psih+dV4psih))/(V^2)
                #cat(dNh,dphih,dpsih,"\n")
                # Mise à jour des bornes
                a1<-dpsih[-L]-dpsih[-1]
                b1<-dphih[-L]-dphih[-1]
                #cat(b1,"\n")
                c1<-dNh[-L]-dNh[-1]
                #cat((b1**2-(4*a1*c1)),"\n")
                if (any((b1**2-(4*a1*c1))<0)) {
                    warning("The algorithm did not converge: square root of a negative number (negative discriminant). Other intial boundaries could solve the problem.")
                    nbh<-NA
                } else {
                    nbh<-ifelse(a1==0&b1==0&c1==0,bhfull[-c(1,L+1)],ifelse(a1==0,-c1/b1,(-b1+sqrt(b1**2-(4*a1*c1)))/(2*a1)))
                    nbh<-pmax(nbh,0)^(1/beta)
                }
                #cat(nbh,"\n")
            }
            if (any(is.na(nbh))) {
                diff<-rep(0,L+1)
                converge<-FALSE
            } else {
                # Une non convergence souvent observée avec les populations étudiées est une dernière strate nulle qui apparaît.
                # Je vais contraindre cette dernière strate à contenir minNh unités en corrigeant la dernière borne si elle est trop grande.
                # C'est une correction inspirée du programme original de Louis-Paul.
                bhmax <- x1[N1-sum(cumsum(rev(wtx1))<minNh)]
                if (nbh[L-1]>bhmax) nbh[L-1] <- bhmax
                # Note : Je pourrais aussi programmer des contraintes semblables sur les autres bornes mais ce serait plus compliqué,
                # sauf pour la borne 1. Je me limite tout à cette correction car c'est la seule qui se trouvait dans le programme original 
                # de Louis-Paul et parce que de toute façon on favorise l'utilisation de Kozak plutôt que Sethi.
                nbhfull<-c(min(x),nbh,max(x)+1)
                iter<-iter+1
                diff<-abs(nbhfull-bhfull)/(nbhfull+epsilon)
                if (max(diff)>=epsilon) bhfull<-nbhfull
            }
        }
        
        # Calcul de n pour la borne finale
        if (converge) {           
            stratumID <- getstratum(xgiven,N,L,certain,bhfull)
            moments<-anticip(xgiven,stratumID,model="loglinear",beta,sig2,ph,pcertain,gamma,epsilon,A,L)
               Nh <- moments$Nh ; Nc <-  moments$Nc; phih <- moments$phih ; psih <- moments$psih ;
               EYh <- moments$EYh ; EYc <- moments$EYc ; VYh <- moments$VYh ;
               TY <- moments$TY ; TAY <- moments$TAY
            out.Alloc <- getAlloc(q1,q2,q3,B,C,Nh,Nc,EYh,VYh)
               gammah <- out.Alloc$gammah ; ah <- out.Alloc$ah ; T1 <- out.Alloc$T1
            out.opti <- optiCriteria(findn,n,CV,q1,q2,q3,bias.penalty,B,C,rh,Nh,VYh,TY,TAY,T1,ah,gammah)
               opti <- out.opti$opti
            desc.iter<-rbind(desc.iter,c(bhfull[-c(1,L+1)],opti,iter,run))
        }
                    
        # Détermination de n et des nh
        out.nh <- getnh(L,A,B,C,n=opti,findn,T1,ah,Nh)
            nh.nonint <- out.nh$nh.nonint; nh <- out.nh$nh
        
        # Vérification nh<Nh sinon -> strates fixées strates-recensement une à une
        if (converge) {
            out.adjust<-takeallauto(nh,Nh,L,A,B,C)
                    B <- out.adjust$B ;  C <- out.adjust$C ; valid <- out.adjust$valid
       } else valid<-TRUE

    }
   
    # Préparation pour la sortie des résultats
    if (iter>=maxiter) {
        warning("the algorithm did not converge: the maximum number of iterations was reached")
        converge <- FALSE
    }
    dimnames(desc.iter)<-list(c(1:dim(desc.iter)[1]),c(paste("b",1:(L-1),sep=""),"n","iter","run"))
    algo.control <- list(maxiter=maxiter)
    bh<-bhfull[-c(1,L+1)]
    takeall.final <- length(C)
}
#############################################################################################

    # Calcul du CV
        # Mon programme restreint chaque strate à contenir minNh individus sauf dans la strate takenone s'il y a lieu. Il se peut alors que 
        # j'arrive à une solution finale nh=minNh pour les (L-1-takenone) dernières strates et un nh négatif dans la première strate takesome. 
        # Si un tel cas d'exception survient, je vais fixer nh=Nh pour toute strate takesome (j'aurai donc n=N).
        # if (any(calculn$nh<0)) nh[-A]<-Nh[-A]
    nh[B] <- pmax(pmin(nh[B],Nh[B]),0) # dernier ajustement, avant c'était pmax(pmin(nh[B],Nh[B]),1)
    out.RMSE <- RMSE(bias.penalty,TAY,Nh,VYh,nh,rh,B,C,TY)
          se <- out.RMSE$se ; prop <- out.RMSE$prop ; bias <- out.RMSE$bias

    # Sortie des résultats
    out<-list(Nh=Nh,nh=nh,n=sum(nh)+Nc,nh.nonint=nh.nonint,certain.info=c(Nc=Nc,meanc=EYc),opti.criteria=opti,bh=as.numeric(bh),
              meanh=ifelse(Nh==0,NA,EYh),varh=ifelse(Nh==0,NA,VYh),mean=TY/N,RMSE=se/N,RRMSE=se/TY,relativebias=bias,propbiasMSE=prop,
              stratumID=stratumID,takeall=takeall.final,call=match.call(),date=date(),
              args=list(x=xgiven,initbh=as.numeric(initbh),n=n,CV=CV,Ls=Ls,certain=certain,alloc=alloc,takenone=takenone,bias.penalty=bias.penalty,
                        takeall=takeall,rh=rh,model=model,model.control=model.control,algo=algo,algo.control=algo.control))
    if(algo=="Kozak"){
          if(nsol<minsol) { out$sol.detail <- sol.detail.check[,-c(3*L+2,3*L+3)] ; out$sol.min <- solmin ; out$nsol <- nsol
          } else {
               out$initbh <- as.numeric(bhi) ; out$iter.detail <- desc.iter ; out$niter <- iter ; out$converge <- converge
               if(rep>1) { out$rep.detail <- desc.rep ; out$rep.min <- repmin } 
          }     
    } else {
          out$initbh <- as.numeric(bhi) ; out$iter.detail <- desc.iter ; out$niter <- iter ; out$converge <- converge
    } 
    class(out)<-"strata"
    return(out)
}
