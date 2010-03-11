`var.strata` <-
function(strata,y=NULL,rh=strata$args$rh,rh.postcorr=FALSE,model=c("none","loglinear","linear","random"),model.control=list())
{
    # Validation des arguments et initialisation de variables locales
    if (!identical(class(strata),"strata")) stop("'strata' must be an object of class 'strata'")
    x<-strata$args$x
    N<-length(x)
    takeall <- strata$takeall
    Ls <- strata$args$Ls
    takenone <- if(is.null(strata$args$takenone)) 0 else strata$args$takenone
    L <- Ls + takenone   
    bias.penalty <- if(is.null(strata$args$bias.penalty)) 1 else strata$args$bias.penalty
    
    beta <- sig2 <- ph <- pcertain <- gamma <- epsilon <- A <- B <- C <- findn <- n <- NULL
    out.check <- checkargs(y=y,rh=rh,rh.postcorr=rh.postcorr,model=model,model.control=model.control,
                           takeall=takeall,takenone=takenone,Ls=Ls,certain=strata$args$certain,testcertain=FALSE)
    for(i in 1:length(out.check)) assign(names(out.check)[i],out.check[[i]])
    
    # Ajustement à postériori pour la non-réponse, si s=demandé
    nh <- if (rh.postcorr) pmin(ceiling(strata$nh*strata$args$rh/rh),strata$Nh) else strata$nh   

    # Calcul de des espérances et variances anticipées
    if (!is.null(y)) {
          moments <- anticip(x=y,stratumID=strata$stratumID,model="none",A=A,L=L)          
    } else {
          moments <- anticip(x=x,stratumID=strata$stratumID,model,beta,sig2,ph,pcertain,gamma,epsilon,A,L)      
    }
    Nh <- moments$Nh ; Nc <- moments$Nc ; EYh <- ifelse(Nh==0,NA,moments$EYh) ; EYc <- moments$EYc ;
    VYh <- ifelse(Nh==0,NA,moments$VYh) ; TY <- moments$TY ; TAY <- moments$TAY
    
    # Pour le calcul du CV
    out.RMSE <- RMSE(bias.penalty,TAY,Nh,VYh,nh,rh,B,C,TY)
    se <- out.RMSE$se ; prop <- out.RMSE$prop ; bias <- out.RMSE$bias
    
    # Sortie des résultats
    out <- list(nh=nh,n=sum(nh)+Nc,certain.info=c(Nc=Nc,meanc=EYc),meanh=EYh,varh=VYh,mean=TY/N,RMSE=se/N,RRMSE=se/TY,
                relativebias=bias,propbiasMSE=prop,call=match.call(),date=date(),
                args=list(strata=strata,y=y,rh=rh,rh.postcorr=rh.postcorr,model=model[1],model.control=model.control))
    class(out)<-"var.strata"
    return(out)
}


`print.var.strata` <-
function(x,...)
{
    # Section des arguments fournis
    cat("Given arguments:\n") 
    cat("strata = "); print(x$call$strata)
    if (!is.null(x$args$y)) { 
        cat("y = "); print(x$call$y)
    }       
    cat("rh.postcorr =",x$args$rh.postcorr)
    if (is.null(x$args$y)) { 
        cat("\nmodel =",x$args$model)
        nparam <- if (identical(x$args$model,"loglinear")) length(x$args$model.control)-2 else length(x$args$model.control)
        if (nparam>0) {
            cat(" : ")
            for (i in 1:nparam) {
                    cat(names(x$args$model.control)[i],"=",x$args$model.control[[i]])
                    if (i<nparam) cat(" , ")
            }
        }
    }

    # Section du tableau de stratification
    takenone <- if (is.null(x$args$strata$args$takenone)) 0 else x$args$strata$args$takenone
    L<-x$args$strata$args$L+takenone
    rh <- if(takenone==0) x$args$rh else c(NA,x$args$rh)
    ph <- c(x$args$model.control$ptakenone,x$args$model.control$ph,x$args$model.control$pcertain)
    tableau<-data.frame(rh,x$args$strata$Nh,x$nh,x$nh/x$args$strata$Nh,rep("|",L),x$meanh,x$varh)
    colnames(tableau) <- c("rh","Nh","nh","fh","|","anticip. mean","anticip. var")
    if(!is.null(x$args$strata$args$certain)) {
        tableau <- rbind(tableau,c(1,x$certain.info["Nc"],x$certain.info["Nc"],1,NA,x$certain.info["meanc"],NA))
        tableau[nrow(tableau),5] <- "|"
    }
    tableau<-rbind(tableau,c(NA,sum(tableau$Nh),x$n,x$n/sum(tableau$Nh),NA,NA,NA))
    rownames(tableau) <- if(!is.null(x$args$strata$args$certain)) c(paste("stratum",1:L),"certain","Total") else c(paste("stratum",1:L),"Total")
    if (identical(x$args$model,"loglinear"))  tableau<-cbind("ph"=c(ph,NA),tableau)
    tableau[,unlist(lapply(tableau,is.numeric))]<-round(tableau[,unlist(lapply(tableau,is.numeric))],2)

    ### Correction pour afficher correctement les NA
    tableauc<-format(tableau)
    for (i in 1:(dim(tableauc)[1]-1)) tableauc[i,] <- ifelse(substr(tableauc[i,],nchar(tableauc[i,])-1,nchar(tableauc[i,]))=="NA","-",tableauc[i,])
    tableauc[dim(tableauc)[1],] <- ifelse(substr(tableauc[dim(tableauc)[1],],nchar(tableauc[dim(tableauc)[1],])-1,
                                          nchar(tableauc[dim(tableauc)[1],]))=="NA","",tableauc[dim(tableauc)[1],])
    ### Fin de la correction
    cat("\n\nStrata information:\n")
    print(tableauc,na.print="")
    cat("\nTotal sample size:",x$n,"\n")
    cat("Anticipated population mean:",x$mean,"\n")

    # Section sur les moments anticipés
    sortie <- if (is.null(x$args$strata$args$takenone)) 1 else { if(0==x$args$strata$args$takenone) 2 else 3 }
    if (sortie%in%c(1,2)) {
        cat("Anticipated CV:",x$RRMSE,"\n")
        if (2==sortie) cat("Note: CV=RRMSE (Relative Root Mean Squared Error) because takenone=0.\n")
    } else {
        est<-cbind(x$relativebias,x$propbiasMSE,x$RRMSE,x$args$strata$args$CV)
        dimnames(est) <- if(length(est)==4)  list(" ",c("relative bias","prop bias MSE","RRMSE","target RRMSE"))
                         else list(" ",c("relative bias","prop bias MSE","RRMSE"))
        cat("\nAnticipated moments of the estimator:\n")
        print.default(est, print.gap = 2, quote = FALSE, right=TRUE)
    }

}
