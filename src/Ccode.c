#include <R.h>
#include <Rmath.h>

/* Exemples de commandes pour afficher dans la console R */
   /*for (j=0; j < *L-1; j++) Rprintf("%d  ",pbh[j]); Rprintf("%d  ",isup ); Rprintf("\n");*/
   /*Rprintf("%f  ",TY ); Rprintf("\n");*/


	
/*******************
 *      getNh      *
 * ok 11 fev 2010  *
 *******************/
 
void getNhC (int *pbh, int *L, int *wtx1, int *N1, 
			 /* sortie */ int *Nh)
{ /* C version of getNh, which calculates Nh, the number of units in each stratum, 
     using the stratum boundaries expressed in terms of data rank (pbh). */
	int j,i,inf,sup;

	for (j=0; j<*L; j++) Nh[j]=0;

	if (pbh[0]>=2) {
		inf=0;
		sup=imin2(pbh[0]-1,*N1);
		for (i=inf; i < sup; i++)	Nh[0] = Nh[0] + wtx1[i];
	}
	
	if (*L>2) {
		for (j=1; j < *L-1; j++){
			if ( (pbh[j]>=2) && (pbh[j-1]<pbh[j]) && (pbh[j-1]<=*N1) ){
				inf=imax2(pbh[j-1],1)-1;
				sup=imin2(pbh[j]-1,*N1);
				for (i=inf; i < sup; i++)	Nh[j] = Nh[j] + wtx1[i];
			} else if ( (pbh[j-1]>=2) && (pbh[j-1]>pbh[j]) && (pbh[j]<=*N1) ) {
				inf=imax2(pbh[j],1)-1;
				sup=imin2(pbh[j-1]-1,*N1);
				for (i=inf; i < sup; i++)	Nh[j] = Nh[j] - wtx1[i];
			}
		}	
	}
	
	if (pbh[*L-2]<=*N1) {
		inf=imax2(pbh[*L-2],1)-1;
		sup=*N1;
		for (i=inf; i < sup; i++)	Nh[*L-1] = Nh[*L-1] + wtx1[i];
	}
	
}

/*******************
 *     pbh2bhC     *
 * ok 11 fev 2010  *
 *******************/
 
void pbh2bhC(int *pbh, double *x1, int *L, int *N1, double *bh)
{ /* makes the conversion from stratum boundaries expressed in terms of data rank (pbh), 
     to stratum boundaries expressed on the scale of the data (bh). */
    int j;
	for (j=0; j < *L-1; j++) {
		if(pbh[j]<=1) {
			bh[j] = x1[0];
		} else if (pbh[j]>=*N1) {
			bh[j] = (x1[*N1-1] + x1[*N1-2])/2;
		} else 
			bh[j] = (x1[pbh[j]-1] + x1[pbh[j]-2])/2;
	}
} 

/*******************
 * strata.internal *
 * ok 11 fev 2010  *
 *******************/
 
void strataC (int *findn, int *n, double *CV, int *L, int *Nc, double *EYc, 
              double *q1, double *q2, double *q3, int *takenone, double *biaspenalty, 
              int *takeall, double *rh, int *Nh, double *EYh, double *VYh, 
              /* sortie */ double *opti, double *nhnonint, int *nh)
{ /* partie de strata commune à tous les modèles */
	int j, ntrunc, nht[*L], index[*L], nbas;
	double TY, TAY, gammah[*L], sgammah, ah[*L], T1, U, V1, V2, V3, V4, napprox, reste[*L], sVh;

	/* Calculs de stat utiles pour le calcul de n ou RRMSE approximatif */
	TY = TAY = sgammah = T1 = U = V3 = V4 = sVh = 0;
	for (j=0; j < *L; j++) {
		TY = TY + Nh[j] * EYh[j];
		if ( VYh[j] < 0 ) VYh[j] = 0;
	}
	TY = TY + *Nc * *EYc;
/*	for (j=*takenone; j < (*L-*takeall); j++) sVh = sVh + VYh[j];
	if (sVh==0) *q3 = 0; */ /* Je ne sais pas encore si c'est pertinent de faire ça */
	for (j=0; j < *L; j++) {
		gammah[j] = R_pow(R_pow(Nh[j],2),*q1) * R_pow(R_pow(EYh[j],2),*q2) * R_pow(VYh[j],*q3);
	}    
	if (*takenone>0)
		for (j=0; j < *takenone; j++) TAY = TAY + Nh[j]*EYh[j];	    
	for (j=*takenone; j < (*L-*takeall); j++) sgammah = sgammah + gammah[j];
	for (j=*takenone; j < (*L-*takeall); j++) {
		ah[j] = gammah[j] / sgammah;
		if ( ah[j] > 0 ) U = U + R_pow(Nh[j],2) * VYh[j] / (ah[j] * rh[j]);
		V3 = V3 + Nh[j] * VYh[j];
	}
	if (*takeall>0) {
		for (j=(*L-*takeall); j < *L; j++) {
			T1 = T1 + Nh[j];
			V4 = V4 + Nh[j] * VYh[j] * ( 1 - 1 / rh[j] );
		}
	}
	T1 = T1 + *Nc;			
    V2 = R_pow(*biaspenalty*TAY,2);
    
    
    /* Calcul de n ou RRMSE approximatif*/
    if (*findn==1) {
	    V1 = R_pow(*CV*TY,2);
	    *opti = T1 + U / (V1 - V2 + V3 + V4);  /*n*/
	    napprox = *opti;
	} else { 
		*opti = R_pow(V2 + (U / (*n - T1)) - V3 - V4 , 0.5) / TY;  /*RRMSE*/
		napprox = *n;
	} 
    

    /* Calcul des nh approximatif (non entiers dans les strates take-some) */
 	for (j=0; j < *takenone; j++)
		nhnonint[j] = nh[j] = 0;
 	for (j=*takenone; j < (*L-*takeall); j++)
		nhnonint[j] = (napprox-T1) * ah[j];
	if (*takeall>0)
		for (j=(*L-*takeall); j < *L; j++)
			nhnonint[j] = nh[j] = Nh[j];
			
	/* Arrondissement des nh */  
	if (*findn==1) {
		for (j=*takenone; j < (*L-*takeall); j++) {
    		nht[j] = floor(nhnonint[j]);
    		reste[j] = nhnonint[j] - nht[j]; 
    		if(reste[j]!=0) nh[j] = nht[j] + 1; else nh[j] = nht[j];
		}
	} else { /* Je travaille sur tes les nhnonint, même ceux des strate take-none et take-all
	            qui sont déjà entiers, car c'est plus simple à programmer. */
		ntrunc = 0;
		for (j=0; j < *L; j++) {
    		nht[j] = floor(nhnonint[j]);
    		reste[j] = nhnonint[j] - nht[j]; 
			ntrunc = ntrunc + nht[j];
			index[j] = j;
		}
    	nbas = *L - (*n - ntrunc - *Nc);
    	/* Si des nh sont < 1 mais supérieur à 0, je les ramène à 1 en premier lieu */
 		for (j=0; j < *L; j++) {
	 		if( nhnonint[j]>0 && nhnonint[j]<1 ) {  
    			nht[j] = 1;
    			reste[j] = 0; 
				nbas = nbas + 1;
			}
		}
		/* Arrondissement final */   	
    	rsort_with_index(reste,index,*L); /* si égalités, fait comme rank avec ties.method="first" */
		for (j=0; j < *L; j++) {
			if ( j < nbas) {
				nh[index[j]] = nht[index[j]];
			} else {
				nh[index[j]] = nht[index[j]] + 1;
			}
		} 	
	}			 
}


void fullbh(double *x, int *N, double *bh, int *L, 
            /* sortie */ double *bhfull)
{
	int j;
	bhfull[0] = x[0];
	for (j=1; j < *L; j++) bhfull[j] = bh[j-1];
	bhfull[*L] = x[*N-1]+1;
}



void strataCnone (double *x, int *N, double *bh, int *findn, int *n, double *CV, int *L, 
			int *Nc, double *EYc, double *q1, double *q2, double *q3, 
			int *takenone, double *biaspenalty, int *takeall, double *rh, 
			/* sortie */ int *Nh, double *EYh, double *VYh, double *opti, double *nhnonint, int *nh)
{ /* strata.internal for the no model */
	int i, j;
	double bhfull[*L+1], EY2h[*L];
	
	fullbh(x,N,bh,L,bhfull);

	/* Calcul de Nh, EYh et VYh */
	for (j=0; j < *L; j++) Nh[j] = EYh[j] = EY2h[j] = 0;
	for (i=0; i < *N; i++) {
		for (j=0; j < *L; j++) {
			if ( (x[i]>=bhfull[j]) && (x[i]<bhfull[j+1]) ) {
				Nh[j] = Nh[j] + 1;
				EYh[j] = EYh[j] + x[i];
				EY2h[j] = EY2h[j] + R_pow(x[i],2);
			}
		}
	}
	for (j=0; j < *L; j++) {
		if ( Nh[j]==0 ) EYh[j] = 0; else EYh[j] = EYh[j] / Nh[j];
		if ( Nh[j]==0 ) EY2h[j] = 0; else EY2h[j] = EY2h[j] / Nh[j];
		if ( Nh[j]==0 ) VYh[j] = 0; else VYh[j] = EY2h[j] - R_pow(EYh[j],2); 
	}
	strataC(findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,Nh,EYh,VYh,opti,nhnonint,nh);
}



void strataCloglinear (double *x, int *N, double *bh, int *findn, int *n, double *CV, int *L, 
			int *Nc, double *EYc, double *q1, double *q2, double *q3, int *takenone, 
			double *biaspenalty, int *takeall, double *rh, double *beta, double *sig2, double *ph, 
			/* sortie */ int *Nh, double *EYh, double *VYh, double *opti, double *nhnonint, int *nh)
{ /* strata.internal for the loglinear model */
	int i, j;
	double bhfull[*L+1], phih[*L], psih[*L];
	
	fullbh(x,N,bh,L,bhfull);

	/* Calcul de Nh, EYh et VYh */
	for (j=0; j < *L; j++) Nh[j] = phih[j] = psih[j] = 0;
	for (i=0; i < *N; i++) {
		for (j=0; j < *L; j++) {
			if ( (x[i]>=bhfull[j]) && (x[i]<bhfull[j+1]) ) {
				Nh[j] = Nh[j] + 1;
				phih[j] = phih[j] + R_pow(x[i],*beta);
				psih[j] = psih[j] + R_pow(x[i],2* *beta);
			}
		}
	}
	for (j=0; j < *L; j++) {
		if ( Nh[j]==0 ) EYh[j] = 0; else EYh[j] = ph[j] * phih[j] / Nh[j]; 
		if ( Nh[j]==0 ) VYh[j] = 0; else VYh[j] = ph[j] * ( (exp(*sig2) * psih[j] / Nh[j]) - (ph[j] * R_pow(phih[j] / Nh[j], 2)) ); 
	}
	strataC(findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,Nh,EYh,VYh,opti,nhnonint,nh);
}


void strataClinear (double *x, int *N, double *bh, int *findn, int *n, double *CV, int *L, 
			int *Nc, double *EYc, double *q1, double *q2, double *q3, int *takenone, 
			double *biaspenalty, int *takeall, double *rh, double *beta, double *sig2, double *gamma, 
			/* sortie */ int *Nh, double *EYh, double *VYh, double *opti, double *nhnonint, int *nh)
{ /* strata.internal for the linear model */
	int i, j;
	double bhfull[*L+1], EXh[*L], EX2h[*L], EXgammah[*L], VXh[*L];
	
	fullbh(x,N,bh,L,bhfull);

	/* Calcul de Nh, EYh et VYh */
	for (j=0; j < *L; j++) Nh[j] = EXh[j] = EX2h[j] = EXgammah[j] = 0;
	for (i=0; i < *N; i++) {
		for (j=0; j < *L; j++) {
			if ( (x[i]>=bhfull[j]) && (x[i]<bhfull[j+1]) ) {
				Nh[j] = Nh[j] + 1;
				EXh[j] = EXh[j] + x[i];
				EX2h[j] = EX2h[j] + R_pow(x[i],2);
				EXgammah[j] = EXgammah[j] + R_pow(x[i],*gamma);
			}
		}
	}
	for (j=0; j < *L; j++) {
		if ( Nh[j]==0 ) EXh[j] = 0; else EXh[j] = EXh[j] / Nh[j];
		if ( Nh[j]==0 ) EX2h[j] = 0; else EX2h[j] = EX2h[j] / Nh[j];
		if ( Nh[j]==0 ) EXgammah[j] = 0; else EXgammah[j] = EXgammah[j] / Nh[j];
		VXh[j] = EX2h[j] - R_pow(EXh[j],2);
		EYh[j] = *beta * EXh[j]; 
		VYh[j] = R_pow(*beta,2) * VXh[j] + *sig2 * EXgammah[j]; 
	}
	strataC(findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,Nh,EYh,VYh,opti,nhnonint,nh);
}


void strataCrandom (double *x, int *N, double *bh, int *findn, int *n, double *CV, int *L, 
			int *Nc, double *EYc, double *q1, double *q2, double *q3, 
			int *takenone, double *biaspenalty, int *takeall, double *rh, double *epsilon, 
			/* sortie */ int *Nh, double *EYh, double *VYh, double *opti, double *nhnonint, int *nh)
{ /* strata.internal for the random replacement model */
	int i, j;
	double bhfull[*L+1], EXh[*L], EX2h[*L], EX, EX2;
	
	fullbh(x,N,bh,L,bhfull);

	/* Calcul de Nh, EYh et VYh */
	EX = EX2 =  0;
	for (j=0; j < *L; j++) Nh[j] = EXh[j] = EX2h[j] = 0;
	for (i=0; i < *N; i++) {
		EX = EX + x[i] / *N;
		EX2 = EX2 + R_pow(x[i],2) / *N;
		for (j=0; j < *L; j++) {
			if ( (x[i]>=bhfull[j]) && (x[i]<bhfull[j+1]) ) {
				Nh[j] = Nh[j] + 1;
				EXh[j] = EXh[j] + x[i];
				EX2h[j] = EX2h[j] + R_pow(x[i],2);
			}
		}
	}
	for (j=0; j < *L; j++) {
		if ( Nh[j]==0 ) EXh[j] = 0; else EXh[j] = EXh[j] / Nh[j];
		if ( Nh[j]==0 ) EX2h[j] = 0; else EX2h[j] = EX2h[j] / Nh[j];
		EYh[j] = (1 - *epsilon) * EXh[j] + *epsilon * EX; 
		VYh[j] = (1 - *epsilon) * EX2h[j] + *epsilon * EX2 - R_pow(EYh[j],2); 
	}
	strataC(findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,Nh,EYh,VYh,opti,nhnonint,nh);
}




/**********************
 * Algo Kozak modifié *
 **********************/

void KozakModif (double *x, double *x1, int *wtx1, int *N, int *N1, int *findn, int *n, double *CV, int *L, 
        int *Nc, double *EYc, double *q1, double *q2, double *q3, int *takenone, double *biaspenalty, 
        int *takeall, double *rh, int *model, double *beta, double *sig2, double *ph, double *gamma, double *epsilon, 
		int *minNh, int *maxstep, int *maxiter, 
		/* sortie */ double *opti, int *pbh, double *desciter, int *iter, int *Nh, double *nhnonint, int *nh)
{
	int j, jbh, step, cond1, cond2, change, still, stepiter, pbhmin[*L-1], npbh[*L-1], nNh[*L], nnh[*L];
	double bh[*L-1], EYh[*L], VYh[*L], nopti, nnhnonint[*L];
	
	*iter = stepiter = still = 0;
	pbh2bhC(pbh,x1,L,N1,bh);
	for (j=0; j < *L-1; j++) {
		pbhmin[j] = pbh[j];
		desciter[j] = bh[j];
	}
	desciter[*L-1] = *opti;
	desciter[*L] = stepiter;
	desciter[*L+1] = *iter;

	while ( (*iter<*maxiter) && (still!=*L-1) ) {
		for (jbh=0; jbh < *L-1; jbh++) {
			for (step=-(*maxstep); step<=*maxstep; step++) {
				if (step!=0) {
					/* Modification des bornes */
					for (j=0; j < *L-1; j++) npbh[j] = pbh[j];
					npbh[jbh] = pbh[jbh] + step;
					/* Vérification du respect de la condition nNh>=minNh pour tout h sauf strates takenone*/
					getNhC(npbh,L,wtx1,N1,nNh);
					cond1=0;
					for (j=0; j<*takenone; j++) if( (nNh[j]>=0) && (npbh[j]>=1) ) cond1 = cond1 + 1; 
					for (j=*takenone; j < *L; j++) if(nNh[j]>=*minNh) cond1 = cond1 + 1;
					/* Si au moins un Nh ne vérifie pas la condition, on ne fait pas les étapes suivantes*/
					if (cond1==*L) {
						pbh2bhC(npbh,x1,L,N1,bh);
						if (*model==0) {
							strataCnone(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,nNh,EYh,VYh,&nopti,nnhnonint,nnh);
						} else if (*model==1) {
							strataCloglinear(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,beta,sig2,ph,nNh,EYh,VYh,&nopti,nnhnonint,nnh);
						} else if(*model==2) {
							strataClinear(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,beta,sig2,gamma,nNh,EYh,VYh,&nopti,nnhnonint,nnh);
						} else if(*model==3) {
							strataCrandom(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,epsilon,nNh,EYh,VYh,&nopti,nnhnonint,nnh);
						}
						/* Vérification de la condition nh>0 pour tout h sauf strates takenone */
						cond2=0;
						for (j=*takenone; j < *L; j++) if(nnh[j]>0) cond2 = cond2 + 1;
						/* test sur n ou RRMSE : a-t-il diminué? */
						/*************/ if ( (nopti<*opti) && (cond2==*L-*takenone) ) change=1; else change=0;
					} else change=0;
					/* Action posée dépendamment de l'étape précédente*/
					if(change==1) {
						for (j=0; j < *L-1; j++) pbhmin[j] = npbh[j];
						*opti=nopti;
						for (j=0; j < *L; j++) {
							Nh[j]=nNh[j];
							nhnonint[j]=nnhnonint[j];
							nh[j]=nnh[j];	
						}
						stepiter=step;
					}
				}
			}
		}
		*iter=*iter+1;
		still=0;
		pbh2bhC(pbhmin,x1,L,N1,bh);
		for (j=0; j < *L-1; j++) {
			pbh[j] = pbhmin[j];
			desciter[(*L+2)* *iter + j] = bh[j];
			if ( floor(desciter[(*L+2)* *iter + j]*1e8)==floor(desciter[(*L+2)*( *iter-1) + j]*1e8) ) still = still + 1;
		}
		desciter[(*L+2)* *iter + *L-1] = *opti;
		desciter[(*L+2)* *iter + *L] = stepiter;
		desciter[(*L+2)* *iter + *L+1] = *iter;
	}
}			
					
     
/***********************
 * Algo Kozak original *
 *  ok 11 fev 2010     *
 **********************/
                 
 void KozakOrig (double *x, double *x1, int *wtx1, int *N, int *N1, int *findn, int *n, double *CV, int *L, 
                 int *Nc, double *EYc, double *q1, double *q2, double *q3, int *takenone, double *biaspenalty, 
                 int *takeall, double *rh, int *model, double *beta, double *sig2, double *ph, double *gamma, 
                 double *epsilon, int *minNh, int *maxstep, int *maxiter, int *maxstill, 
                 /* sortie */ double *opti, int *pbh, double *desciter, int *iter, int *Nh, double *nhnonint, int *nh)
{
	int j, jbh, sign, step, cond1, cond2, change, istill, stepiter, npbh[*L-1], nNh[*L], nnh[*L];
	double bh[*L-1], EYh[*L], VYh[*L], nopti, nnhnonint[*L];
	
	*iter = stepiter = istill = 0;
	pbh2bhC(pbh,x1,L,N1,bh);
	for (j=0; j < *L-1; j++) desciter[j] = bh[j];
	desciter[*L-1] = *opti;
	desciter[*L] = stepiter;
	desciter[*L+1] = *iter;

		while ( (*iter<*maxiter) && (istill<*maxstill) ) {	
		/* Choix aléatoire de la borne à remplacer et du pas fait*/
		GetRNGstate(); 
		jbh = floor(unif_rand()* (*L-1)); /* 0, 1, ..., L-2 */
			if (jbh==(*L-1)) jbh=0; /* si jbh vaut L-1, ce qui a une probabilité presque nulle */
		step =  floor(unif_rand()* *maxstep) + 1; /* 1, 2, ..., *maxstep */
			if (step==(*maxstep+1)) step=1; /* si step vaut *maxstep+1, ce qui a une probabilité presque nulle */
		sign = floor(unif_rand()*2); /* 0 ou 1 */
		if (sign==0) step=-step;
		PutRNGstate();
		
		/* Modification des bornes */
		for (j=0; j < *L-1; j++) npbh[j] = pbh[j];
		npbh[jbh] = pbh[jbh] + step;
		/* Vérification du respect de la condition nNh>=minNh pour tout h sauf strates takenone*/
		getNhC(npbh,L,wtx1,N1,nNh);
		cond1=0;
		for (j=0; j<*takenone; j++) if( (nNh[j]>=0) && (npbh[j]>=1) ) cond1 = cond1 + 1; 
		for (j=*takenone; j < *L; j++) if(nNh[j]>=*minNh) cond1 = cond1 + 1;
		/* Si au moins un Nh ne vérifie pas la condition, on ne fait pas les étapes suivantes*/
		if (cond1==*L) {
			pbh2bhC(npbh,x1,L,N1,bh);
			if (*model==0) {
				strataCnone(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,nNh,EYh,VYh,&nopti,nnhnonint,nnh);
			} else if (*model==1) {
				strataCloglinear(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,beta,sig2,ph,nNh,EYh,VYh,&nopti,nnhnonint,nnh);
			} else if(*model==2) {
				strataClinear(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,beta,sig2,gamma,nNh,EYh,VYh,&nopti,nnhnonint,nnh);
			} else if(*model==3) {
				strataCrandom(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,epsilon,nNh,EYh,VYh,&nopti,nnhnonint,nnh);
			}
			/* Vérification de la condition nh>0 pour tout h sauf strates takenone */
			cond2=0;
			for (j=*takenone; j < *L; j++) if(nnh[j]>0) cond2 = cond2 + 1;
			/* test sur n ou RRMSE : a-t-il diminué? */
			/*************/ if ( (nopti<*opti) && (cond2==*L-*takenone) ) change=1; else change=0;
		} else change=0;
		*iter=*iter+1;
		/* Action posée dépendamment de l'étape précédente */
		if(change==1) {
			istill=0;
			pbh2bhC(npbh,x1,L,N1,bh);
			for (j=0; j < *L-1; j++) { 
				pbh[j] = npbh[j];
				desciter[(*L+2)* *iter + j] = bh[j];
			}
			*opti=nopti;
			for (j=0; j < *L; j++) {
				Nh[j]=nNh[j];
				nhnonint[j]=nnhnonint[j];
				nh[j]=nnh[j];	
			}
			desciter[(*L+2)* *iter + *L-1] = *opti;
			desciter[(*L+2)* *iter + *L] = step;
			desciter[(*L+2)* *iter + *L+1] = *iter;
		} else istill = istill + 1 ;
	}

}


/*******************
 *    checknhC     *
 * inutilisée à partir du 15 février car on adopte une nouvelle approche 
 * pour les bornes initiales (voir fichier internals.R) 
 *******************/


			 
void checknhC (double *x, double *x1, int *N, int *findn, int *n, double *CV, int *L, 
			int *Nc, double *EYc, double *q1, double *q2, double *q3, int *takenone, double *biaspenalty, 
			int *takeall, double *rh, int *model, double *beta, double *sig2, double *ph, 
            double *gamma, double *epsilon, int *minNh, int *Nh, 
            /* sortie */ int *pbh)
{ /* modification des bornes si les nh ne sont pas acceptés */
/* If the initial boundaries are associated with non positive
nh, we correct the boundaries by decreasing the size of the takenone or take-all strata. At each step of the
correction algorithm, we modify the boundaries in order to remove one unit from the takenone stratum (if there is one) or the last take-all
stratum and add a unit in the smallest take-some stratum. The process is stopped when no
more nh are non positive or the takesome strata contains 0 units and every take-all strata contains minNh units.*/

	int i, nhneg, Nhok, NhBmin, imin, iinf, isup, nh[*L];
	double bh[*L-1], EYh[*L], VYh[*L], opti, nhnonint[*L];
	
	nhneg = Nhok = 1;
	while ( (nhneg>=1) && (Nhok==1) )
	{
		/* Afin de trouver la position du Nh minimum pour les strates takesome */
		NhBmin=Nh[*takenone];
		imin=*takenone;
		for (i=*takenone+1; i<*L-*takeall; i++) if ( Nh[i]<NhBmin ) { NhBmin=Nh[i]; imin=i; }
		/* Les bornes à modifier le sont */
		if ( *takenone>0 && Nh[*takenone-1]>0 ) {
			/* S'il y a une strate takesome, on commence par la réduire */
			isup=imin-1;
			iinf=0;
			if (*takenone>1) for (i=0; i<*takenone-1; i++) if(Nh[i]<=0) iinf=iinf+1;
			for (i=iinf; i<=isup; i++) pbh[i]=pbh[i]-1;
		} else if ( *takeall>0 && Nh[*L-*takeall]>*minNh ) {
			/* Ensuite, on réduit les strates takeall. */ 
			iinf=imin;
			isup=*L-1-1; 
			if (*takeall>1) for (i=*L-1; i>*L-*takeall; i--) if(Nh[i]<=*minNh) isup=isup-1;
			for (i=iinf; i<=isup; i++) pbh[i]=pbh[i]+1;
		} else Nhok=0; /* Afin de ne pas avoir une boucle infinie si les bornes ne peuvent être modifiées */ 
		/* Le plan d'échantillonnage est obtenu avec les nouvelles bornes */
		for (i=0; i < *L-1; i++) if(pbh[i]<=1) bh[i] = x1[0]; else bh[i] = (x1[pbh[i]-1] + x1[pbh[i]-2])/2;
		if (*model==0) {
			strataCnone(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,Nh,EYh,VYh,&opti,nhnonint,nh);
		} else if (*model==1) {
			strataCloglinear(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,beta,sig2,ph,Nh,EYh,VYh,&opti,nhnonint,nh);
		} else if(*model==2) {
			strataClinear(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,beta,sig2,gamma,Nh,EYh,VYh,&opti,nhnonint,nh);
		} else if(*model==3) {
			strataCrandom(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,epsilon,Nh,EYh,VYh,&opti,nhnonint,nh);
		}
		/* Les conditions sont évaluées */
		nhneg=0;
		for (i=*takenone; i<*L; i++) if (nh[i]<=0) nhneg=nhneg+1;
		if (Nh[*L-*takeall]<=*minNh) Nhok=0;
		/*for (i=0; i < *L-1; i++) Rprintf("%d  ",pbh[i]); Rprintf("imin: %d  ",imin ); Rprintf("isup: %d  ",isup ); Rprintf("\n");*/
	}	
	
}



