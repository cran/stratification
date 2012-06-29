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
 getnhC     
 créé le 25 avril 2012
 *******************/
 
void getnhC(double *nhnonint, int *findn, int *n, int *L, int *Nh, int *takenone, int *takeall, 
			/* sortie */ int *nh)
 { /* Arrondissement des nhnonint afin d'obtenir les nh */
 
/* Algorithme :
 Pas d'arrondissement à faire pour les strates takenone et takeall.
 Pour un CV cible : les nhnonint des strates takesome sont tous arrondis vers le haut.
 Pour un n cible (il faut que la somme des nh arrondis = n) :
 1- Ramener à 1 les nhnonint >0 et <1.
 2- Identifier les nhnonint qui reste à arrondir.
 3- Calculer le nombre de nhnonint à arrondir de chacune des façons possibles
    (nm1 : partie entière - 1, mp0 : partie entière, np1 : partie entière + 1)
	Un nm1 non nul peut survenir quand des nhnonint entre 0 et 1 sont ramenés à 1.
 4- Ordonner les parties décimales de nhnonint qui reste à arrondir en ordre croissant.
 5- Arrondissement : 
    Les nhnonint avec les nm1 plus petites parties décimales seront arrondis par : partie entière - 1,
    les nhnonint avec les np1 plus grandes parties décimales seront arrondis par : partie entière + 1,
	les autres seront arrondis par : partie entière.*/
	
/* Arguments en entrée :
nhnonint : les tailles d'échantillon par strate non entières (à arrondir)
findn : 1 si on a un CV cible, 0 si on  un n cible
n : n cible (utilisé seulement pour findn=0)
L : le nombre total de strates
Nh : tailles des strates
takenone : le nombre de strates takenone
takeall : le nombre de strates takeall*/

/* Valeur en sortie :
nh : les tailles d'échantillon par strate entières = nhnonint arrondies.*/

	int j; /* pour les boucles */

 	/* Pas d'arrondissement à faire pour les strates takenone et takeall. */
	for (j=0; j < *takenone; j++)
		nh[j] = 0;
	int snhta = 0; /* Somme des nh pour les strates takeall, utile au calcul de nhaut */
	if (*takeall>0) {
		for (j=(*L-*takeall); j < *L; j++){
			nh[j] = Nh[j];
			snhta = snhta + nh[j];
		}
	}

	/* Pour un CV cible : les nhnonint des strates takesome sont tous arrondis vers le haut.*/
	int Lts = *L - *takenone - *takeall; /* le nombre de strates takesome */
	if (*findn==1) {
		int nht[Lts]; /* parties entières des nhnonint à arrondir */
		double reste[Lts]; /* parties décimales des nhnonint à arrondir  */	
		for (j=0; j < Lts; j++) {
		/* Ce code fait la même chose que la fonction ceiling en R. */
    		nht[j] = floor(nhnonint[*takenone + j]);
    		reste[j] = nhnonint[*takenone + j] - nht[j]; 
    		if(reste[j]!=0) nh[*takenone + j] = nht[j] + 1; else nh[*takenone + j] = nht[j];
		}
	} else { /* Pour un n cible : */
	
 		/* 1- Ramener à 1 les nhnonint >0 et <1 */
		int Ltsr = Lts; /* nombre de nhnonint à arrondir parmi les strates takesome */
		for (j=0; j < Lts; j++) {
	 		if( nhnonint[*takenone + j] > 0 && nhnonint[*takenone + j] <= 1 ) {
				nh[*takenone + j] = 1; 
				Ltsr = Ltsr - 1;
			} else if ( nhnonint[*takenone + j] <= 0) { 
			/* Si un nhnonint est négatif (impossible en principe, j'ai mis ça pour m'assurer
			   de couvrir toute la droite des réels et m'assurer de bien calculer Ltsr), on le ramène à 0 */
				nh[*takenone + j] = 0; 
				Ltsr = Ltsr - 1;
			}
		}
		
		if ( Ltsr > 0 ) { /* Si on a quelques chose à arrondir seulement */

			/* 2- Identifier les nhnonint qui reste à arrondir */
			int Id[Ltsr];
			int i=0; /* Indice pour le vecteur des valeurs à arrondir uniquement */
			for (j=0; j < Lts; j++) {
				if( nhnonint[*takenone + j] > 1 ) {
					Id[i] = *takenone + j;
					i = i + 1;
				}
			}
   
			/* 3- Calculer le nombre de nhnonint à arrondir de chacune des façons possibles */
			int snht = 0; /* somme des parties entières pour les nhnonint à arrondir */
			int nht[Ltsr]; /* parties entières des nhnonint à arrondir */
			double reste[Ltsr]; /* parties décimales de nhnonint à arrondir  */	
			int index[Ltsr]; /* vecteur d'indices utile pour ordonner les parties décimales */
			for (j=0; j < Ltsr; j++) {
				nht[j] = floor(nhnonint[Id[j]]);
				reste[j] = nhnonint[Id[j]] - nht[j]; 
				snht = snht + nht[j];
				index[j] = j;			
			}
			int nhaut; /* nombre de nhnonint à arrondir vers le haut (peut être négatif) */
			nhaut = (*n - (Lts - Ltsr) - snht - snhta );
   
			int nm1, np0; /* (nm1 strates : partie entière - 1, mp0 strates : partie entière, 
							  np1 strates : partie entière + 1 (pas besoin de le calculer)) */
			if ( nhaut < 0) {
				nm1 = -nhaut;
				np0 = Ltsr - nm1;
			} else {
				nm1 = 0;
				np0 = Ltsr - nhaut;
			}
			
			/* 4- Ordonner les parties décimales de nhnonint qui reste à arrondir en ordre croissant. */
			rsort_with_index(reste,index,Ltsr); /* si égalités, fait comme rank avec ties.method="first" */
			
			/* 5- Arrondissement final */   	
			for (j=0; j < Ltsr; j++) {
				if ( j < nm1) {
					nh[Id[index[j]]] = nht[index[j]] - 1;
				} else if ( j >= nm1 && j < nm1 + np0) {
					nh[Id[index[j]]] = nht[index[j]];
				} else {
					nh[Id[index[j]]] = nht[index[j]] + 1;
				}
				/* Si un nh est inférieur à 0 (improbable mais possible), je vais le ramener à 0 */
				nh[Id[index[j]]] = fmax2(0, nh[Id[index[j]]]);
			}
		}
	}			 

 }

/* J'ai étudié les valeurs possibles de nhaut.
   *n est le ncible. Il est en fait égal à  :
	0 [strates takenone] + sum(nhnonint) [strates takesome] + snhta [strates takeall].
	On peut briser sum(nhnonint) en deux parties : une pour les nhnonint entre 0 et 1, 
	et une autre partie pour les autres nhnonint.
	On soustrait à ça 
	0 [strates takenone] + ((Lts - Ltsr)*1 + snht) [strates takesome] + snhta [strates takeall].
	On voit que les termes pour les strates takenone et takeall s'annulent.
	Écrivons le résultat de la soustraction comme la somme des deux termes suivants :
	t1 = sum(nhnonint) [strates takesome avec nhnonint entre 0 et 1] - (Lts - Ltsr)
	t2 = sum(nhnonint) [strates takesome avec nhnonint pas entre 0 et 1] - snht
	où, rappelons-le, snht est la somme des parties entères des nhnonint pas entre 0 et 1
	pour une strate takesome.
	Les valeurs max et min de ces deux termes sont :
	min(t1) tend vers - (Lts - Ltsr), max(t1) tend vers 0;
	min(t2) tend vers 0, max(t2) tend vers Ltsr.
	Ainsi, min(nhaut) = min(t1) + min(t2) : tend vers - (Lts - Ltsr)
		   max(nhaut) = max(t1) + max(t2) : tend vers Ltsr.
	C'est donc dire qu'il est impossible d'avoir à arrondir vers le haut plus de nombres
	qu'on en a, ce qui est une excellente chose.
	Cependant, on pourrait avoir des cas limite comme
	nhnonhint = 0.3333 , 0.3333, 0.3333 : serait arrondi par 1, 1, 1, mais ici ncible = 1.
	Donc les nh fournis par notre algorithme ne somment pas au ncible.
	Je pense qu'on n'a pas ici à essayer de modifier l'algo d'arrondissement. 
	C'est simplement un signe que le ncible demandé n'est pas réaliste. */
			


/*******************
 RMSEC     
 créé le 27 avril 2012
 *******************/

void RMSEC (double *biaspenalty, double *TAY, int *Nh, double *VYh, double *nh, double *rh, int *takenone, int *L,
            /* sortie */ double *RMSE)
{ /* Calcul du RMSE */

	int j; /* pour les boucles */

    /* Je peux avoir des nh nuls ou négatifs (improbable mais possible)*/
	/* Si c'est le cas, RMSE prendra une valeur manquante */
	for (j=*takenone; j < *L; j++) {
		if (nh[j] <= 0)
		  *RMSE = NA_REAL;	
	}
	if (!ISNA(*RMSE)) {
		double sum = 0; /* Somme à calculer sur les strates takesome et takeall */
		for (j=*takenone; j < *L; j++) {
			sum = sum + R_pow(Nh[j],2) * VYh[j] * (R_pow(nh[j]*rh[j],-1) - R_pow(Nh[j],-1));
		}
		*RMSE = R_pow(fmax2(0, R_pow(*biaspenalty * *TAY,2) + sum), 0.5);
	}
}
  
/* Valeurs manquantes en C, source = Wrinting R extensions, section 5.10.3 */
 

/*******************
 * strata.internal *
 * ok 11 fev 2010  *
 * modif avril 2012*
 *******************/
 
void strataC (int *findn, int *n, double *CV, int *L, int *Nc, double *EYc, 
              double *q1, double *q2, double *q3, int *takenone, double *biaspenalty, 
              int *takeall, double *rh, int *Nh, double *EYh, double *VYh, 
              /* sortie */ double *optinh, double *optinhnonint, double *nhnonint, int *nh)
{ /* partie de strata commune à tous les modèles */
	int j;
	double TY, TAY, gammah[*L], sgammah, ah[*L], T1, U, V1, V2, V3, V4, napprox, sVh;

	/* Calculs de stat utiles pour le calcul de n ou RRMSE approximatif */
	TY = TAY = sgammah = T1 = U = V3 = V4 = sVh = 0;
	for (j=0; j < *L; j++) {
		TY = TY + Nh[j] * EYh[j];
		if ( VYh[j] < 0 ) VYh[j] = 0;
	}
	TY = TY + *Nc * *EYc; /* Ajout à TY pour la strate certain */
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
    
    
    /* Calcul du critère d'optimisation (n ou RRMSE) avec nhnonint*/
    if (*findn==1) {
	    V1 = R_pow(*CV*TY,2);
	    *optinhnonint = T1 + U / (V1 - V2 + V3 + V4);  /*n*/
	    napprox = *optinhnonint;
	} else { 
		*optinhnonint = R_pow(fmax2(0, V2 + (U / (*n - T1)) - V3 - V4), 0.5) / TY;  /*RRMSE*/
		napprox = *n;
	} 
    

    /* Calcul des nh approximatif (non entiers dans les strates take-some) */
 	for (j=0; j < *takenone; j++)
		nhnonint[j] = 0;
 	for (j=*takenone; j < (*L-*takeall); j++)
		nhnonint[j] = (napprox-T1) * ah[j];
	if (*takeall>0)
		for (j=(*L-*takeall); j < *L; j++)
			nhnonint[j] = Nh[j];

			
	/* Arrondissement des nh */
    getnhC(nhnonint,findn,n,L,Nh,takenone,takeall,nh);	


    /* Calcul du critère d'optimisation (n ou RRMSE) avec nh entier*/
    if (*findn==1) {
	    *optinh = 0;
	    for (j=0; j < *L; j++)
		    *optinh = *optinh + nh[j];
	} else {
		double RMSE = 0;
		double nhd[*L];        /* La fonction RMSEC prend un nh de type double car éventuellement     */
		for (j=0; j < *L; j++) /* je l'utiliserai pour calculer le RMSE avec nh entier ou non entier. */
		    nhd[j] = nh[j];    /* Je dois donc faire une petite convertion de type ici.               */
		RMSEC(biaspenalty, &TAY, Nh, VYh, nhd, rh, takenone, L, &RMSE);
		*optinh = RMSE / TY;
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
			/* sortie */ int *Nh, double *EYh, double *VYh, double *optinh, double *optinhnonint, double *nhnonint, int *nh)
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
	strataC(findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,Nh,EYh,VYh,optinh,optinhnonint,nhnonint,nh);
}



void strataCloglinear (double *x, int *N, double *bh, int *findn, int *n, double *CV, int *L, 
			int *Nc, double *EYc, double *q1, double *q2, double *q3, int *takenone, 
			double *biaspenalty, int *takeall, double *rh, double *beta, double *sig2, double *ph, 
			/* sortie */ int *Nh, double *EYh, double *VYh, double *optinh, double *optinhnonint, double *nhnonint, int *nh)
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
	strataC(findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,Nh,EYh,VYh,optinh,optinhnonint,nhnonint,nh);
}


void strataClinear (double *x, int *N, double *bh, int *findn, int *n, double *CV, int *L, 
			int *Nc, double *EYc, double *q1, double *q2, double *q3, int *takenone, 
			double *biaspenalty, int *takeall, double *rh, double *beta, double *sig2, double *gamma, 
			/* sortie */ int *Nh, double *EYh, double *VYh, double *optinh, double *optinhnonint, double *nhnonint, int *nh)
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
	strataC(findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,Nh,EYh,VYh,optinh,optinhnonint,nhnonint,nh);
}


void strataCrandom (double *x, int *N, double *bh, int *findn, int *n, double *CV, int *L, 
			int *Nc, double *EYc, double *q1, double *q2, double *q3, 
			int *takenone, double *biaspenalty, int *takeall, double *rh, double *epsilon, 
			/* sortie */ int *Nh, double *EYh, double *VYh, double *optinh, double *optinhnonint, double *nhnonint, int *nh)
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
	strataC(findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,Nh,EYh,VYh,optinh,optinhnonint,nhnonint,nh);
}




/**********************
 * Algo Kozak modifié *
 * modif avril 2012   *
 **********************/

void KozakModif (double *x, double *x1, int *wtx1, int *N, int *N1, int *findn, int *n, double *CV, int *L, 
        int *Nc, double *EYc, double *q1, double *q2, double *q3, int *takenone, double *biaspenalty, 
        int *takeall, double *rh, int *model, double *beta, double *sig2, double *ph, double *gamma, double *epsilon, 
		int *minNh, int *maxstep, int *maxiter, int *idoptinh,
		/* sortie */ double *optinh, double *optinhnonint, int *pbh, double *desciter, int *iter, int *Nh, double *nhnonint, int *nh)
{
	int j, jbh, step, cond1, cond2, change, still, stepiter, pbhmin[*L-1], npbh[*L-1], nNh[*L], nnh[*L];
	double bh[*L-1], EYh[*L], VYh[*L], noptinh, noptinhnonint, nnhnonint[*L];
	
	*iter = stepiter = still = 0;
	pbh2bhC(pbh,x1,L,N1,bh);
	for (j=0; j < *L-1; j++) {
		pbhmin[j] = pbh[j];
		desciter[j] = bh[j];
	}
	desciter[*L-1] = *optinh;
	desciter[*L]   = *optinhnonint;
	desciter[*L+1] = stepiter;
	desciter[*L+2] = *iter;

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
							strataCnone(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,nNh,EYh,VYh,
							&noptinh,&noptinhnonint,nnhnonint,nnh);
						} else if (*model==1) {
							strataCloglinear(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,beta,sig2,ph,nNh,
							EYh,VYh,&noptinh,&noptinhnonint,nnhnonint,nnh);
						} else if(*model==2) {
							strataClinear(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,beta,sig2,gamma,nNh,
							EYh,VYh,&noptinh,&noptinhnonint,nnhnonint,nnh);
						} else if(*model==3) {
							strataCrandom(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,epsilon,nNh,EYh,VYh,
							&noptinh,&noptinhnonint,nnhnonint,nnh);
						}
						/* Vérification de la condition nh>0 pour tout h sauf strates takenone */
						cond2=0;
						for (j=*takenone; j < *L; j++) if(nnh[j]>0) cond2 = cond2 + 1;
						if (cond2==*L-*takenone) {
							/* test sur le critère à optimiser (n ou RRMSE) : a-t-il diminué? */
							if (*idoptinh) { /* Si on veux calculer le critère sur les nh entiers */
								if ( noptinh<*optinh ) { /* si noptinh est plus petit que *optinh, */
									change=1;             /* on change les bornes */
								} else if (noptinh==*optinh) { /* si noptinh est égale à *optinh, */
									/* on va comparer les critère calculé sur les nhnonint */       
									if (noptinhnonint<*optinhnonint) change=1; else change=0;
								} else {       /* si noptinh est plus grand que *optinh, */					
									change=0;  /* c'est certain qu'on ne change pas les bornes */
								} 
							} else { /* Si on veux calculer le critère sur les nhnonint (non entiers) */
								if ( noptinhnonint<*optinhnonint ) change=1; else change=0;
							}
							/* fin du test sur le critère à optimiser */
						} else {
							change=0;
						}
					} else change=0;
					/* Action posée dépendamment de l'étape précédente*/
					if(change==1) {
						for (j=0; j < *L-1; j++) pbhmin[j] = npbh[j];
						*optinh=noptinh;
						*optinhnonint=noptinhnonint;
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
			desciter[(*L+3)* *iter + j] = bh[j];
			if ( floor(desciter[(*L+3)* *iter + j]*1e8)==floor(desciter[(*L+3)*( *iter-1) + j]*1e8) ) still = still + 1;
		}
		desciter[(*L+3)* *iter + *L-1] = *optinh;
		desciter[(*L+3)* *iter + *L]   = *optinhnonint;
		desciter[(*L+3)* *iter + *L+1] = stepiter;
		desciter[(*L+3)* *iter + *L+2] = *iter;
	}
}			
					
     
/***********************
 * Algo Kozak original *
 *  ok 11 fev 2010     *
 * modif avril 2012    *
**********************/
                 
 void KozakOrig (double *x, double *x1, int *wtx1, int *N, int *N1, int *findn, int *n, double *CV, int *L, 
                 int *Nc, double *EYc, double *q1, double *q2, double *q3, int *takenone, double *biaspenalty, 
                 int *takeall, double *rh, int *model, double *beta, double *sig2, double *ph, double *gamma, 
                 double *epsilon, int *minNh, int *maxstep, int *maxiter, int *maxstill, int *idoptinh, 
                 /* sortie */ double *optinh, double *optinhnonint, int *pbh, double *desciter, int *iter, int *Nh, double *nhnonint, int *nh)
{
	int j, jbh, sign, step, cond1, cond2, change, istill, stepiter, npbh[*L-1], nNh[*L], nnh[*L];
	double bh[*L-1], EYh[*L], VYh[*L], noptinh, noptinhnonint, nnhnonint[*L];
	
	*iter = stepiter = istill = 0;
	pbh2bhC(pbh,x1,L,N1,bh);
	for (j=0; j < *L-1; j++) desciter[j] = bh[j];
	desciter[*L-1] = *optinh;
	desciter[*L]   = *optinhnonint;
	desciter[*L+1] = stepiter;
	desciter[*L+2] = *iter;

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
				strataCnone(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,nNh,EYh,VYh,
				&noptinh,&noptinhnonint,nnhnonint,nnh);
			} else if (*model==1) {
				strataCloglinear(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,beta,sig2,ph,nNh,EYh,VYh,
				&noptinh,&noptinhnonint,nnhnonint,nnh);
			} else if(*model==2) {
				strataClinear(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,beta,sig2,gamma,nNh,EYh,VYh,
				&noptinh,&noptinhnonint,nnhnonint,nnh);
			} else if(*model==3) {
				strataCrandom(x,N,bh,findn,n,CV,L,Nc,EYc,q1,q2,q3,takenone,biaspenalty,takeall,rh,epsilon,nNh,EYh,VYh,
				&noptinh,&noptinhnonint,nnhnonint,nnh);
			}
			/* Vérification de la condition nh>0 pour tout h sauf strates takenone */
			cond2=0;
			for (j=*takenone; j < *L; j++) if(nnh[j]>0) cond2 = cond2 + 1;
			if (cond2==*L-*takenone) {
				/* test sur le critère à optimiser (n ou RRMSE) : a-t-il diminué? */
				if (*idoptinh) { /* Si on veux calculer le critère sur les nh entiers */
					if ( noptinh<*optinh ) { /* si noptinh est plus petit que *optinh, */
						change=1;             /* on change les bornes */
					} else if (noptinh==*optinh) { /* si noptinh est égale à *optinh, */
					    /* on va comparer les critère calculé sur les nhnonint */       
						if (noptinhnonint<*optinhnonint) change=1; else change=0;
					} else {       /* si noptinh est plus grand que *optinh, */					
						change=0;  /* c'est certain qu'on ne change pas les bornes */
					} 
				} else { /* Si on veux calculer le critère sur les nhnonint (non entiers) */
					if ( noptinhnonint<*optinhnonint ) change=1; else change=0;
				}
				/* fin du test sur le critère à optimiser */
			} else {
				change=0;
			}
		} else change=0;
		*iter=*iter+1;
		/* Action posée dépendamment de l'étape précédente */
		if(change==1) {
			istill=0;
			pbh2bhC(npbh,x1,L,N1,bh);
			for (j=0; j < *L-1; j++) { 
				pbh[j] = npbh[j];
				desciter[(*L+3)* *iter + j] = bh[j];
			}
			*optinh=noptinh;
			*optinhnonint=noptinhnonint;
			for (j=0; j < *L; j++) {
				Nh[j]=nNh[j];
				nhnonint[j]=nnhnonint[j];
				nh[j]=nnh[j];	
			}
			desciter[(*L+3)* *iter + *L-1] = *optinh;
			desciter[(*L+3)* *iter + *L]   = *optinhnonint;
			desciter[(*L+3)* *iter + *L+1] = step;
			desciter[(*L+3)* *iter + *L+2] = *iter;
		} else istill = istill + 1 ;
	}

}

