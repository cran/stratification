#include <R.h>
#include <Rmath.h>

/* Note : le code contient des Rprintf en commentaires. Ceux-ci servent a debugger le code au besoin

Exemples de commandes pour afficher dans la console R
   for (j=0; j < *L-1; j++) Rprintf("%d  ",pbh[j]); Rprintf("%d  ",isup ); Rprintf("\n");
   Rprintf("%f  ",TY ); Rprintf("\n");
   
%d = entier (integer), %f = reel (double)  */


/* Definition des parametres communs a plusieurs fonctions 

Un nom de parametre represente toujours la meme chose.
Deux suffixes reviennent a quelques reprises :
...noc : no certainty stratum (n'inclut pas la strate certain)
...nonint : non integer (non entier)

## arguments en lien avec les donnees et crees en R :

@param xnoc observations de la variable de stratification donnees en entree, excluant les  
            observations pour la certainty stratum
@param Nnoc nombre total d'observations apres avoir retire les observations pour la certainty stratum
            (longueur de xnoc et de stratumIDnoc) 
@param x1noc vecteur des valeurs uniques, ordonnees de la plus petite a la plus grande, 
          dans les observations de la variable de stratification donnees en entree, excluant les  
            observations pour la certainty stratum (xnoc)
@param N1noc nombre de valeurs uniques dans les observations xnoc (longueur de x1noc)
@param Nc le nombre d'observations dans la strate certain
@param EYc l'esperance anticipee de Y dans la strate certain
@param EX moyenne de toutes les observations, peu importe leur strate, incluant les obs de la strate certain
          (utile au modele random)
@param EX2 moyenne de toutes les observations au carre, peu importe leur strate, incluant les obs de la strate certain
           (utile au modele random)

## arguments qui sont en fait des arguments donnees en entree a une
   fonction R publique (et possiblement formate) :
   
@param bhfull les bornes des strates, incluant b0 et bL (vecteur de longueur L + 1)
@param L nombre total de strates (ne compte pas la strate certain) 
         -> L doit etre >=2 pour assurer le bon fonctionnement des fonctions
@param takenone le nombre de strate take-none : 0 ou 1
@param takeall le nombre de strate take-all : 0, 1, ou plus rarement un nombre entier > 1
@param q1 le premier exposant definissant l'allocation
@param q2 le deuxieme exposant definissant l'allocation
@param q3 le troisieme exposant definissant l'allocation
@param findn indicatrice : 1 si on a un CV cible, 0 si on  un n cible
@param n le n cible donne en entree
@param CV RRMSE cible donne en entree a une fonction externe
@param rhL taux de reponse : vecteur de longueur L dont la premiere valeur est 1 s'il y a une strate takenone,
           donc qui n'est pas necessairement identique au rh donne en entree a une fonction externe.
           (La valeur 1 est sans importance car elle n'est pas utilisee. C'est la position des autres valeurs 
            qui est importante.) 
@param biaspenalty argument donne en entree a une fonction externe
@param nmodelID un identifiant du modele 0 = none, 1 = loglinear, 2 = linear, 3 =  random
@param beta
@param sig2
@param ps 
@param ph vecteur de longueur L pour les taux de mortalite dans les L strates 
@param gamma
@param epsilon
@param minNh parametre de l'algo donne en entree a la fonction externe : 
             valeur minimale de Nh acceptee pour les strates echantillonnees (take-some et take-all)

## arguments pour des fonctions C internes

@param xs le vecteur double des observations dans une state specifique (notee strate s)
@param Ns le nombre d'observations dans la strate s (longueur de xs)
@param nhcalcul les nh choisis pour faire le calcul, parfois des reels, d'autres fois des entiers

## arguments qui sont obtenus d'une fonction C

@param stratumIDnoc tel que cree par get_stratumIDnoc_C
@param Nh vecteur de longueur L des tailles de populations obtenu de get_Nh_C
@param VYh la variance anticipee de Y dans chaque strate obtenue de get_momentY_C
@param ah vecteur de longueur L : ah pour les strates take-some, 0 pour les autres strates, obtenu de get_ah_C
@param TCN la somme des tailles de populations Nh pour les strates take-all (C), obtenue de get_nnonint_C
@param nhnonint un vecteur de L nombres reels : les tailles d'echantillon dans les L strates, 
                issues de l'application de la regle d'allocation, obtenu de get_nhnonint_C                
@param TAY la somme anticipee de Y dans les strates take-none, obtenue de get_momentY_C

*/

/* Definition des valeurs retournees (ou plutot modifiees) communes a plusieurs fonctions

@return stratumIDnoc un vecteur de longueur Nnoc contenant les chiffres 1 a L, identifiant la strate
                     a laquelle chaque observation appartient
@return Nh vecteur de longueur L des tailles de populations
@return EYh Un vecteur de longueur L : l'esperance anticipee de Y dans chaque strate
@return VYh Un vecteur de longueur L : la variance anticipee de Y dans chaque strate
@return TY Un nombre de type double : la somme gloable anticipee de Y (incluant la certainty stratum)
@return nhnonint un vecteur de L nombre reel : les tailles d'echantillon dans les L strates, 
                        issues de l'application de la regle d'allocation
@return nh les tailles d'echantillon entieres pour les L strates = nhnonint arrondies.
*/

/* ******************************************************************************************** */

/* Determination des strates en ne considerant pas la strate certain

Pour chaque observation, cette fonction identifie sa strate.

Aucun parametre n'est unique a cette fonction, 
voir la description generale des parametres au debut de ce fichier.

Aucune valeur retournee (ou modifiee) n'est unique a cette fonction, 
voir la description generale des sorties au debut de ce fichier.
  
Creee et verifiee le 26 septembre 2012
voir wrapper R
*/
void get_stratumIDnoc_C (double *xnoc, int *Nnoc, double *bhfull, int *L,
                      /* Sortie */ int *stratumIDnoc)
{
  int i, j;
  for (i=0; i < *Nnoc; i++) {
		for (j=0; j < *L; j++) {
			if ( (xnoc[i] >= bhfull[j]) && (xnoc[i] < bhfull[j+1]) ) {
				stratumIDnoc[i] = j + 1;
			}
		}
	}
  /* La fonction retourne un resultat meme si bhfull est illogique (voir testDevel.R) */
}

/* ******************************************************************************************** */

/* Determination des tailles de population dans les strates (sans considerer la strate certain)

Aucun parametre n'est unique a cette fonction, 
voir la description generale des parametres au debut de ce fichier.

Aucune valeur retournee (ou modifiee) n'est unique a cette fonction, 
voir la description generale des sorties au debut de ce fichier.

Creee et verifiee le 26 septembre 2012
voir wrapper R
*/
void get_Nh_C (int *stratumIDnoc, int *Nnoc, int *L,
                /* Sortie */ int *Nh)
{
  int i, j;
  
  for (j=0; j < *L; j++) Nh[j] = 0;
  
  for (i=0; i < *Nnoc; i++) {
  	for (j=0; j < *L; j++) {
			if ( stratumIDnoc[i] == j + 1 ) {
				Nh[j] = Nh[j] + 1;
			}
		}
	}
}

/* ******************************************************************************************** */

/* Fonction qui extrait les observations d'une strate

@param nstatum le numero de la strate que l'on souhaite extraire
voir la description generale au debut de ce fichier pour les autres parametres.

@return xs observations pour la strate nstratum

Creee et verifiee avec valeurs par defaut le 27 septembre 2012
pas de wrapper R, fonction C seulement
*/
void extract_stratum_C (int *nstratum, double *xnoc, int *stratumIDnoc, int *Nnoc,
                         /* Sortie */ double *xs)
{
  int i, j=0;
  
  for (i=0; i < *Nnoc; i++) {
    if (stratumIDnoc[i] == *nstratum) {
      xs[j] = xnoc[i];
      j = j + 1;
    }
  }
}

/* ******************************************************************************************** */

/* Obtention de E[Y] pour une strate specifique s

Calcul l'esperance anticipee de Y pour les observations d'une seule strate, la state s.

Aucun parametre n'est unique a cette fonction, 
voir la description generale des parametres au debut de ce fichier.

@return \item{EYs}{Un nombre de type double : l'esperance anticipee de Y pour les observations de la state s}
@return \item{EXs}{Un nombre de type double : l'esperance de X pour les observations de la state s 
                   (uniquement utile a get_VYs_C)}
@return \item{phis}{Un nombre de type double : la valeur de phi pour la state s 
                    (utile a get_VYs_C et get_momentY_C)}
                    
Precondition : Ns est positif

Creee et verifiee avec valeurs par defaut le 26 septembre 2012
voir wrapper R
*/
void get_EYs_C (double *xs, int *Ns, int *nmodel, 
                  double *beta, double *sig2, double *ps, double *gamma, double *epsilon, double *EX,
                  /* Sortie */ double *EYs, double *EXs, double *phis)
{
  int i;  
  
  /* Premiere etape : sommes */
  *EXs = *phis = 0;
  for (i = 0; i < *Ns; i++){ 
    if (*nmodel == 1) {
      *phis = *phis + R_pow(xs[i], *beta);
    } else {
      *EXs = *EXs + xs[i] / *Ns;
    }
  }
  
  /* Deuxieme etape : multiplications ou divisions */
  if (*nmodel == 1) {
    *EYs = *ps * *phis / *Ns;
  } else {
    if (*nmodel == 0) *EYs = *EXs;
    if (*nmodel == 2) *EYs = *beta * *EXs;
    if (*nmodel == 3) *EYs = (1 - *epsilon) * *EXs + *epsilon * *EX;
  }
}

/* ******************************************************************************************** */

/* Obtention de Var[Y] pour une strate specifique s

Calcul la variance anticipee de Y pour les observations d'une seule strate, la state s.

@param EYs l'esperance anticipee de Y pour les observations de la state s obtenue de get_EYs_C 
@param EXs l'esperance de X pour les observations de la state s obtenue de get_EYs_C
@param phis la valeur de phi pour la state s obtenue de get_EYs_C
voir la description generale au debut de ce fichier pour les autres parametres.

@return \item{VYs}{Un nombre de type double : la variance anticipee de Y pour les observations de la state s}
@return \item{psis}{Un nombre de type double : la valeur de psi pour la state s (uniquement utile a get_momentY_C)}

Precondition : Ns est positif

Creee et verifiee avec valeurs par defaut le 27 septembre 2012
pas de wrapper R, fonction C seulement
*/
void get_VYs_C (double *xs, int *Ns, double *EYs, double *EXs, double *phis, int *nmodel, 
                  double *beta, double *sig2, double *ps, double *gamma, double *epsilon, double *EX, double *EX2,
                  /* Sortie */ double *VYs, double *psis)
{
  int i;  
  double EX2s = 0, EXgammas = 0, VXs;
  
  /* Premiere etape : sommes */
  *psis = 0;
  for (i = 0; i < *Ns; i++){ 
    if (*nmodel == 1) {
      *psis = *psis + R_pow(xs[i], 2 * *beta);
    } else {
      EX2s = EX2s + R_pow(xs[i], 2) / *Ns;
      if (*nmodel == 2) EXgammas = EXgammas + R_pow(xs[i], *gamma) / *Ns;
    }
  }
  
  /* Deuxieme etape : multiplications ou divisions */
  if (*nmodel == 1) {
    *VYs = *ps * ( (exp(*sig2) * *psis / *Ns) - (*ps * R_pow(*phis / *Ns, 2)) );
  } else {
    VXs = EX2s - R_pow(*EXs, 2);
    if (*nmodel == 0) *VYs = VXs;
    if (*nmodel == 2) *VYs = R_pow(*beta,2) * VXs + *sig2 * EXgammas;
    if (*nmodel == 3) *VYs = (1 - *epsilon) * EX2s + *epsilon * *EX2 - R_pow(*EYs, 2);
  }
  
  /* Derniere etape : ajustement cas extreme */
  if (*VYs < 0) *VYs = 0;
  /* Cet ajustement est utile dans un cas ou tous les individus d'une strate ont la meme valeur,
     donc la variance est nulle. Mais a cause de l'ecriture de notre formule pour VYs, on obtient
     plutot VYh[j]=-0. Cette valeur negative peut ensuite faire planter R dans get_gammah_C lors
     de l'execution de la fonction R_pow(VYh[j],*q3) si q3 n'est pas entier. */
}

/* ******************************************************************************************** */

/* Obtention de E[Y], Var[Y] pour toutes les strates (excluant la strate certain) + TY et TAY

Calcul des esperances et variances anticipees de Y pour les observations de chacune des L strates.
Aussi, calcul de TY et TAY.

voir la description generale au debut de ce fichier pour les autres parametres.

@return \item{phih}{Un vecteur de longueur L : la valeur de phi dans chaque strate 
                   (uniquement utile a l'algo de Sethi dans la fonction R strata.LH)}
@return \item{psih}{Un vecteur de longueur L : la valeur de psi dans chaque strate
                   (uniquement utile a l'algo de Sethi dans la fonction R strata.LH)}
@return \item{TAY}{Un nombre de type double : la somme anticipee de Y dans les strates take-none}
voir la description generale au debut de ce fichier pour les autres valeurs en sortie.

Creee et verifiee avec valeurs par defaut le 27 septembre 2012
voir wrapper R
*/
void get_momentY_C (double *xnoc, int *stratumIDnoc, int *Nnoc, int *Nh, int *L, int *Nc, double *EYc, int *takenone,
                      int *nmodel, double *beta, double *sig2, double *ph, double *gamma, double *epsilon, 
                      double *EX, double *EX2,
                      /* Sortie */ double *EYh, double *VYh, double *phih, double *psih, double *TY, double *TAY)
{
  int j;
  
  *TY = 0;
  for (j=0; j < *L; j++){
    int Ns = Nh[j];
    if (Ns > 0) {
      int nstratum = j + 1;
      double xs[Ns], EXs;
      extract_stratum_C(&nstratum, xnoc, stratumIDnoc, Nnoc, xs);
      get_EYs_C(xs, &Ns, nmodel, beta, sig2, &ph[j], gamma, epsilon, EX, &EYh[j], &EXs, &phih[j]);
      get_VYs_C(xs, &Ns, &EYh[j], &EXs, &phih[j], nmodel, beta, sig2, &ph[j], gamma, epsilon, EX, EX2, 
                    &VYh[j], &psih[j]);
      *TY = *TY + Ns * EYh[j];
    } else {
      /* valeurs non utilisees je crois, mais je ne prends pas de chance */
      *EYh = *VYh = *phih = *psih = 0;
    }
  }
  if (*Nc > 0) *TY = *TY + *Nc * *EYc; /* Ajout a TY pour la strate certain */
  
  *TAY = 0;
  if (*takenone>0){
		for (j=0; j < *takenone; j++) {
      if (Nh[j] > 0) *TAY = *TAY + Nh[j] * EYh[j];	      
		}
  }
} 

/* ******************************************************************************************** */

/* Obtention de gammah pour les L strates, dans le but de faire l'allocation

@param EYh l'esperance anticipee de Y dans chaque strate obtenue de get_momentY_C
voir la description generale au debut de ce fichier pour les autres parametres.

@return \item{gammah}{Un vecteur de longueur L : gammah pour toutes les strates}

Creee et verifiee le 27 septembre 2012
*/
void get_gammah_C (int *Nh, double *EYh, double *VYh, int *L, double *q1, double *q2, double *q3,
                      /* Sortie */ double *gammah)
{
  int j;  
  for (j=0; j < *L; j++) {
		gammah[j] = R_pow(R_pow(Nh[j],2),*q1) * R_pow(R_pow(EYh[j],2),*q2) * R_pow(VYh[j],*q3);
	}    
}

/* ******************************************************************************************** */

/* Obtention de ah pour les strates take-some uniquement, dans le but de faire l'allocation

Cette fonction calcule aussi TCN car cette valeur est utile a get_nnoc_C, qu'on utilise seulement
avec un CV cible, mais aussi a la fonction get_nhnonint_C qui est utilise peu importe le type
de cible (CV ou n).

@param gammah vecteur gammah pour les L strates, obtenu de get_gammah_C
voir la description generale au debut de ce fichier pour les autres parametres.

@return \item{ah}{Un vecteur de longueur L : ah pour les strates take-some, 0 pour les autres strates}
@return \item{TCN}{un nombre entier : la somme des tailles de populations Nh pour les strates take-all}
@return \item{U2}{Un nombre reel : U2 = sommes des gamma des strates take-some
                  (uniquement utile a l'algo de Sethi dans la fonction R strata.LH)}
}

Creee et verifiee le 27 septembre 2012
*/
void get_ah_C (double *gammah, int *Nh, int *L, int *takenone, int *takeall,
                   /* Sortie */ double *ah, int *TCN, double *U2)
{
  int j;
  double sgammah = 0;
  
  for (j=0; j < *L; j++) ah[j] = 0;
  for (j=*takenone; j < (*L-*takeall); j++) sgammah = sgammah + gammah[j];
  *U2 = sgammah;
	for (j=*takenone; j < (*L-*takeall); j++) ah[j] = gammah[j] / sgammah;

  *TCN = 0;
  if (*takeall>0)
		for (j=(*L-*takeall); j < *L; j++)
			*TCN = *TCN + Nh[j];		
}

/* ******************************************************************************************** */

/* Obtention de nnoc reel, uniquement utile dans le cas d'un CV cible

@param TY la somme gloable anticipee de Y (incluant la certainty stratum), obtenue de get_momentY_C
voir la description generale au debut de ce fichier pour les autres parametres.

@return \item{nnoc}{un nombre reel : la taille d'echantillon permettant d'atteindre le CV cible 
                       (excluant la strate certain)}
@return \item{U}{Un nombre reel : U (uniquement utile a l'algo de Sethi dans la fonction R strata.LH)}
@return \item{V}{Un nombre reel : V (uniquement utile a l'algo de Sethi dans la fonction R strata.LH)}

Creee et verifiee le 27 septembre 2012
*/
void get_nnoc_C (double *CV, int *Nh, double *ah, double *VYh, double *rhL, int *L, 
                      int *takenone, int *takeall, double *biaspenalty, double *TAY, double *TY, int *TCN, 
                      /* Sortie */ double *nnoc, double *U, double *V)
{
  int j;
  double V1, V2, V3 = 0, V4 = 0;
  
  V1 = R_pow(*CV * *TY, 2);
  V2 = R_pow(*biaspenalty * *TAY, 2);
  
  *U = 0;
  for (j=*takenone; j < (*L-*takeall); j++) {
    if (Nh[j] != 0) {
		  if (Nh[j] != 0 && VYh[j] != 0) *U = *U + R_pow(Nh[j],2) * VYh[j] / (ah[j] * rhL[j]);
		  V3 = V3 + Nh[j] * VYh[j];
    }
	}
  
	if (*takeall>0){
		for (j=(*L-*takeall); j < *L; j++){
			if (Nh[j] != 0) V4 = V4 + Nh[j] * VYh[j] * ( 1 - 1 / rhL[j] );
		}
	}
      
  *V = V1 - V2 + V3 + V4;
	if (*U == 0) *nnoc = *TCN; else *nnoc = *TCN + *U / *V;
}

/* ******************************************************************************************** */

/* Allocation = obtention des taille d'echantillon nh reel pour les L strates

@param ntargetnoc le n cible de l'allocation, en excluant la strate certain
                  Cet argument est de type double car dans le cas d'un CV cible on fournira un nombre
                  reel obtenu de la fonction get_nnoc_C. Cependant, dans le cas d'un n cible,
                  le ntargetnoc est le n cible auquel on soustrait Nc. Il s'agira donc a l'origine
                  d'un entier, que l'on devra convertir en double afin de pouvoir le donner en entree
                  a get_nhnonint_C.
voir la description generale au debut de ce fichier pour les autres parametres.

Aucune valeur retournee (ou modifiee) n'est unique a cette fonction, 
voir la description generale des sorties au debut de ce fichier.

Creee et verifiee le 27 septembre 2012
*/
void get_nhnonint_C (double *ntargetnoc, int *TCN, int *Nh, double *ah, int *L, int *takenone, int *takeall, 
                         /* Sortie */ double *nhnonint)
{
  int j;
  if (*takenone > 0){
    for (j=0; j < *takenone; j++){
		  nhnonint[j] = 0;
    }
  }
 	for (j=*takenone; j < (*L-*takeall); j++){
		if (Nh[j] == 0) nhnonint[j] = 0; else nhnonint[j] = (*ntargetnoc - *TCN) * ah[j];
 	}
  if (*takeall>0){
  	for (j=(*L-*takeall); j < *L; j++){
			nhnonint[j] = Nh[j];
  	}
  }
}

/* ******************************************************************************************** */

/* simple verification : doit-on faire un ajustement pour strate take-all?

Aucun parametre n'est unique a cette fonction, 
voir la description generale des parametres au debut de ce fichier.

@return takeall ici la valeur de takeall peut etre modifiee, en etant augmentee de 1
@return valid vaut 0 si on a besoin de faire un ajustement (la valeur de takeall a ete modifiee), 
                   1 si aucun ajustement n'est necessaire

Creee et verifiee le 28 septembre 2012
*/
void verif_takeall_C (double *nhnonint, int *Nh, int *L, int *takenone, 
                         /* Sortie et entree */ int *takeall,
                         /* Sortie */ int *valid)
{
  int j, Tta = 0;
  
  for (j=*takenone; j < (*L-*takeall); j++)
		if (nhnonint[j] > Nh[j]) Tta = Tta + 1;
    /* Attention, si j'ai des problemes, ca pourrait provenir de cette comparaison entre un int et un double */
    
  if ((Tta > 0) && (*takeall < *L - 1 - *takenone)){
    *takeall = *takeall + 1; 
    *valid = 0;     
  } else {
    *valid = 1; 
  }
/* J'ai change la condition nh[B]>=Nh[B] pour nh[B]>Nh[B] car on ne veux pas de depassement. Si nh[B]=Nh[B]
   on atteint le CV cible ou le n cible, c'est correct. Ce changement regle certains cas problematiques que j'avais
   rencontre avec Kozak qui restait pris dans un minimum global a cause d'un ajustement pour strate takeall
   qui n'avait pas lieu d'etre fait (dans premieres strates par exemple). Par contre, ca peut creer des resultats
   surprenant pour lesquels une strate qui parait recensement est presente en dessous d'une strate echantillonnee. */
}

/* ******************************************************************************************** */

/* Arrondissement des nhnonint afin d'obtenir les nh

Algorithme :
 Pas d'arrondissement a faire pour les strates takenone et takeall.
 Pour un CV cible : les nhnonint des strates takesome sont tous arrondis vers le haut.
 Pour un n cible (Il faut que la somme des nh arrondis = n, mais attention si on a une
 strate certain, le n cible pour l'arrondissement, donc la valeur donnee a l'agument n
 en entree, doit etre le n cible total - Nc, la taille de la strate certain.) :
 1- Ramener a 1 les nhnonint >0 et <1.
 2- Identifier les nhnonint qui reste a arrondir.
 3- Calculer le nombre de nhnonint a arrondir de chacune des facons possibles
    (nm1 : partie entiere - 1, mp0 : partie entiere, np1 : partie entiere + 1)
  Un nm1 non nul peut survenir quand des nhnonint entre 0 et 1 sont ramenes a 1.
 4- Ordonner les parties decimales de nhnonint qui reste a arrondir en ordre croissant.
 5- Arrondissement : 
    Les nhnonint avec les nm1 plus petites parties decimales seront arrondis par : partie entiere - 1,
    les nhnonint avec les np1 plus grandes parties decimales seront arrondis par : partie entiere + 1,
	les autres seront arrondis par : partie entiere.

@param ntargetround Utilise seulement si on a un n cible : il s'agit du n cible pour l'arrondissement.
                    C'est en fait ntargetnoc, le n cible demande - le nombre d'observations dans la strate certain,
                    donc un entier.
voir la description generale au debut de ce fichier pour les autres parametres.

Aucune valeur retournee (ou modifiee) n'est unique a cette fonction, 
voir la description generale des sorties au debut de ce fichier.

Creee le 25 avril 2012, ajustee et verifiee le 28 septembre 2012
*/ 
void get_nh_C(double *nhnonint, int *findn, int *ntargetround, int *Nh, int *TCN, 
                  int *L, int *takenone, int *takeall, 
  		            /* sortie */ double *nh)
{ 
  /* Note : Logiquement, nh devrait etre un entier. Cependant, les entiers ne peuvent pas prendre la
     valeur nan en C. Alors la programmation est plus simple si je stocke nh en double, em gardant en tete
     qu'en theorie il s'agit bien d'un entier */
  
	int j; /* pour les boucles */

 	/* Pas d'arrondissement a faire pour les strates takenone et takeall. */
  if (*takenone>0)
	  for (j=0; j < *takenone; j++) 
      nh[j] = 0;
	if (*takeall>0)
		for (j=(*L-*takeall); j < *L; j++)
			nh[j] = Nh[j];

	int Lts = *L - *takenone - *takeall; /* le nombre de strates takesome */
	if (*findn == 1) {
    /* Pour un CV cible : les nhnonint des strates takesome sont tous arrondis vers le haut.*/
		int nht[Lts]; /* parties entieres des nhnonint a arrondir */
		double reste[Lts]; /* parties decimales des nhnonint a arrondir  */	
		for (j=0; j < Lts; j++) {
		/* Ce code fait la meme chose que la fonction ceiling en R. */
    		nht[j] = floor(nhnonint[*takenone + j]);
    		reste[j] = nhnonint[*takenone + j] - nht[j]; 
    		if(reste[j]!=0) nh[*takenone + j] = nht[j] + 1; else nh[*takenone + j] = nht[j];
		}
	} else { 
    /* Pour un n cible : */
	
 		/* 1- Ramener a 1 les nhnonint >0 et <1 */
		int Ltsr = Lts; /* nombre de nhnonint a arrondir parmi les strates takesome */
		for (j=0; j < Lts; j++) {
	 		if( nhnonint[*takenone + j] > 0 && nhnonint[*takenone + j] < 1 ) {
				nh[*takenone + j] = 1; 
				Ltsr = Ltsr - 1;
			} 
		}
    
    /* Note : on laisse aller les nhnonint nuls ou negatifs. Il vont simplement ete arrondis comme les autres.
       Ces valeurs ne sont pas valides logiquement, mais elles sont possibles mathematiquement dans les cas suivants :
       Situation menant a un nhnonint nul : variance nulle dans une strate take-some
       Situation menant a des nhnonint negatifs : strate take-all plus grande que le n cible */
		
		if ( Ltsr > 0 ) { /* Si on a quelque chose a arrondir seulement */

			/* 2- Identifier les nhnonint qui doivent etre arrondis */
			int Id[Ltsr];
			int i=0; /* Indice pour le vecteur des valeurs a arrondir uniquement */
			for (j=0; j < Lts; j++) {
				if( !(nhnonint[*takenone + j] > 0 && nhnonint[*takenone + j] < 1) ) { 
					Id[i] = *takenone + j;
					i = i + 1;
				}
			}
   
			/* 3- Calculer le nombre de nhnonint a arrondir de chacune des facons possibles */
			int snht = 0; /* somme des parties entieres pour les nhnonint a arrondir */
			int nht[Ltsr]; /* parties entieres des nhnonint a arrondir */
			double reste[Ltsr]; /* parties decimales de nhnonint a arrondir  */	
			int index[Ltsr]; /* vecteur d'indices utile pour ordonner les parties decimales */
			for (j=0; j < Ltsr; j++) {
				nht[j] = floor(nhnonint[Id[j]]);
				reste[j] = nhnonint[Id[j]] - nht[j]; 
				snht = snht + nht[j];
				index[j] = j;			
			}
			int nhaut; /* nombre de nhnonint a arrondir vers le haut (peut etre negatif) */
			nhaut = (*ntargetround - (Lts - Ltsr) - snht - *TCN );
   
			int nm1, np0; /* (nm1 strates : partie entiere - 1, mp0 strates : partie entiere, 
							  np1 strates : partie entiere + 1 (pas besoin de le calculer)) */
			if ( nhaut < 0) {
				nm1 = -nhaut;
				np0 = Ltsr - nm1;
			} else {
				nm1 = 0;
				np0 = Ltsr - nhaut;
			}
			
			/* 4- Ordonner les parties decimales de nhnonint qui reste a arrondir en ordre croissant. */
			rsort_with_index(reste,index,Ltsr); /* si egalites, fait comme rank avec ties.method="first" */
			
			/* 5- Arrondissement final */   	
			for (j=0; j < Ltsr; j++) {
				if ( j < nm1) {
					nh[Id[index[j]]] = nht[index[j]] - 1;
				} else if ( j >= nm1 && j < nm1 + np0) {
					nh[Id[index[j]]] = nht[index[j]];
				} else {
					nh[Id[index[j]]] = nht[index[j]] + 1;
				}
				/* Si un nh est inferieur a 0, je vais le ramener a 0 */
				nh[Id[index[j]]] = fmax2(0, nh[Id[index[j]]]);
        /* Seul cas ou c'est possible : nnoc negatif, ce qui est vraiment tire par les chaveux.
           il est certain ue les ah sont positifs etant donne que les Nh et les VhY sont toujours des 
           nombres positifs, et que EYh est eleve au carre dans la formule des gammah. */
			}
		} 
	}
  
  /* Si les nhnonint prennaient la valeur nan, les nh cacules par le code ci-dessus ne prendront
     pas la valeur nan, mais plutot une valeur extreme. Alors je vais faire une petite correction ici */
  for (j=0; j < Lts; j++) {
    if (!R_FINITE(nhnonint[*takenone + j])) nh[*takenone + j] = R_NaN;
  }
  /* Valeurs manquantes en C, source = Wrinting R extensions, sections 5.10.3 et 6.4 */
  
  
/* J'ai etudie les valeurs possibles de nhaut.
   *ntargetround est le ncible. Il est en fait egal a  :
  0 [strates takenone] + sum(nhnonint) [strates takesome] + TCN [strates takeall].
	On peut briser sum(nhnonint) en deux parties : une pour les nhnonint entre 0 et 1, 
	et une autre partie pour les autres nhnonint.
	On soustrait a ca 
	0 [strates takenone] + ((Lts - Ltsr)*1 + snht) [strates takesome] + snhta [strates takeall].
	On voit que les termes pour les strates takenone et takeall s'annulent.
	ecrivons le resultat de la soustraction comme la somme des deux termes suivants :
	t1 = sum(nhnonint) [strates takesome avec nhnonint entre 0 et 1] - (Lts - Ltsr)
	t2 = sum(nhnonint) [strates takesome avec nhnonint pas entre 0 et 1] - snht
	ou, rappelons-le, snht est la somme des parties enteres des nhnonint pas entre 0 et 1
	pour une strate takesome.
	Les valeurs max et min de ces deux termes sont :
	min(t1) tend vers - (Lts - Ltsr), max(t1) tend vers 0;
	min(t2) tend vers 0, max(t2) tend vers Ltsr.
	Ainsi, min(nhaut) = min(t1) + min(t2) : tend vers - (Lts - Ltsr)
		   max(nhaut) = max(t1) + max(t2) : tend vers Ltsr.
	C'est donc dire qu'il est impossible d'avoir a arrondir vers le haut plus de nombres
	qu'on en a, ce qui est une excellente chose.
	Cependant, on pourrait avoir des cas limite comme
	nhnonhint = 0.3333 , 0.3333, 0.3333 : serait arrondi par 1, 1, 1, mais ici ncible = 1.
	Donc les nh fournis par notre algorithme ne somment pas au ncible.
	Je pense qu'on n'a pas ici a essayer de modifier l'algo d'arrondissement. 
	C'est simplement un signe que le ncible demande n'est pas realiste. */
}			

/* ******************************************************************************************** */

/* Calcul du RRMSE

Aucun parametre n'est unique a cette fonction, 
voir la description generale des parametres au debut de ce fichier.

@return RRMSE Le Relative Root Mean Squared Error pour le plan donne en entree. 

Creee le 27 avril 2012, ajustee et verifiee le 28 septembre 2012
*/
void get_RRMSE_C (double *biaspenalty, double *TY, double *TAY, int *Nh, double *VYh, double *nhcalcul, 
                  double *rhL, int *L,int *takenone,
                  /* sortie */ double *RRMSE)
{
  int j; /* pour les boucles */

  /* Si les nh sont nuls, RMSE prendra une valeur manquante */
  /* Cette situation est impossible dans strata_bh_opti_C avec dotests = 1 */
	for (j=*takenone; j < *L; j++) {
		if (nhcalcul[j] < 0)
		  *RRMSE = NA_REAL;	
	}
	if (!ISNA(*RRMSE)) { 
		double sum = 0; /* Somme a calculer sur les strates takesome et takeall */
		for (j=*takenone; j < *L; j++) {
      if (Nh[j] != 0 && VYh[j] != 0) sum = sum + R_pow(Nh[j],2) * VYh[j] * (R_pow(nhcalcul[j]*rhL[j],-1) - R_pow(Nh[j], -1));     
      /* Ici, on utilise R_pow(Nh[j], -1) plutot que 1/Nh[j] parce qu'avec le symbole /, etant
         donne que Nh[j] et 1 sont des entiers, le resultat est aussi stocke comme un entier. Ainsi la
         fraction Nh[j] est arrondie a 1 (si Nh[j]=1) ou 0 (sinon). Alors que R_pow transforme Nh[j] en 
         double avant de faire le calcul (def. de R_pow : double R_pow (double x, double y)).      
         R_pow(nhcalcul[j]*rhL[j],-1) pourrait, pour sa part, etre remplace par (1/(nhcalcul[j]*rhL[j]))
         parce que nhcalcul et rhL sont deja des double. */
		}
		*RRMSE = R_pow(fmax2(0, R_pow(*biaspenalty * *TAY, 2) + sum), 0.5) / *TY;
	} 
/* Valeurs manquantes en C, source = Wrinting R extensions, sections 5.10.3 et 6.4 */
}

/* ******************************************************************************************** */

/* Calcul de n, la taille d'echantillon totale, incluant la strate certain
   (il s'agit du critere a optimiser dans le cas d'un CV cible)

Aucun parametre n'est unique a cette fonction, 
voir la description generale des parametres au debut de ce fichier.

@return n la somme des elements de nhcalcul + Nc

Creee et verifiee le 28 septembre 2012.
*/
void get_n_C (double *nhcalcul, int *L, int *Nc,
                /* sortie */ double *n)
{
  int j;
  *n = *Nc; /* Il faut compter les unites de la strates certain dans cette somme */
  for (j=0; j < *L; j++) *n = *n + nhcalcul[j];
}

/* ******************************************************************************************** */

/* Teste les conditions sur les Nh 

Conditions : Nh >= minNh pour toutes les Ls strates take-some et take-all,
             Nh >= 0 pour les strate take-none

Aucun parametre n'est unique a cette fonction, 
voir la description generale des parametres au debut de ce fichier.

@return NhOK 1 si les conditions sur les Nh sont respectees pour toutes les strates, 0 sinon

Creee et verifiee le 12 octobre 2012.
*/
void test_Nh_C (int *Nh, int *L, int *takenone, int *minNh,
                  /* sortie */ int *NhOK)
{
  int j, cond = 0;
	for (j=0; j<*takenone; j++)    if( Nh[j] >= 0 )      cond = cond + 1; 
	for (j=*takenone; j < *L; j++) if( Nh[j] >= *minNh ) cond = cond + 1;
  if ( cond == *L) *NhOK = 1; /* 1 = les conditions sont respectees pour toutes les strates*/
  else *NhOK = 0;             /* 0 = les conditions ne sont pas respectee */
}


/* ******************************************************************************************** */

/* Teste les conditions sur les nh 

Conditions : nh > 0 pour toutes les Ls strates take-some et take-all
             (ce test sur les nh entiers revient a nh !=0 car on a ramene a 0 tous les nhnonint negatifs, 
              ce qui devrait survenir tres rarement, uniquement dans des cas limite),
             rien a tester sur les strates take-none car dans ce cas par definition nhnonint = nh = 0.

Aucun parametre n'est unique a cette fonction, 
voir la description generale des parametres au debut de ce fichier.

@return nhOK 1 si les conditions sur les nh sont respectees pour toutes les strates, 0 sinon

Creee et verifiee le 12 octobre 2012.
*/
void test_nh_C (double *nh, int *L, int *takenone,
                  /* sortie */ int *nhOK)
{
  int j, cond = 0;
	for (j=*takenone; j < *L; j++) if( nh[j] > 0 ) cond = cond + 1;
  if ( cond == *L - *takenone) *nhOK = 1; /* 1 = les conditions sont respectees pour toutes les strates*/
  else *nhOK = 0;                         /* 0 = les conditions ne sont pas respectee */
}


/* ******************************************************************************************** */

/* Fonction qui fait tous les calculs necessaires a l'algorithme de Kozak lorsqu'il "essaie une borne"

Cette fonction fournit les details d'un plan d'echantillonnage stratifie, pour des bornes donnees.
Elle travaille sur le vecteur des observations excluant la strate certain.

@param dotests indicatrice : 1 si on doit faire les tests des conditions sur les Nh et les nh 
                             (c'est le cas pour l'algo de Kozak, incluant l'enumeration complete),
                             0 si on n'a pas besoin de faire ces tests (strata.bh, .geo et .cumrootf)
@param takealladjust indicatrice : 1 si on demande a faire un ajustement automatique pour les strates recensement, 
                                   0 sinon
voir la description generale au debut de ce fichier pour les autres parametres.

@result takeallout un entier : valeur de takeall a la fin des calculs.
@result optinhnonint le critere a optimiser (n si CV cible, RRMSE si n cible) calcule a partir
                     des tailles d'echantillon non entieres nhnonint
@result optinh le critere a optimiser (n si CV cible, RRMSE si n cible) calcule a partir
               des tailles d'echantillon entieres nh
voir la description generale au debut de ce fichier pour les autres sorties.

Creee et verifiee le 28 septembre 2012
voir fonction R strata.bh.internal qui appelle cette fonction
*/
void strata_bh_opti_C(double *xnoc, int *Nnoc, double *bhfull, int *L, int *takenone, int *takeall,
            int *Nc, double *EYc, double *q1, double *q2, double *q3,
            int *nmodel, double *beta, double *sig2, double *ph, double *gamma, double *epsilon, double *EX, double *EX2,
            int *findn, int *n, double *CV, double *rhL, double *biaspenalty, int *takealladjust, int *dotests, int *minNh,
            /* Sortie */ int *NhOK, int *nhOK, /* uniquement utile a l'algo de Kozak*/ 
                         double *phih, double *psih, double *gammah, double *ah, double *U2, double *U, double *V, 
                         /* uniquement utile a l'algo de Sethi*/
                         int *stratumIDnoc, int *Nh, double *EYh, double *VYh, double *TY, double *TAY,  
                         double *nhnonint, int *takeallout, double *nh, double *optinhnonint, double *optinh)
{
  int j, TCN, valid = 0, arret1, arret2;
  double ntargetnoc;
  *takeallout = *takeall;
  
  /* Delimitation des strates */
  get_stratumIDnoc_C(xnoc, Nnoc, bhfull, L, stratumIDnoc);
  
  /* Calcul des tailles de population dans les strates*/
  get_Nh_C(stratumIDnoc, Nnoc, L, Nh);

/*Rprintf("Nh : "); for (j=0; j < *L; j++) Rprintf("%d  ",Nh[j]);  Rprintf("\n");*/

  /* Si besoin, test sur les Nh ici */
  if (*dotests == 1) {
    test_Nh_C(Nh, L, takenone, minNh, NhOK);
    if (*NhOK == 0) arret1 = 1; else arret1 = 0; /* On arrete l'execution de la fonction si NhOK=0 */
  } else {
    arret1 = 0;
  }

/*Rprintf("arret1 : %d  ",arret1);  Rprintf("\n");*/
 
  if (arret1 == 0) {
    /* Tout le reste de la fonction est dans ce if car si on a teste les Nh et qu'ils ne respectaient pas les
       conditions, il est inutile de faire les calculs subsequents (donc on sauve du temps en ne les faisant pas,
       mais attention, les valeurs en sortie autres que NhOK ne seront pas bonnes, il ne faut pas les considerer) */
    
    /* Calcul des moments anticipees (EYh et VYh) et de sommes dont on a besoin dans les calculs subsequents (TY et TAY)*/
    get_momentY_C(xnoc, stratumIDnoc, Nnoc, Nh, L, Nc, EYc, takenone, nmodel, beta, sig2, ph, gamma, epsilon, 
                      EX, EX2, EYh, VYh, phih, psih, TY, TAY);

/*Rprintf("EYh : "); for (j=0; j < *L; j++) Rprintf("%f  ",EYh[j]);  Rprintf("\n");
Rprintf("VYh : "); for (j=0; j < *L; j++) Rprintf("%f  ",VYh[j]);  Rprintf("\n");*/

    /* Allocation : repartition de n parmi les L strates. */
    /* Etape 1 alloc : calcul des gammah */
    get_gammah_C(Nh, EYh, VYh, L, q1, q2, q3, gammah);
    
/*Rprintf("gammah : "); for (j=0; j < *L; j++) Rprintf("%f  ",gammah[j]);  Rprintf("\n");*/

    while (valid == 0) {
      /* Etape 2 alloc : calcul des ah et de TCN*/
      get_ah_C(gammah, Nh, L, takenone, takeallout, ah, &TCN, U2);
      
/*Rprintf("ah : "); for (j=0; j < *L; j++) Rprintf("%f  ",ah[j]);  Rprintf("\n");*/

      /* Etape 3 alloc : calcul du ntargetnoc*/
      if (*findn == 1){
        /* Dans le cas d'un CV cible, on doit appliquer une formule */
        get_nnoc_C(CV, Nh, ah, VYh, rhL, L, takenone, takeallout, biaspenalty, TAY, TY, &TCN, &ntargetnoc, U, V);
      } else {
        ntargetnoc = *n - *Nc; 
      }
      
/*Rprintf("ntargetnoc : %f  ",ntargetnoc);  Rprintf("\n");*/

      /* Etape 4 alloc : application de la regle d'allocation 
                         (qui donnera en resultat les tailles d'echantillon non entieres nhnonint) */
      get_nhnonint_C(&ntargetnoc, &TCN, Nh, ah, L, takenone, takeallout, nhnonint);

/*Rprintf("nhnonint : "); for (j=0; j < *L; j++) Rprintf("%f  ",nhnonint[j]);  Rprintf("\n");*/

      /* Etape 5 alloc : verification strates recensement, si besoin */
      if (*takealladjust == 1){
        verif_takeall_C(nhnonint, Nh, L, takenone, takeallout, &valid);
      } else {
        valid = 1;
      }

/*Rprintf("valid : %d  ",valid);  Rprintf("\n");*/


    }
    /* Etape 6 alloc : Arrondissement des nhnonint pour obtenir les nh entiers */
    int ntargetround = ntargetnoc;
    get_nh_C(nhnonint, findn, &ntargetround, Nh, &TCN, L, takenone, takeallout, nh);

/*Rprintf("nh : "); for (j=0; j < *L; j++) Rprintf("%d  ",nh[j]);  Rprintf("\n");*/

       /* Ici, si on fait un ajustement automatique pour strates recensement, on ne ramene pas a Nh les nh > Nh
       (nh entiers ou reels), car le but de cette fonction est de calculer les criteres pour 
       l'optimisation. La fonction sera utilisee par l'algorithme de Kozak. 
       Si on rameneait a Nh les nh > Nh dans l'optimisation, on aurait le probleme suivant.
       Dans le cas d'un CV cible, le critere a minimiser et n. Si on tronque a Nh les valeurs de nh, l'apparition
       d'une strate recensement causerait une diminution du critere n, qu'on cherche justement a minimiser. Alors
       les plans faisant apparaitre une strate recensement seraient favorises. L'optimisation prendrait alors une 
       mauvaise direction. De toute facon, l'ajustement automatique pour strate recensement est obligatoire a 
       toutes les fonctions qui calculent des bornes. Cet ajustement nous assure que tous les nh (entiers ou reels)
       seront <= Nh au terme des calculs. */     
    if (*takealladjust == 0 || *takeallout == *L - 1 - *takenone)
      for (j=*takenone; j < (*L-*takeall); j++)
        if (nh[j] > Nh[j]) nh[j] = Nh[j];
    /* Si l'ajustement automatique n'est pas fait ou s'il y a une seule strate take-some, il se pourrait qu'un 
       nh soit > Nh a la sortie de la presente fonction. Dans ce cas, on va ramener a Nh les nh entiers > Nh, 
       mais pas les nh reels pour que l'utilisateur puisse comprendre a partir de la sortie pourquoi le CV cible 
       n'est pas necessairement atteint sans correction pour strate recensement. */
        
    /* Si besoin, test sur les nh ici */
    if (*dotests == 1) {
      test_nh_C(nh, L, takenone, nhOK);
      if (*nhOK == 0) arret2 = 1; else arret2 = 0; /* On arrete l'execution de la fonction si nhOK=0 */
    } else {
      arret2 = 0;
    }
    
/*Rprintf("arret2 : %d  ",arret2);  Rprintf("\n");*/
  
    if (arret2 == 0) {
      /* Tout le reste de la fonction est dans ce if car si on a teste les nh et qu'ils ne respectaient pas les
         conditions, il est inutile de faire les calculs subsequents (donc on sauve du temps en ne les faisant pas, mais
         attention, les valeurs en sortie optinhnonint et optinh ne seront pas bonnes, il ne faut pas les considerer) */
      
      /* Calcul du critere a optimiser, sur les nh et sur les nhnonint */
      if (*findn == 1) { /* Pour un CV cible, le critere a optimiser est n, la taille d'echantillon totale */
        get_n_C(nhnonint, L, Nc, optinhnonint);
        get_n_C(nh, L, Nc, optinh);
      } else { /* Pour un n cible, le critere a optimiser est le RRMSE */
        get_RRMSE_C(biaspenalty, TY, TAY, Nh, VYh, nhnonint, rhL, L, takenone, optinhnonint);
        get_RRMSE_C(biaspenalty, TY, TAY, Nh, VYh, nh, rhL, L, takenone, optinh);
      }
      
/*Rprintf("optinh : %f  ",optinh);  Rprintf("optinhnonint : %f  ",optinhnonint); Rprintf("\n");*/
    }
  }
}


/* ******************************************************************************************** */

/* Makes the conversion from stratum boundaries expressed in terms of data rank (pbh), 
   to stratum boundaries expressed on the scale of the data (bhfull)

@param pbh vecteur de longueur L-1 representant des bornes de strates, mais sur l'echelle des rangs des donnees :
           chaque element de pbh est un entier representant la position dans le vecteur x1noc d'une borne.
voir la description generale au debut de ce fichier pour les autres parametres.

@return bhfull le vecteur L+1 des bornes pleines, sur l'echelle des donnees, equivalentes aux pbh

Creee le 12 octobre 2012
*/
void pbh2bhfull_C(int *pbh, int *L, double *x1noc, int *N1noc, double *bhfull)
{ 
  int j;
  bhfull[0] = x1noc[0];
  for (j=0; j < *L-1; j++) {
		if(pbh[j]<=1) {
			bhfull[j+1] = x1noc[0];
		} else if (pbh[j]>*N1noc) {
			bhfull[j+1] = x1noc[*N1noc - 1] + 1;
		} else 
			bhfull[j+1] = (x1noc[pbh[j]-1] + x1noc[pbh[j]-2])/2;
	}
  bhfull[*L] = x1noc[*N1noc - 1] + 1;  
} 


/* ******************************************************************************************** */

/* Enumeration complete des ensembles de bornes possibles pour touver le meilleur 

On essaie tous les ensembles de bornes possibles et on recherche la meilleure, celle qui minimise le critere (RRMSE si n cible, n si CV cible).

@param pbhsol vecteur de longueur (L-1)*nsol qui contient tous les ensembles de bornes a essayer, 
              representees sur l'echelle des rangs : les L-1 premiers elements forment le premier ensemble de bornes,
              les L-1 elements suivants le deuxieme ensemble, etc.
@param nsol le nombre d'ensemble de bornes a essayer
voir la description generale au debut de ce fichier pour les autres parametres.

@return soldetail un enorme vecteur de longueur ((L-1)+2*L+5)*nsol, qui devriendra une matrice en R,
                  comprennant les elements dont on veut garder la trace (bh, Nh, nh, takeall, Nhok, nhok, 
                  optinh, optinhnonint) pour chaque ensemble de bornes a essayer. Ordre dans le vecteur :
                  tous les elements pour le premier ensemble de bornes, suivi de tous les elements pour le 
                  deuxieme, etc.

Creee le 12 octobre 2012
*/
void complete_enum_C(int *pbhsol, int *nsol, int *L, double *x1noc, int *N1noc, double *xnoc, int *Nnoc, 
                         int *takenone, int *takeall, int *Nc, double *EYc, double *q1, double *q2, double *q3,
                         int *nmodel, double *beta, double *sig2, double *ph, double *gamma, double *epsilon, 
                         double *EX, double *EX2, int *findn, int *n, double *CV, double *rhL, 
                         double *biaspenalty, int *minNh,
                         /* Sortie */ double *soldetail)
{
  int i, j, pbh[*L-1], NhOK, nhOK, stratumIDnoc[*Nnoc], Nh[*L], takeallout;
  double bhfull[*L+1], EYh[*L], VYh[*L], TY, TAY, nh[*L] , nhnonint[*L], optinhnonint, optinh;
  double phih[*L], psih[*L], gammah[*L], ah[*L], U2, U, V;
  int takealladjust = 1; /* Ici on fait l'ajustement pour strates takeall, */
  int dotests = 1;       /* ainsi que les tests sur Nh et nh. */

  for (i=0; i < *nsol; i++){
    /* Determiner les bornes pleines sur l'echelle des donnees a partir des combinaisons fournies en entree */
    for (j=0; j < *L - 1; j++) pbh[j] = pbhsol[j + i * (*L - 1)];
    pbh2bhfull_C(pbh, L, x1noc, N1noc, bhfull);

    /* Calculs pour la stratification */
    strata_bh_opti_C(xnoc, Nnoc, bhfull, L, takenone, takeall, Nc, EYc, q1, q2, q3,nmodel, beta, sig2, ph, 
           gamma, epsilon, EX, EX2, findn, n, CV, rhL, biaspenalty, &takealladjust, &dotests, minNh,
           &NhOK, &nhOK, phih, psih, gammah, ah, &U2, &U, &V, stratumIDnoc, Nh, EYh, VYh, &TY, &TAY, 
           nhnonint, &takeallout, nh, &optinhnonint, &optinh);
           
    /* Enregistrement des resultats (bh, Nh, nh, takeall, Nhok, nhok, optinh, optinhnonint)*/
    for (j=0; j < *L - 1; j++) soldetail[j + i * (3 * *L + 4)] = bhfull[j + 1];
    for (j=0; j < *L; j++) soldetail[*L - 1 + j + i * (3 * *L + 4)] = Nh[j];
    for (j=0; j < *L; j++) soldetail[2 * *L - 1 + j + i * (3 * *L + 4)] = nh[j];
    soldetail[3 * *L - 1 + i * (3 * *L + 4)] = optinh;
    soldetail[3 * *L + i * (3 * *L + 4)]     = optinhnonint;
    soldetail[3 * *L + 1 + i * (3 * *L + 4)] = takeallout;
    soldetail[3 * *L + 2 + i * (3 * *L + 4)] = NhOK;
    soldetail[3 * *L + 3 + i * (3 * *L + 4)] = nhOK;
  }
  
}


/* ******************************************************************************************** */

/* Algorithme de Kozak (original) 

@param combin2try vecteur de longueur (2*(L-1)+6)*ncombin qui contient toutes les combinaisons a essayer.
                  Le vecteur contient dans l'ordre
                  elements : (ibhtype,    ibh,    ipbh, optinh, optinhnonint, takeall, maxstep, maxstill)
                  positions: (      0, 1->L-1, L->2L-2,   2L-1,           2L,    2L+1,    2L+2,     2L+3)
                  definissant la premiere combinaison, tous les memes elements pour la deuxieme combinaison, etc.
@param ncombin le nombre de combinaisons a essayer (le nombre de lignes dans la version matrice de combin2try)
voir la description generale au debut de ce fichier pour les autres parametres.

@return rundetail  un enorme vecteur de longueur (2*(L-1)+8)*(ncombin*rep), qui devriendra une matrice en R,
                   comprennant les elements dont on veut garder la trace pour chaque combinaison a essayer :
                   elements : (    bh, optinh, optinhnonint, takeall, niter, ibhtype,       ibh,  rep)
                   positions: (0->L-2,    L-1,            L,     L+1,   L+2,     L+3, L+4->2L+2, 2L+3)
                   Ordre dans le vecteur : tous les elements pour la premiere combinaison, 
                   suivi de tous les elements pour la deuxieme, etc.
@return iterdetail un enorme vecteur de longueur ((L-1)+6)*nrowiterdetail, qui devriendra une matrice en R,
                   comprennant les elements dont on veut garder la trace pour chaque iteration de l'algo :
                   elements : (    bh, optinh, optinhnonint, takeall, step, iter, run)
                   positions: (0->L-2,    L-1,            L,     L+1,  L+2,  L+3, L+4)
                   Ordre dans le vecteur : tous les elements pour la premiere iteration, 
                   suivi de tous les elements pour la deuxieme, etc.
@return nrowiterdetail le nombre de lignes dans la version matrice du vecteur iterdetail

Creee le 17 et 18 octobre 2012
*/
void algo_Kozak_C(double *combin2try, int *ncombin, int *L, double *x1noc, int *N1noc, double *xnoc, int *Nnoc, 
                      int *takenone, int *takeall, int *Nc, double *EYc, double *q1, double *q2, double *q3,
                      int *nmodel, double *beta, double *sig2, double *ph, double *gamma, double *epsilon, 
                      double *EX, double *EX2, int *findn, int *n, double *CV, double *rhL, 
                      double *biaspenalty, int *minNh, int *maxiter, int *idoptinh, int *rep,
                      /* Sortie */ double *rundetail, double *iterdetail, int *nrowiter)
{
  int r, i, j, iter, istill, accept, maxstep, maxstill, jbh, step, sign;
  int pbh[*L-1], npbh[*L-1], NhOK, nhOK, stratumIDnoc[*Nnoc], Nh[*L], takeallout;
  double bhfull[*L+1], EYh[*L], VYh[*L], TY, TAY, phih[*L], psih[*L], gammah[*L], ah[*L], U2, U, V;
  double nh[*L], nhnonint[*L], optinhnonint, optinh, noptinhnonint, noptinh, diffrelopti;
  int takealladjust = 1; /* Ici on fait l'ajustement pour strates takeall, */
  int dotests = 1;       /* ainsi que les tests sur Nh et nh. */
  int irow = 0;  /* Le compteur du nombre de lignes dans iterdetail est initialise a zero. */
  int rrow = 0;  /* Le compteur du nombre de lignes dans rundetail est initialise a zero. */
  int ncol_combin = 2 * *L + 4; /* Le nombre de colonnes dans la version matrice du vecteur combin2try */
  int ncol_run = 2 * *L + 4;    /* Le nombre de colonnes dans la version matrice du vecteur rundetail */
  int ncol_iter = *L + 5;       /* Le nombre de colonnes dans la version matrice du vecteur iterdetail */
  
  for (i=0; i < *ncombin; i++){  /* boucle sur les combinaisons a essayer */
  
    for (r=0; r < *rep; r++){  /* boucle sur les repetitions de l'algorithme */   
      
      /* Variables a initialiser pour l'algorithme */
      /* Rappel utile a la programmation ici -> combin2try : elements et leur positions
         (ibhtype,    ibh,    ipbh, optinh, optinhnonint, takeall, maxstep, maxstill)
         (      0, 1->L-1, L->2L-2,   2L-1,           2L,    2L+1,    2L+2,     2L+3) */
      iter = 0;   /* Le compteur d'iterations par repetitions est initialise a zero. */
      istill = 0; /* Le compteur d'iterations sans acceptation de nouvelles bornes est initialise a zero. */
      for (j=0; j < *L - 1; j++)  /* bornes initiales sur l'echelle des rangs */
        pbh[j] = combin2try[(j + *L) + i * ncol_combin];
      optinh       = combin2try[(2 * *L - 1) + i * ncol_combin];
      optinhnonint = combin2try[(2 * *L + 0) + i * ncol_combin];
      maxstep      = combin2try[(2 * *L + 2) + i * ncol_combin];
      maxstill     = combin2try[(2 * *L + 3) + i * ncol_combin];
      

      /* Info a mettre dans iterdetail pour les bornes initiales */
      /* Rappel utile a la programmation ici -> iterdetail : elements et leur positions
         (    bh, optinh, optinhnonint, takeall, step, iter, run)
         (0->L-2,    L-1,            L,     L+1,  L+2,  L+3, L+4) */
      for (j=0; j < *L - 1; j++)                     /* bornes sur l'echelle des donnees */
        iterdetail[j + irow * ncol_iter] = combin2try[(j + 1) + i * ncol_combin];
      for (j=*L - 1; j < *L + 2; j++)                /* optinh, optinhnonint et takeall */
        iterdetail[j + irow * ncol_iter] = combin2try[(j + *L) + i * ncol_combin];  
      iterdetail[(*L + 2) + irow * ncol_iter] = 0;    /* pas fait par l'algo de Kozak */
      iterdetail[(*L + 3) + irow * ncol_iter] = iter; /* numero de l'iteration */
      iterdetail[(*L + 4) + irow * ncol_iter] = rrow + 1; /* numero du run de l'algo */
      /* On vient de remplir une ligne de iterdetail, il faut donc incrementer de 1 irow. */
      irow = irow + 1;

      /* On peut deja mettre les info pour les bornesinitiales dans dans rundetail  */
      /* Rappel utile a la programmation ici -> rundetail : elements et leur positions
         (    bh, optinh, optinhnonint, takeall, niter, ibhtype,       ibh,  rep)
         (0->L-2,    L-1,            L,     L+1,   L+2,     L+3, L+4->2L+2, 2L+3) */
      rundetail[(*L + 3) + rrow * ncol_run] = combin2try[0 + i * ncol_combin];
      for (j=0; j < *L-1; j++) rundetail[(j + *L + 4) + rrow * ncol_run] = combin2try[(j + 1) + i * ncol_combin];      
      rundetail[(2 * *L + 3) + rrow * ncol_run] = r + 1;      

      /* Iterations de l'algorithme */      
    	while ( (iter < *maxiter) && (istill < maxstill) ) {	
        
  		  /* Identifier le nouvel ensemble de bornes a essayer */
        /* Choix aleatoire de la borne a remplacer*/
    		GetRNGstate(); 
    		jbh = floor(unif_rand()* (*L-1)); /* 0, 1, ..., L-2 */
    			if (jbh==(*L-1)) jbh=0; /* si jbh vaut L-1, ce qui a une probabilite presque nulle */
        /* Choix aleatoire du pas a faire*/
    		step =  floor(unif_rand() * maxstep) + 1; /* 1, 2, ..., maxstep */
    			if (step==(maxstep+1)) step=1; /* si step vaut maxstep+1, ce qui a une probabilite presque nulle */
    		sign = floor(unif_rand()*2); /* 0 ou 1 */
    		if (sign==0) step=-step;
    		PutRNGstate();    		
    		/* Obtention des nouvelles bornes en modifiant les anciennes */
    		for (j=0; j < *L-1; j++) npbh[j] = pbh[j];
    		npbh[jbh] = pbh[jbh] + step;
        
/*Rprintf("%d  : ",iter); 
for (j=0; j < *L-1; j++) Rprintf("%d  ",npbh[j]); */
        
        /* Determination des bornes pleines sur l'echelle des donnees associees a ces bornes modifiees */
        pbh2bhfull_C(npbh, L, x1noc, N1noc, bhfull);
    
        /* Calculs pour la stratification */
        strata_bh_opti_C(xnoc, Nnoc, bhfull, L, takenone, takeall, Nc, EYc, q1, q2, q3,nmodel, beta, sig2, ph, 
               gamma, epsilon, EX, EX2, findn, n, CV, rhL, biaspenalty, &takealladjust, &dotests, minNh,
               &NhOK, &nhOK, phih, psih, gammah, ah, &U2, &U, &V, stratumIDnoc, Nh, EYh, VYh, &TY, &TAY, 
               nhnonint, &takeallout, nh, &noptinhnonint, &noptinh);
               
        /* Acceptation ou rejet des nouvelles bornes */
        if(NhOK == 1 && nhOK == 1){
    			/* test sur le critere a optimiser (n ou RRMSE) : a-t-il diminue? */
  				if (*idoptinh) { /* Si on veux calculer le critere sur les nh entiers */
            if (!R_FINITE(noptinh)){
              accept = 0;
            } else if ( noptinh < optinh) { /* si noptinh est positif et plus petit que *optinh, */
  						accept = 1;             /* on change les bornes */
  					} else if (noptinh == optinh) { /* si noptinh est egale a *optinh, */
  					    /* on va comparer les critere calcule sur les nhnonint */       
  						if (noptinhnonint < optinhnonint) accept = 1; else accept = 0;
  					} else {     /* si noptinh est plus grand que optinh, */					
  						accept = 0;  /* c'est certain qu'on ne change pas les bornes */
  					} 
  				} else { /* Si on veux calculer le critere sur les nhnonint (non entiers) */
            if (!R_FINITE(noptinhnonint)){
              accept = 0;
            } else {
              if ( noptinhnonint < optinhnonint) accept = 1; else accept = 0;
            }
  				}
  				/* fin du test sur le critere a optimiser */          
        } else {
          /* Si les conditions sur les Nh et les nh ne sont pas respectees, on n'accepte pas les nouvelles bornes */
          accept = 0;
        }

/*Rprintf(" : %d ",accept); 
Rprintf("\n");*/

      	/* Actions posees dependamment de l'acceptation ou du rejet des nouvelles bornes */
    		if(accept==1) {
          /* On verifie d'abord si on doit changer la valeur de maxstep et maxstill,
             car lorsque l'algorithme se raproche de la convergence, on veut un petit maxstep
             afin de minimiser le temps de calcul. */
            diffrelopti = (optinhnonint - noptinhnonint)/optinhnonint;
            if (istill >= 50 && (step <= 3 && step >= -3) && diffrelopti < 0.001){ 
            /* La regle arbitraire que je me suis donne est donc : 
               si on a du faire plusieurs (plus de 50) iterations sans changement avant d'accepter 
               les nouvelles bornes et
               si le pas accepte est petit (inferieur a 3) en valeur absolue et
               si le changement relatif dans le critere d'optimisation sur les nhnonint est petit 
               (inferieur a 0.001),
               alors on change le maxstep pour 3 et le maxstill pour 50 (suggestions initiales de Kozak
               dans son article de 2006 dans Techniques d'enquete).*/
            maxstep = 3;
            maxstill = 50;
          }
          /* Le compteur d'iterations sans acceptation de nouvelles bornes est reinitialise a zero. */
    			istill = 0;
          /* Les bornes et les criteres d'optimisation doivent etre mis a jour */
    			for (j=0; j < *L-1; j++) pbh[j] = npbh[j];
    			optinh = noptinh;
    			optinhnonint = noptinhnonint;
          /* On doit enregistrer les resultats dans iterdetail */
          /* Rappel utile a la programmation ici -> iterdetail : elements et leur positions
             (    bh, optinh, optinhnonint, takeall, step, iter, run)
             (0->L-2,    L-1,            L,     L+1,  L+2,  L+3, L+4) */
      		for (j=0; j < *L-1; j++) iterdetail[j + irow * ncol_iter] = bhfull[j+1];
    			iterdetail[(*L - 1) + irow * ncol_iter] = optinh;
    			iterdetail[(*L + 0) + irow * ncol_iter] = optinhnonint;
      		iterdetail[(*L + 1) + irow * ncol_iter] = takeallout;
    			iterdetail[(*L + 2) + irow * ncol_iter] = step;
    			iterdetail[(*L + 3) + irow * ncol_iter] = iter + 1;
      		iterdetail[(*L + 4) + irow * ncol_iter] = rrow + 1;
          /* On vient de remplir une ligne de iterdetail, il faut donc incrementer de 1 irow. */
        	irow = irow + 1;
    		} else {
          /* Le compteur d'iterations sans acceptation de nouvelles bornes est incremente de 1. */
      	  istill = istill + 1 ;
          /* Il n'y a pas d'autres actions a poser. */
    		}
        
        /* Fin de l'iteration */
        iter = iter + 1;
        
    	}
      
      /* On enregistre les resultats pour les bornes finales dans rundetail  */
      /* Rappel utile a la programmation ici -> rundetail : elements et leur positions
         (    bh, optinh, optinhnonint, takeall, niter, ibhtype,       ibh,  rep)
         (0->L-2,    L-1,            L,     L+1,   L+2,     L+3, L+4->2L+2, 2L+3) */
      /* En premier on tire  bh, optinh, optinhnonint et takeall de iterdetail */
      for (j=0; j < *L+2; j++) rundetail[j + rrow * ncol_run] = iterdetail[j + (irow - 1) * ncol_iter];
      rundetail[(*L + 2) + rrow * ncol_run] = iter;      

/*for (j=0; j < ncol_run; j++) Rprintf("%f  ",rundetail[j + rrow * ncol_run]); Rprintf("\n");*/

      /* On vient de remplir une ligne de rundetail, il faut donc incrementer de 1 rrow. */
      rrow = rrow + 1;

    }
    
  }
  
  /* Au terme des calculs on connait le nombre de lignes de la version matrice de iterdetail */
  *nrowiter = irow;
  
}


