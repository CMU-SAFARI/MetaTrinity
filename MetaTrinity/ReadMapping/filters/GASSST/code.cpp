/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file code.h
 * \brief Module Code, responsable du codage des séquences de nucléotides
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */
 
 #include "code.h"
 
 #include <math.h>
 #include "constants.h"


/** 
 * Méthode qui calcule la valeur du ScoreBit d'un alignement en fonction de son score
 * \param s le score de l'alignement
 * \return la valeur ScoreBit de l'alignement
 */
double ScoreBit(double s)
{
  double x,y,z;

  x = s * LAMBDA;
  y = log(CONST_K);
  z = x-y;
  z = z / log(2);
  return z;
}

/**
 * Méthode calculant la E-value des banques pour un certain score
 * \param s un entier représentant un score
 * \param n1 la taille de la première banque
 * \param i1 le nombre de sequences dans la première banque
 * \param n2 la taille de la ceconde banque
 * \param i2 le nombre de séquences dans la seconde banque
 */
double ComputeEvalue(int s,long long  int n1, int i1, int n2, int i2)
{
  double e,x,y,x1,x2,mn;
  int    l;

  /// Calcul du effectif HSP length
  /// l = log(Kmn)/H

  x1 = (double) n1;
  x2 = (double) n2;
  x = x1 * x2;
  x1 = (double) CONST_K;
  x = x * x1;
  x = log (x);
  x1 = CONST_H;
  x = x / x1;
  l = (int) x;

  /// On calcule l'espace de recherche

  x1 = (double) n1 - i1*l;
  x2 = (double) n2 - i2*l;
  mn = x1*x2;

  /// On calcule la E-value
  /// e = K mn exp (-LAMBDA * s)

  x  = (double) mn;
  x1 = (double) CONST_K;
  x = x * x1;
  x2 = (double) LAMBDA;
  x1 = (double) s;
  y = x1 * x2;
  y = -y;
  y = exp(y);

  e = x * y;
  
  return e;
}



//donne estimation de 'bonne' taille de graine
//len : taille de la sequence
//nerr : nombre derreurs max autorisees
//neval : nb de simu pour evaluation stat
//param quantile: on vire les 'nquantile ' premiers res :
// correspond au quantile  (nquantile/neval) %
int find_seed_len(int len, int nerr, int neval,int nquantile)
{


  int posalea,i,j,minseed,eprec,seedval;
  
  char * taberr;
  taberr = (char *) malloc(sizeof(char)*(len+1));
  std::list<int> listtir;

  seedval=10*len;
  for (j=0; j<neval; j++)
    {
      memset(taberr, 0, len+1);
      for(i=0; i< nerr; i++)
	{
	  posalea = 1 + (int) ( ((float)(len) )* (rand() / (RAND_MAX + 1.0)));  // tirage dans [1; len]
	  taberr[posalea]=1;
	}
      minseed = 0;
      for(i=1; i<= len; i++)
	{
	  eprec = i-1;
    
	  while (taberr[i]==0 && i <len) i++ ;
	  minseed = max(minseed,i-eprec-1) ;
	}

      listtir.push_back(minseed);
      seedval =  min(seedval,minseed);
    }


  listtir.sort();

  std::list <int>::iterator ii = listtir.begin();
  ii= listtir.begin(); 
  std::advance(ii, nquantile);


  free(taberr); 
  
  return(*ii);
  
}

