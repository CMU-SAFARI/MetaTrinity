/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Alignment.cpp
 * \brief Class Alignment, d�finissant un alignement entre deux s�quences
 * \author Damien Fleury
 * \author Guillaume Rizk
 * \date 15/07/10
 */


#include "Alignment.h"

#include <iostream>
#include "constants.h"
#include "misc.h"
#include "code.h"
#include "display.h"
#include "withgap.h"
#include <string.h>


/**
 * Le consructeur par d�faut d'Alignment
 */
Alignment::Alignment() : length(0), nb_mismatches(0), nb_gaps(0),etat_cigar(0),cpt_cigar(0), cigar_index(0)
{
	sequence1 = NULL;
	sequence2 = NULL;
	cigar = NULL;
}

/**
 * Un consructeur d'Alignment
 * \param size la taille de l'alignement
 */
Alignment::Alignment(int size) : length(0), nb_mismatches(0), nb_gaps(0),etat_cigar(0),cpt_cigar(0), cigar_index(0),e_value(0)
{
	sequence1 = new char[size+1];
	sequence2 = new char[size+1];
	cigar = new char[4*size+1];
}

/**
 * La fonction d'initialisation d'un Alignment
 */
void Alignment::init()
{
	length = 0;
	nb_mismatches = 0;
	nb_gaps = 0;
	start1 = 0;
	end1 = 0;
	start2 = 0;
	end2 = 0;
	e_value = 0;
	j1 = NULL;
	j2 = NULL;
	rev_comp = 0;
	sequence1[length] = '\0';
	sequence2[length] = '\0';
	cigar[length] = '\0';
	cpt_cigar = 0;etat_cigar=0; cigar_index=0;

}


/**
 * La fonction d'initialisation d'un Alignment
 */
void Alignment::init(char rever)
{
	length = 0;
	nb_mismatches = 0;
	nb_gaps = 0;
	start1 = 0;
	end1 = 0;
	start2 = 0;
	end2 = 0;
	e_value = 0;
	j1 = NULL;
	j2 = NULL;
	rev_comp = 0;
	sequence1[length] = '\0';
	sequence2[length] = '\0';
	rev_comp=rever;
	cigar[length] = '\0';
	cpt_cigar = 0;etat_cigar=0;cigar_index=0;
}


/**
 * Le consructeur par recopie d'Alignment
 * \param al un objet Alignment
 */
Alignment::Alignment(const Alignment& al)
{
	*this = al;
}

/**
 * L'op�rateur d'affectation d'Alignment
 */
Alignment& Alignment::operator=(const Alignment& al)
{
	if(this != &al)
	{
		length = al.length;
		nb_mismatches = al.nb_mismatches;
		nb_gaps = al.nb_gaps;
		sequence1 = new char[length+20];
		sequence2 = new char[length+20];
		strncpy(sequence1, al.sequence1, length); 
		strncpy(sequence2, al.sequence2, length);
		sequence1[length] = '\0';
		sequence2[length] = '\0';

		int lencigar = strlen(al.cigar) +1; // +1 pour '\0'
		cigar = new char[lencigar];
		strncpy(cigar, al.cigar, lencigar);

		start1 = al.start1;
		start2 = al.start2;
		end1 = al.end1;
		end2 = al.end2;

		e_value = al.e_value;
		j1 = al.j1;
		j2 = al.j2;
		n1 = al.n1;
		n2 = al.n2;
		rev_comp = al.rev_comp;
		cpt_cigar = al.cpt_cigar;
		etat_cigar = al.etat_cigar;
		cigar_index = al.cigar_index;

	}
	return *this;
}

/**
 * Le destructeur par d�faut d'Alignment
 */
Alignment::~Alignment()
{
	delete [] sequence1;
	delete [] sequence2;
	delete [] cigar;
}

/**
 * M�thode d'ajout d'un couple d'un caract�re aux s�quences
 * \param c1 le caract�re � ajouter � la premi�re s�quence
 * \param c2 le caract�re � ajouter � la seconde s�quence
 */
void Alignment::addPair(char c1, char c2)
{

	if(OUTPUT_FORMAT==SAM_READY_FORMAT) 
	  {
	    //ajout  calcul string cigar
	    char  etat_courant;
	    if(c1==CHAR_GAP) etat_courant='I';
	    else if(c2==CHAR_GAP) etat_courant='D';
	    else etat_courant='M';

	    // cas initialisation
	    if(etat_cigar == 0) {cpt_cigar=1; etat_cigar=etat_courant;}
	    else if (etat_cigar == etat_courant) cpt_cigar++;
	    else // chgt d'etat
	      {
		cigar_index += sprintf(cigar + cigar_index ,"%i%c",cpt_cigar,etat_cigar);
		cpt_cigar=1; etat_cigar = etat_courant;
	      }
	
	    if(c2!=CHAR_GAP)
	      {
		sequence2[length] = c2;
		length++;
		sequence2[length] = '\0';
	      }
	    return;
	  }

	sequence1[length] = c1;
	sequence2[length] = c2;
	++length;
	sequence1[length] = '\0';
	sequence2[length] = '\0';
	
}

/**
 * M�thode permettant d'initialiser les adresses de d�but et de fin de l'alignement dans les s�quences
 * \param deb1 l'indice de d�part du fragment de la premi�re s�quence
 * \param deb2 l'indice de d�part du fragment de la seconde s�quence
 * \param len la longueur de l'alignement
 * \param sizeSeq1 la longueur de la plus longue s�quence
 * \rev_comp un bool�en indiquant si la premi�re s�quence est invers�e et compl�ment�e ou non
 */
void Alignment::setOffsets(int deb1, int deb2, int len, int sizeSeq1, bool rev_comp)
{
	/// Calcul des index r�els de d�but et de fin de l'alignement
	start2 = deb2 + 1;
	end2 = deb2 + len;

	start1 = deb1 + 1;
	end1 = deb1 + len;
		
}

void Alignment::adjust_rev_comp( int sizeSeq1, bool rev_comp)
{
  if(rev_comp==REV_COMP_ALIGN)
    {
      start1 = sizeSeq1 - (start1-1);
      end1 =  sizeSeq1 - (end1-1);
    }
}


//suivant format de sortie, pas la meme chose a faire
// pour seq reverse

// a la base si align sur strand reverse, le read a ete rev complemente
// donc on re rev comp ref et read pour STD_OUTPUT_FORMAT ou M8  car on veut ref reverse et read normal
// pour sam cest le contraire, on affiche toujours align sur forward ref,
// et on affiche le complem du read ( donc rien a faire pour read,ni pour cigar, bien calcule)
void Alignment::apply_rev_comp( int output_format)
{

  if(rev_comp)
    {
      if(output_format==STD_OUTPUT_FORMAT || output_format==M8_OUTPUT_FORMAT)
	{
	  char * tmp; 
	  tmp = new char [length]; 
	  for(int p=0;p<length;p++)
	    {
	      tmp[p] = sequence1[p];
	    }
	  for(int p=0;p<=length-1;p++)
	    {
	      sequence1[p] = complNTG(tmp[length-1-p]);
	    }

	  for(int p=0;p<length;p++)
	    {
	      tmp[p] = sequence2[p];
	    }
	  for(int p=0;p<=length-1;p++)
	    {
	      sequence2[p] = complNTG(tmp[length-1-p]);
	    }

	  int stemp=start1; start1=end1;  end1=stemp;
	  delete [] tmp;
	}
      else if (output_format==SAM_READY_FORMAT)
	{
	  // rien a faire 
	  ;

	}
    }
}


/**
 * M�thode permettant d'obtenir le caract�re � l'indice en param�tre dans la premi�re s�quence
 * \param indice l'indice du caract�re dans le tableau s�quence1
 * \return le caract�re � l'indice voulu dans la premi�re s�quence
 */
char Alignment::getChar1(int indice)
{
	if(indice < length && indice >= 0)
	{
		return sequence1[indice];
	}
	else return '#';
}

/**
 * M�thode permettant d'obtenir le caract�re � l'indice en param�tre dans la seconde s�quence
 * \param indice l'indice du caract�re dans le tableau s�quence2
 * \return le caract�re � l'indice voulu dans la seconde s�quence
 */
char Alignment::getChar2(int indice)
{
	if(indice < length && indice >= 0)
	{
		return sequence2[indice];
	}
	else return '#';
}
		
/**
 * M�thode d'affichage d'un Aligment
 */
void Alignment::affichage()
{
	printf("seq1 :  ");
	printf("%s ",getSeq1());
	printf("\nseq2 :  ");
	printf("%s ",getSeq2());
	printf("\n");
}








		
/**
 * M�thode de fin de  calcul du cigar
 */
void Alignment::finish_cigar()
{

  if(etat_cigar==0) return;
  sprintf(cigar+cigar_index,"%i%c",cpt_cigar,etat_cigar);
  
}
