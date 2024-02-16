#ifndef ALIGNMENT_H
#define ALIGNMENT_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Alignment.h
 * \brief Class Alignment, d�finissant un alignement entre deux s�quences
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */


/**
 * \class Alignment, Un alignement entre deux s�quences
 * \brief Un alignement contient toutes les informations identifiant l'alignement et permettant de l'enregistrer
 * Cette classe est principalement utilis�e pour le cas des alignements avec gap
 */


class Alignment
{
private:
	/**
	 * La premi�re s�quence de l'alignement, c'est-�-dire celle de la base de donn�es
	 * La cha�ne est termin�e par un \0
	 */
	char* sequence1;
	
	/**
	 * La seconde s�quence de l'alignement, c'est-�-dire celle de la base de requ�tes
	 * La cha�ne est termin�e par un \0
	 */
	char* sequence2;

	/**
	 * La longueur de l'alignement en nombre de caract�res
	 */
	int length;
	
	/**
	 * Le nombre de m�sappariements dans l'alignement
	 */
	int nb_mismatches;
	
	/**
	 * Le nombre de gaps dans l'alignement
	 */
	int nb_gaps;
	
	/**
	 * L'indice de d�part de l'alignement dans la premi�re s�quence
	 */
	int start1;
	
	/**
	 * L'indice de d�part de l'alignement dans la seconde s�quence
	 */
	int start2;
	
	/**
	 * L'indice de fin de l'alignement dans la premi�re s�quence
	 */
	int end1;
	
	/**
	 * L'indice de fin de l'alignement dans la seconde s�quence
	 */
	int end2;
	

	//2variables utiles pour le calcul du cigar
	/**
	 * operation cigar en cours : 'M'  'I' ou 'D'
	 */
	char etat_cigar;

	/**
	 *  longueur de l'operation en cours du cigar
	 */
	unsigned int cpt_cigar; 

	/**
	 *  nombre de char deja imprimes du cigar
	 */
	unsigned int cigar_index;

public:
	
	/**
	 * indique si alignement sur strand -, dans ce cas il faudra reverse la sortie
	 */
	char rev_comp;

	//score evalue de l alignement
	double e_value ;

	//index des commentaires des sequences
	char * j1;

	char * j2;

	//le numero des sequences,relatif a la partition
	unsigned int  n1,n2;

	//description de l'align en cigar, (pour format sam)
	char * cigar;


	void finish_cigar();


	/**
	 * Le consructeur par d�faut d'Alignment
	 */
	Alignment();
	
	/**
	 * Un consructeur d'Alignment
	 * \param size la taille de l'alignement
	 */
	Alignment(int size);
	
	/**
	 * La fonction d'initialisation d'un Alignment
	 * \param size la taille de l'alignement
	 */
	void init();
	
	/**
	 * La fonction d'initialisation d'un Alignment
	 * \param rever si align sur strand -
	 */
	void init(char rever);
	/**
	 * Le consructeur par recopie d'Alignment
	 * \param al un objet Alignment
	 */
	Alignment(const Alignment& al);
	
	/**
	 * L'op�rateur d'affectation d'Alignment
	 * \param al un objet Alignment
	 * \return l'objet Alignment affect�
	 */
	Alignment& operator=(const Alignment& al);
	
	/**
	 * Le destructeur par d�faut d'Alignment
	 */
	virtual ~Alignment();
	
	/**
	 * M�thode d'obtention de la longueur de l'alignement
	 * \return la longueur de l'objet Alignment appelant
	 */
	inline int getLength();
	
	/**
	 * M�thode d'obtention du nombre de m�sappariements
	 * \return le nombre de m�sappariements de l'objet Alignment appelant
	 */
	inline int getMis();
	
	/**
	 * M�thode d'obtention du nombre de gaps
	 * \return le nombre de gaps de l'objet Alignment appelant
	 */
	inline int getGaps();
	
	/**
	 * M�thode d'incr�mentation du nombre de m�sappariements
	 */
	inline void addMis();
	
	/**
	 * M�thode d'incr�mentation du nombre de gaps
	 */
	inline void addGap();

	/**
	 * M�thode permettant d'initialiser les adresses de d�but et de fin de l'alignement dans les s�quences
	 * \param deb1 l'indice de d�part du fragment de la premi�re s�quence
	 * \param deb2 l'indice de d�part du fragment de la seconde s�quence
	 * \param len la longueur de l'alignement
	 * \param sizeSeq1 la longueur de la plus longue s�quence
	 * \rev_comp un bool�en indiquant si la premi�re s�quence est invers�e et compl�ment�e ou non
	 */
	void setOffsets(int deb1, int deb2, int len, int sizeSeq1, bool rev_comp);
	
	//complement the start position if rev _comp
	void adjust_rev_comp( int sizeSeq1, bool rev_comp);



	//complement the alignment for correct output
	void apply_rev_comp(int output_format);
	/**
	 * Fonction utilis�e pour incr�menter la valeur de l'indice de d�part dans la premi�re s�quence
	 * On peut ainsi tenir compte des gaps utilis�s
	 */
	inline void incAlign();
	
	/**
	 * Fonction utilis�e pour d�cr�menter la valeur de l'indice de d�part dans la premi�re s�quence
	 * On peut ainsi tenir compte des gaps utilis�s
	 */
	inline void decAlign();
	
	inline void incEnd1();
	inline void decStart1();

	inline void decEnd1();
	inline void incStart1();
	/**
	 * M�thode d'acc�s � l'attribut start1
	 * \return la valeur de start1
	 */
	inline int getStart1();
	
	/**
	 * M�thode d'acc�s � l'attribut end1
	 * \return la valeur de end1
	 */
	inline int getEnd1();
	
	/**
	 * M�thode d'acc�s � l'attribut start2
	 * \return la valeur de start2
	 */
	inline int getStart2();
	
	/**
	 * M�thode d'acc�s � l'attribut end2
	 * \return la valeur de end2
	 */
	inline int getEnd2();
	
	/**
	 * M�thode d'ajout d'un couple d'un caract�re aux s�quences
	 * \param c1 le caract�re � ajouter � la premi�re s�quence
	 * \param c2 le caract�re � ajouter � la seconde s�quence
	 */
	 void addPair(char c1, char c2);

	/**
	 * M�thode permettant de r�cup�rer la premi�re s�quence
	 * \return la premi�re s�quence de l'Alignment
	 */
	inline char* getSeq1();
	
	/**
	 * M�thode permettant de r�cup�rer la seconde s�quence
	 * \return la seconde s�quence de l'Alignment
	 */
	inline char* getSeq2();
	 
	/**
	 * M�thode permettant de r�cup�rer le string cigar
	 * \return le string cigar
	 */
	inline char* getCigar();
	/**
	 * M�thode permettant d'obtenir le caract�re � l'indice en param�tre dans la premi�re s�quence
	 * \param indice l'indice du caract�re dans le tableau s�quence1
	 * \return le caract�re � l'indice voulu dans la premi�re s�quence
	 */
	char getChar1(int indice);
	
	/**
	 * M�thode permettant d'obtenir le caract�re � l'indice en param�tre dans la seconde s�quence
	 * \param indice l'indice du caract�re dans le tableau s�quence2
	 * \return le caract�re � l'indice voulu dans la seconde s�quence
	 */
	char getChar2(int indice);	
	
	/**
	 * M�thode d'affichage d'un Aligment
	 */
	void affichage();
};


/**
 * M�thode d'obtention de la longueur de l'alignement
 * \return la longueur de l'objet Alignment apelant
 */
int Alignment::getLength()
{
	return length;
}

/**
 * M�thode d'obtention du nombre de m�sappariements
 * \return le nombre de m�sappariements de l'objet Alignment appelant
 */
int Alignment::getMis()
{
	return nb_mismatches;
}

/**
 * M�thode d'obtention du nombre de gaps
 * \return le nombre de gaps de l'objet Alignment appelant
 */
int Alignment::getGaps()
{
	return nb_gaps;
}

/**
 * M�thode d'incr�mentation du nombre de m�sappariements
 */
void Alignment::addMis()
{
	++nb_mismatches;
}

/**
 * M�thode d'incr�mentation du nombre de gaps
 */
void Alignment::addGap()
{
	++nb_gaps;
}

/**
 * Fonction utilis�e pour incr�menter la valeur de l'indice de d�part dans la premi�re s�quence
 * On peut ainsi tenir compte des gaps utilis�s
*/
void Alignment::incAlign()
{
	++start1;
}

   void Alignment::incEnd1()
{
	++end1;

}
   void Alignment::decStart1()
{
	--start1;

}

 void Alignment::decEnd1()
{
	--end1;

}
   void Alignment::incStart1()
{
	++start1;

}
/**
 * Fonction utilis�e pour d�cr�menter la valeur de l'indice de d�part dans la premi�re s�quence
 * On peut ainsi tenir compte des gaps utilis�s
 */
void Alignment::decAlign()
{
	--start1;
}

/**
 * M�thode d'acc�s � l'attribut start1
 * \return la valeur de start1
 */
inline int Alignment::getStart1()
{
	return start1;
}

/**
 * M�thode d'acc�s � l'attribut end1
 * \return la valeur de end1
 */
inline int Alignment::getEnd1()
{
	return end1;
}

/**
 * M�thode d'acc�s � l'attribut start2
 * \return la valeur de start2
 */
inline int Alignment::getStart2()
{
	return start2;
}
	
/**
 * M�thode d'acc�s � l'attribut end2
* \return la valeur de end2
 */
inline int Alignment::getEnd2()
{
	return end2;
}

/**
 * M�thode permettant de r�cup�rer la premi�re s�quence
 * \return la premi�re s�quence de l'Alignment
 */
inline char* Alignment::getSeq1()
{
	return sequence1;
}

/**
 * M�thode permettant de r�cup�rer la seconde s�quence
 * \return la seconde s�quence de l'Alignment
 */
inline char* Alignment::getSeq2()
{
	return sequence2;
}


/**
 * M�thode permettant de r�cup�rer le string cigar
 * \return le string cigar
 */
inline char* Alignment::getCigar()
{
  return cigar;
}


#endif
