#ifndef HIT_H
#define HIT_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Hit.h
 * \brief Classe Hit, d�finissant un alignement entre 2 s�quences
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */

class Bank;

/**
 * \class Hit, Une position d'une graine dans une s�quence
 * \brief Cette classe correspond � une position d'une graine dans une s�quence d'une banque
 */
class Hit{
	 public:

	/**
	 * Index de la position de la graine
	 */
	int offhit;
 	/**
	 * Index de la s�quence
	 */
	int offseq;
	/**
	 * Taille de la s�quence
	 */
	int sizeseq;
	/** 
	 * Num�ro de la s�quence
	 */
	int numseq;
	
	//	unsigned char seqleft;
	//	unsigned char seqright;
	/**
	 * Constructeur par d�faut
	 */
	Hit(); 
	
	/**
	 * Constructeur de Hit
	 * \param BK un pointeur vers la banque o� est cr�� l'alignement
	 * \param num_sequence le num�ro de la s�quence o� est situ� l'alignement
	 * \param offset_sequence la position de l'alignement dans la s�quence
	 */
	Hit(Bank *BK, int num_sequence, int offset_sequence);
	
	/**
	 * Constructeur de Hit par recopie
	 * \param h un objet Hit
	 */
	Hit(const Hit& h);
	
	/**
	 * Destructeur de Hit
	 */
	~Hit();
	
	/**
	 * Op�rateur d'affectation de Hit
	 * \param h un objet Hit
	 * \return l'objet Hit affect�
	 */
	Hit& operator=(const Hit& h);

	//operateur ordre utilis� pour tri
	int operator<(const Hit &h) const;

};

#endif
