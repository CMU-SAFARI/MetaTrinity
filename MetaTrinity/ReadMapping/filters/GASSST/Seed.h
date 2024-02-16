#ifndef SEED_H
#define SEED_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Seed.h
 * \brief Classe Seed, d�finissant une graine
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \author Guillaume Rizk
 * \version 6.1
 * \date 19/12/2008
 */

/**
 * Le nombre de graines diff�rentes
 */
long long  NB_DIFF_SEED;

/**
 * La taille de la graine
 */
int SIZE_SEED;

/**
 * \class Seed, Une graine correspond � une petite s�quence d'ADN
 * \brief Cette class d�finit une graine, c'est � dire une courte s�quence qui permettra la recherche des alignements
 */
class Seed{
	public:
	
	/**
	 * Num�ro de la s�quence
	 */
	int num_seq;
	/**
	 * Index de la graine dans la s�quence
	 */
	int off_seq;
	/**
	 *  code s�quence des 16 carac a gauche de la graine
	 */
	unsigned int left;

	/**
	 *  code s�quence des 16 carac a droite de la graine
	 */
	unsigned int right;
	
	/**
	 * Constructeur par d�faut
	 */
	Seed();
	
	/**
	 * Constructeur de Seed
	 * \param num le num�ro de la s�quence
	 * \param offset, l'index dans s�quence
	 */
	Seed(int num, int offset, unsigned int seqleft, unsigned int seqright);
	/**
	 * Constructeur par recopie
	 * \param s un objet Seed
	 */
	Seed(const Seed& s);
	
	/**
	 * Destructeur
	 */
	~Seed();
	
	/**
	 * Op�rateur d'affectation
	 * \param s un objet Seed
	 * \return l'objet Seed affect�
	 */
	Seed& operator=(const Seed& s);
};

#endif
