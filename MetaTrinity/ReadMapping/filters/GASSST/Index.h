#ifndef INDEX_H
#define INDEX_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Index.h
 * \brief Classe Index, responsable de l'indexage des banques de s�quences
 * \author Dominique Lavenier
 * \author Guillaume Rizk
 * \author Damien Fleury
 * \date 28/08/2008
 */

#include "Seed.h"
#include "Pool.h"


class Bank;

/**
 * \class Index, Un index des graines d'une banque
 * \brief Cette classe permet de r�pertorier les graines d�tect�es dans une banque de s�quences d'ADN
 */
class Index{
	
	public:
	/**
	 * Tableau contenant les index des graines dans le tableau seed
	 */
	int* offset_seed;
	/**
	 * Tableau contenant le nombre d'occurences de chaque graine
	 */
	int* nb_seed;
	/**
	 * Tableau des graines
	 */
	Seed* seed;

	/**
	 * pool memoire pour graine > 14
	 */
	Pool* storage;


	/**
	 * Table de hachage pour les graines de taille >14
	 */	
	cell ** table_hachage;
	/**
	 * Constructeur d'Index par d�faut
	 */
	Index();
	
	/**
	 * Destructeur d'Index
	 */
	~Index();
	
	/**
	 * Constructeur d'Index par recopie
	 * \param i un objet Index
	 */
	Index(const Index& i);
	
	/**
	 * Op�rateur d'affectation de Hit
	 * \param i un objet Index
	 * \return l'objet Index affect�
	 */
	Index& operator=(const Index& i);
	
	/**
	 * M�thode permettant d'indexer une banque
	 * \param B le pointeur de la banque � indexer
	 */
	void indexBank(Bank *B,int index_stride);
	/**
	 * fonction de hachage, renvoit clef
	 * \param  code  la grande graine 
	 * \retrun  clef 
	 */
	int hashCode(long long code, int max);

	void get_seed_info_through_hashTable(long long seed, int * nb_occur, int * offset_seed );


	int dec_haschcode;
	long long mask1;
	long long mask2;
	long long mask3;

	void free_hashtable();


	//for debug purposes
	void printIndex();
	void printstat();


};

#endif
