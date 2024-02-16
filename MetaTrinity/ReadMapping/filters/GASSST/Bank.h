#ifndef BANK_H
#define BANK_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Bank.h
 * \brief Classe Bank, responsable de la r�cup�ration des donn�es dans les banques de s�quences
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */

#include <iostream>

/// Nombre maximal de partitions
#define MAX_PART 1024
/// Taille d'une ligne d'un fichier au format fasta
#define SIZE_LINE 1024
/// TAille du nom des banques de s�quences
#define TAILLE_NOM 1024

/**
 * \class Bank, Une banque de s�quences d'ADN
 * \brief Cette classe d�finit une banque de s�quences d'ADN, ainsi que les m�thodes permettant d'en r�cup�rer les informations
 */
class Bank{

	public:
	/**
	 * Nom de la banque
	 */
	char fileBank[TAILLE_NOM];
	/**
	 * Nombre de s�quences total
	 */
	int  nb_tot_seq;
	/**
	 * Nombre de r�sidus total
	 */
	long long  int  nb_tot_res;
	/**
	 * Nombre de partitions
	 */
	int  nb_part;
	/**
	 * Num�ro de la partition en cours de traitement
	 */
	int  num_part;
	/**
	 * Num�ro de la partition suivante
	 */
	int  next_part;
	/**
	 * Adresses de d�part des partitions dans le fichier
	 */
	long* start_offset;
	/**
	 * Adresses de fin des partitions dans le fichier
	 */
	long* stop_offset;
	/**
	 * Nombre de s�quences dans chaque partition
	 */
	int* nb_seq;
	/**
	 * Image m�moire de la partie de banque en cours de traitement
	 */
	char* data;
	/**
	 * Tableau des index des sequences , en int donc taille max partition = 2G
	 */
	int* seq;
	/**
	 * Tableau des index des commentaires
	 */
	int* com;
	/**
	 * Tableau des tailles des s�quences
	 */
	long long* size;
	
	/**
	 * Tableau des positions globales deu d�but des s�quences
	 */
	long long* pos_seq;
	/**
	 * Longueur maximale des s�quences de la partition index�e de la banque
	 */
	int tailleMaxSeq;

	int tSeq;


	/**
	 * Constructeur de banque de s�quences
	 * \param fname, un pointeur de caract�res contenant le nom de la banque de s�quences
	 * \param size_max, la taille maximale de la banque de s�quences
	 * \param FILE, fichier de sortie, pour mettre header avec noms des contig
	 * \param bankref indique qu on lit la banque de reference
	 */
	Bank(char *fname,long long size_max, FILE *ff, char bankref);
	
	/**
	 * Constructeur de banque par recopie
	 * \param bk, une banque de s�quences
	 */
	Bank(const Bank& bk);
	
	/**
	 * Op�rateur d'affectation
	 * \param bk une banque de s�quences
	 * \return l'objet Bank affect�
	 */
	Bank& operator=(const Bank& bk);
	
	/**
	 * Destructeur de Bank
	 */
	~Bank();
	
	/**
	* M�thode d'affichage d'une banque de s�quences
	*/
	void writeInfoBank();
	
	/**
	* M�thode de r�initialisation des partitions d'une banque de s�quences
	*/
	void resetBank();
	
	/**
	 * M�thode permettant d'indexer une partition de la banque
	 * \param lx, un bool�en qui indique si on utilise le filtre Low Complexity
	 * \return 1 si une partition a �t� index�e, 0 si toutes les partitions sont d�j� �t� index�es
	 */
	int readBank(bool lx);
	/**
	 * M�thode permettant de faire le "reverse complement" des s�quences de la partition courante
	 * de la banque
	 * Toutes les s�quences de la partition courante de la banque seront invers�es et leurs bases
	 * seront remplac�es par les bases compl�mentaires
	 */
	void reverseComplement();
};

#endif
