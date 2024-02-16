#ifndef WITHGAP_H
#define WITHGAP_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file withgap.h
 * \brief Module Withgap, responsable de la r�alisation de l'alignement avec gap
 * \author Damien Fleury
 * \author Guillaume Rizk
 * \version 6.1
 * \date 21/01/2009
 */

#include <pthread.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/sysinfo.h>
#include <sys/time.h>

class Bank;
class Index;
class Hit;
class Alignment;
class Stat;
class Doublon;
#define LEFT 0
#define RIGHT 1

#define AL_VALID 1
#define AL_DOUBLON 2
#define AL_ERR 0

/**
 * La structure caseMatrix permet une m�morisation du chemin utilis� pour d�terminer
 * le meilleur alignement dans la matrice de l'algorithme
 */
typedef struct
{
	/**
	 * Le score correspondant � la case
	 */
	int score;
	/**
	 * L'origine de la valeur du score
	 */
	int path;
  	/**
	 * Le nombre de gap
	 */
  //	int ngaps;
	/**
	 * Le nombre de mismatch
	 */
  //	int nmis;
}caseMatrix;


/// Les diff�rentes valeurs pour l'attribut path de la structure caseMatrix
#define APP 0
#define MES 1
#define GA1 2
#define GA2 3


/**
 * Fonction de recherche du meilleur alignement selon le mod�le de Wunsch et Needleman
 * \param al l'objet Alignment qui constituera le r�sultat
 * \param seq1 le pointeur du tableau des caract�res de la premi�re banque
 * \param seq2 le pointeur du tableau des caract�res de la seconde banque
 * \param size le nombre de bases � parcourir
 * \param num_gaps le nombre maximal de gaps autoris�s dans la partie de l'alignement
 * \param side indiquant si on fait l'alignement sur le c�t� gauche ou droit de la graine
 */
//void WN_alignment(Alignment* al, char* seq1, char* seq2, int size, int num_gaps, int side);
int WN_alignment(Alignment* al, char* seq1, char* seq2, int size1, int size2, int num_gaps, int side, int num_seed, int max_mis);

int WN_big_alignment(Alignment* al, char* seq1, char* seq2, int size1, int size2, int num_gaps, int side, int num_seed, int max_mis);

int WN_very_big_alignment(Alignment* al, char* seq1, char* seq2, int size1, int size2, int num_gaps, int side, int num_seed, int max_mis);

/**
 * Fonction permettant de calculer un score d'alignement obtenu entre 2 s�quences
 * \param al l'objet Alignment qui constituera le r�sultat
 * \param seq1 le pointeur du tableau des caract�res de la premi�re banque
 * \param seq2 le pointeur du tableau des caract�res de la seconde banque
 * \param left1 le nombre de bases pr�sentes � gauche dans la premi�re s�quence
 * \param left2 le nombre de bases pr�sentes � gauche dans la seconde s�quence
 * \param right1 le nombre de bases pr�sentes � droite dans la premi�re s�quence
 * \param right2 le nombre de bases pr�sentes � droite dans la seconde s�quence
 * \param max_mis le nombre maximal de m�sappariement pour conserver l'alignement
 * \param num_gaps le nombre maximal de gaps autoris�s dans les alignements
 * \return Le nombre de m�sappariements obtenus dans l'alignement
 */
int withgap_scoreHit(Alignment* al, char* seq1, char* seq2, int left1, int left2, int right1, int right2, int max_mis, int num_gaps);


/**
 * M�thode qui � partir des positions d'une graine va effectuer tous les alignements possibles avec les s�quences des banques
 * \param ff le fichier de sortie, o� seront affich�s les alignements obtenus
 * \param BK1 le pointeur de la premi�re banque de s�quences
 * \param BK2 le pointeur de la seconde banque de s�quences
 * \param I1 le pointeur de l'index de la premi�re banque de s�quences
 * \param start num premiere query a faire
 * \param end num derniere query a faire
 * \param idpc le pourcentage de ressemblances minimal entre les 2 s�quences pour qu'un alignement soit conserv�
 * \param num_gaps le nombre maximal de gaps autoris�s dans les alignements
 * \return le nombre d'alignements obtenus
 */
int withgap_align(FILE *ff, Bank *BK1, Bank *BK2, Index * I1, int idpc, bool rev_comp, int num_gaps, char ** tabprec, char ** tabnt, char **tabprec7, Stat * St,Doublon * Doub,int start,int end);

/*
 * Max Rumpf 2020.
 *
 * Permits access of the function by our heuristics suite.
 *
 */

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#definee EXTERNC
#endif

EXTERNC int filtre_align_NT_vec(int left2, int right2, int max_mis
        ,int al,int cl,int tl, int gl, int ar,int cr,int tr, int gr
        ,unsigned int fleft2,unsigned int fright2);

EXTERNC int filtre_align_NT(char* seq1, char* seq2, int left1, int left2, int right1, int right2, int max_mis,
                    char ** tabnt
        ,unsigned int fleft1,unsigned int fright1
        ,unsigned int fleft2,unsigned int fright2);

EXTERNC int filtre_align_DC(char* seq1, char* seq2, int left2, int right2, int max_mis,char ** tabprec
        ,unsigned int fleft1,unsigned int fright1
        ,unsigned int fleft2,unsigned int fright2,int num_gaps);

EXTERNC  int filtre_align( int left2, int right2, int max_mis,
                           char ** tabprec
        ,unsigned int fleft1,unsigned int fright1
        ,unsigned int fleft2,unsigned int fright2,int num_gaps);

EXTERNC inline int pre_filtre_align( int left2,int right2, int max_mis, char ** tabprec ,unsigned int fleft1,unsigned int fright1
        ,unsigned int fleft2,unsigned int fright2);


//SIMD improvement
//EXTERNC int filtre_align_NT_vec_sse( int left, int right, int max_mis,unsigned int fleft,unsigned int fright,__m128i nt_query_left,__m128i nt_query_right,
//                                     __m128i vec_5,__m128i vec_3,
//                                     __m128i vec_0f,__m128i vec_h1,__m128i vec_nt
//        /*	,int * a1,int * a2,int * c1,int * c2,int * t1,int * t2,int * g1,int * g2,int * minerr*/);

#undef EXTERNC


#endif
