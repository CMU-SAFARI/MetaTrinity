#ifndef MISC_H
#define MISC_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file misc.h
 * \brief Module Misc, contenant des fonctions utilitaires, et de gestion des options et des erreurs
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */


#include <iostream>
#include <stdlib.h>
#include <string.h>


/// Variables globales utiles � la gestion des options
/// Elles permettent la m�morisation des options sp�cifi�es
char KEY0[64];
char KEY1[64];
char PARAM0[64][64];
char PARAM1[64][64];
extern int  NB_KEY1;
extern int  NB_KEY0;
char PROG_NAME[1024];
extern int NBTHREADS;
extern int MAXHITS;
extern int MAXPOS;
extern int MAXPOS_AMONT;
extern int BESTAL;
extern int SLEVEL;
extern int BITSTAT;
extern int NUMGAPS_AUTO;


/*
 * Fontion minimum de deux valeurs
 * \param a un entier
 * \param b un entier
 * \return la valeur minimale entre a et b
 */
inline int min(int a, int b)
{
	return a < b ? a : b;
}

/*
 * Fontion maximum de deux valeurs
 * \param a un entier
 * \param b un entier
 * \return la valeur maximale entre a et b
 */
inline int max(int a, int b)
{
	return a > b ? a : b; 
}

/**
 * Fonction d'impession d'un message d'erreur et de fermeture du programme
 * \param msg la cha�ne du message d'erreur
 */
void ExitError(char *msg);


/**
 * Fonction d'affichage d'un erreur dans les param�tres du programme et de fermeture du programme
 */
void SyntaxError();


/**
 * M�thode de g�n�ration des diff�rentes options du programme
 * \param key le caract�re permettant d'identifier l'option
 * \param strg_in le nom de l'option
 * \param opt indique si l'option est obligatoire, 1 si oui, 0 sinon
 */
void option(char key, char *strg_in, int opt);


/**
 * M�thode permettant de v�rifier les options sp�cifi�es lors de l'ex�cution du programme et d'en r�cup�rer les param�tres
 * \param key le caract�re permettant d'identifier l'option
 * \param strg_out la valeur utilis�e pour l'option
 * \param argc le nombre de param�tres utilis�s lors de l'appel au programme
 * \param argv un pointeur vers le tableau des param�tres du programme
 * \param opt indique si l'option est obligatoire, 1 si oui, 0 sinon
 * \return 1 si l'option a �t� sp�cifi�e, 0 si l'option est facultative et non utilis�e
 */
int getoption(char key, char *strg_out, int argc, char *argv[], int opt);


/**
 * M�thode de v�rification des options
 * \param argc le nombre de param�tres utilis�s lors de l'appel au programme
 * \param argv un pointeur vers le tableau des param�tres du programme
 */
void checkoption(int argc, char *argv[]);

#endif
