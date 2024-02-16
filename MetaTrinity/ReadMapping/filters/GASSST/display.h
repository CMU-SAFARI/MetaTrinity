#ifndef DISPLAY_H
#define DISPLAY_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file display.h
 * \brief Module Display, responsable de l'enregistrement du r�sultat dans un fichier de sortie
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */

#include <iostream>

class Alignment;

/// Variable contenant le type de format de sortie du programme
int OUTPUT_FORMAT;


/**
 * M�thode de v�rification du format de sortie
 * \param k le type de format demand�
 */
void checkOutputFormat(int k);


/**
 * Fonction permettant d'enregistrer un alignement dans le fichier de sortie
 * \param ff le pointeur vers l'objet FILE du fichier de sortie
 * \param s1 le pointeur vers la premi�re s�quence
 * \param s2 le pointeur vers la seconde s�quence
 * \param c1 le pointeur vers le d�but du commentaire de la s�quence de la premi�re banque
 * \param c2 le pointeur vers le d�but du commentaire de la s�quence de la seconde banque
 * \param start1 la position de d�part de l'alignement pour la premi�re s�quence
 * \param start2 la position de d�part de l'alignement pour la seconde s�quence
 * \param len la longueur de l'alignement (longueur de la plus courte s�quence)
 * \param l_seq1 la longueur de la premi�re s�quence, la plus longue
 * \param eval la valeur de la E-value
 * \param rev_comp un bool�en indiquant si on a un alignement sur la s�quence normale (false) ou invers�e et compl�ment�e (true)
 */
void display_ungap(FILE *ff, char* s1, char* s2, char* c1, char* c2, int start1, int start2, int len, int l_seq1, double eval, bool rev_comp);


/**
 * Fonction permettant d'enregistrer un alignement avec gap(s) dans le fichier de sortie
 * \param ff le pointeur vers l'objet FILE du fichier de sortie
 * \param al l'objet Alignment contenant les informations sur l'alignement
 * \param c1 le pointeur vers le d�but du commentaire de la s�quence de la premi�re banque
 * \param c2 le pointeur vers le d�but du commentaire de la s�quence de la seconde banque
 * \param eval la valeur de la E-value
 */
void display_withgap(FILE *ff, Alignment* al, char* c1, char* c2, double eval);

#endif
