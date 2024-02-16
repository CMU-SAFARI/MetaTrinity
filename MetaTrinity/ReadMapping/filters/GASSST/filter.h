#ifndef FILTER_H
#define FILTER_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file filter.h
 * \brief Module Filter, d�finit le filtre Low Complexity qui d�tecte les zones non pertinentes des s�quences afin de ne pas les indexer
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */

/// Fen�tre du filtre Low complexity
#define WS 12


using namespace std;


/**
 * M�thode permettant d'appliquer le filtre Low Complexity
 * \param data le tableau des donn�es de la banque de s�quences
 * \param id l'index du d�but de la s�quence
 * \param lenseq la longueur de la s�quence
 * \return 1 si l'op�ration s'est d�roul�e correctement
 */
int filterLowComplexity (char* data, int lenseq);

#endif

