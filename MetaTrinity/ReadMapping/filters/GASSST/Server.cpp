/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Server.cpp
 * \brief Class Server, responsable du partage des t�ches entre les threads
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */

//modification , parallelisatio nsur Nquery au lieu N graines
#include "Server.h"
#include "misc.h"

/**
 * Le constructeur par d�faut de Server
 */
Server::Server() : last_seed(0), size_task(0), index(0)
{
	pthread_rwlock_init(&verrou_serveur,NULL);
}

/**
 * Le destructeur par d�faut de Server
 */
Server::~Server()
{
}

/**
 * L'op�rateur d'affectation de Server
 * \param s un objet Server
 * \return l'objet server
 */
Server& Server::operator=(const Server& s)
{
	if(this != &s)
	{
		last_seed = s.last_seed;
		size_task = s.size_task;
		index = s.index;
		nb_part = s.nb_part;
	}
	return *this;
}

/**
 * Le constructeur par recopie de Server
 * \param s un objet Server
 */
Server::Server(const Server& s)
{
	*this = s;
}

/**
 * Constructeur de Server
 * `param nbSeeds un entier correspondant au nombre total de graines � traiter
 * \param nb_partitions un entier correspondant au nombre de partitions � faire avec
 * l'ensemble des graines
 */
Server::Server(int nbSeeds, int nb_partitions) : last_seed(nbSeeds-1), size_task(nb_partitions), index(0), nb_part(nbSeeds / nb_partitions)
{
	pthread_rwlock_init(&verrou_serveur,NULL);
}

/**
 * M�thode permettant d'obtenir un sous-ensemble de graines � traiter
 * \param t un pointeur vers une structure task, correspondant � une partition de l'ensemble des graines
 * \return 1 si il reste des graines � traiter, on a dans ce cas affect� en cons�quence les champs
 * de la structure task point�e par t, 0 si toutes les graines ont �t� attribu�es, les champs de la
 * structure point�e par t sont mis � -1
 */
int Server::give_task(task* t)
{
	int id;
	int tai=0;

	/// Section critique
	pthread_rwlock_wrlock(&verrou_serveur);
	id = index;
	if(index <= last_seed)
	{
	  tai = min (size_task,last_seed+1-index);
	  index+=tai;
	}
	pthread_rwlock_unlock(&verrou_serveur);
	/// Fin de la section critique
	
	if(id <= last_seed)
	{
		t->start = id;
		t->end = id + tai -1;
		return 1;
	}
	else
	{
		t->start = -1;
		t->end = -1;
		return 0;
	}

}

/**
 * M�thode permettant de r�initialiser l'index du server
 */
void Server::reset()
{
	index = 0;
}
