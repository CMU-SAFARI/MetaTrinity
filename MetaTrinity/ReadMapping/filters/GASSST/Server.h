#ifndef SERVER_H
#define SERVER_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Server.h
 * \brief Class Server, responsable du partage des t�ches entre les threads
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */


#include <pthread.h>

/**
 * Structure d�finissant une portion de graines � traiter 
 */
typedef struct
{
	/**
	 * Le num�ro de la premi�re graine de la portion
	 */
	int start;
	
	/**
	 * Le num�ro de la derni�re graine de la portion
	 */
	int end;
} task;

/**
 * \class Server, Un objet permettant de r�partir le travail de recherche d'alignements entre diff�rents threads
 * \brief Un objet Server g�re l'ensemble des graines � traiter lors de la recherche, les threads vont 
 * demander une t�che � l'objet Server qui va leur attribuer une partie de l'ensemble des graines.
 * Le Server va donc r�partir le travail le plus �quitablement possible entre les diff�rents threads 
 */
class Server
{
private:
	/**
	 * Le num�ro de la derni�re graine � traiter
	 */
	int last_seed;
	
	/**
	 * La taille d'un ensemble de graines attribu� lors d'une demande par un thread
	 */
	int size_task;
	
	/**
	 * Un indice permettant de r�pertorier les graines qui n'ont pas encore �t� attribu�es � un thread
	 */
	int index;
	
	/**
	 * Le nombre de sous-ensembles de graines � attribuer aux threads
	 */
	int nb_part;
	
	/**
	 * Un verrou emp�chant les acc�s concurentiels aux variables du Server par les threads
	 */
	pthread_rwlock_t verrou_serveur;

public:
	/**
	 * Le constructeur par d�faut de Server
	 */
	Server();
	
	/**
	 * Le destructeur par d�faut de Server
	 */
	virtual ~Server();
	
	/**
	 * Le constructeur par recopie de Server
	 * \param s un objet Server
	 */
	Server(const Server& s);
	
	/**
	 * L'op�rateur d'affectation de Server
	 * \param s un objet Server
	 * \return l'objet server
	 */
	Server& operator=(const Server& s);
	
	/**
	 * Constructeur de Server
	 * `param nbSeeds un entier correspondant au nombre total de graines � traiter
	 * \param nb_partitions un entier correspondant au nombre de partitions � faire avec l'ensemble des graines
	 */
	Server(int nbSeeds, int nb_partitions);
	
	/**
	 * M�thode permettant d'obtenir un sous-ensemble de graines � traiter
	 * \param t un pointeur vers une structure task, correspondant � une partition de l'ensemble des graines
	 * \return 1 si il reste des graines � traiter, on a dans ce cas affect� en cons�quence les champs de la structure task point�e par t, 0 si toutes les graines ont �t� attribu�es, les champs de la structure point�e par t sont mis � -1
	 */
	int give_task(task* t);
	
	/**
	 * M�thode permettant de r�initialiser l'index du server
	 */
	void reset();
};

#endif
