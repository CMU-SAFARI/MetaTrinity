#ifndef resu_H
#define resu_H

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file resu.h
 * \author Guillaume Rizk
 * \date 01/02/2011
 */


/**
 * \class resu, Une graine correspond � une petite s�quence d'ADN
 * \brief Cette class d�finit une graine, c'est � dire une courte s�quence qui permettra la recherche des alignements
 */
class resu{
 public:
  
  
  unsigned int score_al;
  
  unsigned int index;




	
	/**
	 * Constructeur par d�faut
	 */
	resu();
	
	/**
	 * Constructeur de resu
	 * \param num le num�ro de la s�quence
	 * \param offset, l'index dans s�quence
	 */
	resu(  unsigned int score_in,
	       unsigned int index_in);

	
	resu(const resu &);
	~resu(){};
	resu &operator=(const resu &rhs);
	int operator==(const resu &rhs) const;
	int operator<(const resu &rhs) const;


};

#endif
