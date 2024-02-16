/* This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use, 
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info". 
 */

/**
 * Logiciel Gassst (Global Alignment Short Sequence Search Tool)
 * \file Bank.cpp
 * \brief Classe Bank, responsable de la r�cup�ration des donn�es dans les banques de s�quences
 * \author Dominique Lavenier
 * \author Damien Fleury
 * \version 5.2
 * \date 28/08/2008
 */

#include "Bank.h"

#include <fstream>
#include "code.h"
#include "filter.h"
#include "display.h"
#include "constants.h"


/**
 * Constructeur de banque de s�quences
 * \param fname, un pointeur de caract�res contenant le nom de la banque de s�quences
 * \param size_max, la taille maximale de la banque de s�quences
 * \param FILE, fichier de sortie, pour mettre header avec noms des contig
 * \param bankref indique qu on lit la banque de reference
 */
//probleme si taille sequence > taille max partition
// on met au moins un contig dans chaque partition, meme si trop gros pour size_max
Bank::Bank(char *fname,long long size_max, FILE *ff, char bankref)
{
	FILE *fbank;
	long offset_line, last_offset_line;
	char line[SIZE_LINE];
	char buf[SIZE_LINE];

	int  nb_res;

  	/// Recopie du nom de la banque
	strcpy(fileBank,fname);
 
  	/// Ouverture du fichier de la banque au format FASTA
	if ((fbank=fopen(fileBank,"r"))==NULL)
	{
		fprintf (stderr,"cannot open %s\n",fileBank);
		exit (0);
	}

  	/// Initialisation
	nb_tot_seq = 0;
	nb_tot_res = 0;
	nb_part    = 0;
	stop_offset = new long[MAX_PART];
	start_offset = new long[MAX_PART];
	nb_seq = new int[MAX_PART];
	tSeq=0;

	data = new char[0];
  	seq = new int[0];
 	com = new int[0];
  	size = new long long[0];
	pos_seq = new long long[0];

  	tailleMaxSeq = 0;

  	/// Comptage du nombre de s�quences / r�sidus
  	/// D�termination du nombre de partitions - m�morisation des offsets

	last_offset_line = 0;
	nb_res = 0;
	
	start_offset[0]=0;
	/// On r�cup�re chaque ligne du fichier
	while (fgets(line,SIZE_LINE,fbank)!=NULL)
	{
     	/// On m�morise l'adresse dans le fichier
		offset_line = ftello64(fbank);
		if (line[0]=='>')
		{

		  // sortie header pour SAM
		  if(OUTPUT_FORMAT==SAM_READY_FORMAT && nb_tot_seq!=0 && bankref)
		    fprintf (ff,"LN:%i\n",nb_res); // taille sequence precedente 

		  // sortie header pour SAM
		  if(OUTPUT_FORMAT==SAM_READY_FORMAT && bankref)
		    {
		      sscanf(&line[1],"%s",buf); // pour n'avoir que debut jusqua premier espace
		      fprintf (ff,"@SQ\tSN:%s\t",buf);  //sortie taille seq ensuite
		    }


	  		/// Cas d'une ligne de description de la s�quence
			++nb_tot_seq;
			nb_tot_res = nb_tot_res + nb_res;
			if(tSeq<nb_res)  tSeq = nb_res;
			nb_res = 0;
	  		/// On v�rifie si la partition courante est inf�rieure � size_max
			if (last_offset_line - start_offset[nb_part] < size_max)
			{
				stop_offset[nb_part]=last_offset_line;
				++nb_seq[nb_part];
			}
			else if (start_offset[nb_part] != stop_offset[nb_part]) // partition assez grande, et non vide : on en cree une nouvelle
			  {
			  	      		/// On fait une nouvelle partition
			    ++nb_part;
			    nb_seq[nb_part] = 1;
			    start_offset[nb_part] = stop_offset[nb_part-1];
			    stop_offset[nb_part] = stop_offset[nb_part-1];
			    /// On v�rifie s'il n'y a pas trop de partitions
			    if (nb_part>=MAX_PART)
				{
				  fprintf (stderr,"number of partitions >= %d\n",MAX_PART);
				  exit (0);
				}  
			  }
			else // partition trop grande, mais vide sinon: on met quand meme le contig dedans
			{
			  stop_offset[nb_part]=last_offset_line;
			  ++nb_seq[nb_part];
			}
		}
		else
		{
	  		/// Cas d'une ligne de nucl�otides
			nb_res = nb_res + offset_line - last_offset_line - 1;
		}
		last_offset_line = offset_line;
	}

	// sortie header pour SAM
	if(OUTPUT_FORMAT==SAM_READY_FORMAT && nb_tot_seq!=0 && bankref)
	  fprintf (ff,"LN:%i\n",nb_res); // taille sequence precedente 

	if(tSeq<nb_res)  tSeq = nb_res;
  	/// On finit de param�trer la derni�re partition
	nb_tot_res = nb_tot_res + nb_res;
	stop_offset[nb_part] = ftello64(fbank);
	++nb_seq[nb_part];
	++nb_part;
  	/// On pr�cise la prochaine partition � traiter (la premi�re)
	num_part = 0;
	fclose(fbank);
}



/**
 * Op�rateur d'affectation
 * \param bk une banque de s�quences
 * \return l'objet Bank affect�
 */
Bank& Bank::operator=(const Bank& bk)
{
	if(this!=&bk)
	{
		strcpy(fileBank,bk.fileBank);
		nb_tot_seq = bk.nb_tot_seq;
		nb_tot_res = bk.nb_tot_res;
		nb_part = bk.nb_part;
		num_part = bk.num_part;
		next_part = bk.next_part;
	
		start_offset = bk.start_offset;
		stop_offset = bk.stop_offset;
		nb_seq = bk.nb_seq;
		data = bk.data;
		seq = bk.seq;
		com = bk.com;
		size = bk.size;
		pos_seq = bk.pos_seq;
		tailleMaxSeq = bk.tailleMaxSeq;
	}
	return *this;
}

/**
 * Constructeur de banque par recopie
 * \param bk une banque de s�quences
 */
Bank::Bank(const Bank& bk)
{
	*this = bk;
}
	
/**
 * Destructeur de Bank
 */
Bank::~Bank()
{
	delete [] start_offset;
	delete [] stop_offset;
	delete [] nb_seq;
	delete [] data;
	delete [] seq;
	delete [] com;
	delete [] size;
	delete [] pos_seq;

}


/**
 * M�thode d'affichage d'une banque de s�quences
 */
void Bank::writeInfoBank()
{
  printf ("%s:  %d sequences / %lli residues\n",fileBank,nb_tot_seq,nb_tot_res);
  // for (int i=0; i< nb_part; i++)
  // {printf("partition %i debut %li  fin %li \n",i,start_offset[i],stop_offset[i]);}
}


/**
 * M�thode de r�initialisation des partitions d'une banque de s�quences
 */
void Bank::resetBank()
{
  num_part = 0;
  next_part = 0;
}




/**
 * M�thode permettant d'indexer une partition de la banque
 * \param lx, un bool�en qui indique si on utilise le filtre Low Complexity
 * \return 1 si une partition a �t� index�e, 0 si toutes les partitions sont d�j� �t� index�es
 */
int Bank::readBank(bool lx)
{
  /// Descripteur du fichier de la banque
  FILE *fbank;

  // int i,j,l;
  int nb_sequences;
  long long  size_part,i,j,l;
  num_part = next_part;
  /// Toutes les partitions ont �t� analys�es
  if (num_part >= nb_part) return 0;
  /// On ouvre le fichier de la banque de s�quences
  if ((fbank=fopen(fileBank,"r"))==NULL) 
  {
  	fprintf (stderr,"cannot open %s\n",fileBank);
  	exit(0);
  }

  size_part = stop_offset[num_part]-start_offset[num_part];

  fseeko64(fbank,start_offset[num_part],SEEK_SET);

  /// Redimensionnement des vecteurs en fonction de la taille de la partition � traiter
  delete [] data;
  delete [] seq;
  delete [] com;
  delete [] size;
  delete [] pos_seq;

  data = new char[size_part];
  seq = new int[nb_seq[num_part]+8];
  com = new int[nb_seq[num_part]+8];
  size = new long long[nb_seq[num_part]+8];
  pos_seq = new long long[nb_seq[num_part]+10];
  pos_seq[0]=0;
 
  /// On copie les caract�res de la partition dans le tableau data

  fread(data,1,size_part,fbank);
  fclose(fbank);


  /// M�morisation des s�quences et des commentaires
  nb_sequences=0;
  i=0;

  /// On parcourt la partition
  while (i<size_part) 
  {
    if (data[i]=='>')
	{
	  /// On lit le commentaire
	  com[nb_sequences]=i;
	
	  while (data[i]!='\n') ++i;
	  //si fichier dos il faut enlever aussi le carriage return
	  if (data[i-1]=='\r') data[i-1]='\0';
	  data[i]='\0';
	  ++i;
	}
    else
	{
	  /// On lit la s�quence de nucl�otides
	  seq[nb_sequences] = i;
	  j=i; l=0;
	  /// On boucle jusqu'� la fin de la partition ou de la s�quence
	  while ((i<size_part)&&(data[i]!='>'))
	  {
	      if (((data[i]>='a')&&(data[i]<'z'))||((data[i]>='A')&&(data[i]<'Z')))
		  /// Recopie des lettres en majuscules
	      {
	      	data[j]=majuscule(data[i]);
	      	++j;
	      	++l;
	      }
	      ++i;
	  }
	  /// Filtre low Complexity
	  if (lx) filterLowComplexity(data+seq[nb_sequences],l);
	  /// On enregistre la taille de la s�quence
	  size[nb_sequences] = l;
	  pos_seq[nb_sequences+1] =  pos_seq[nb_sequences] + l;

	  ++nb_sequences;
	  
	  if(tailleMaxSeq<l)
	  {
	      tailleMaxSeq = l;
	  }
	  
	}
  }

  /// On enregistre le nombre de s�quences
  nb_seq[num_part]=nb_sequences;

  /// On passe � la partition suivante
  ++next_part;
  
  return 1;
}

/**
 * M�thode permettant de faire le "reverse complement" des s�quences de la partition courante
 * de la banque
 * Toutes les s�quences de la partition courante de la banque sont invers�es et leurs bases
 * sont remplac�es par les bases compl�mentaires
 */
void Bank::reverseComplement()
{
	int n,i,j,k,l;
	char* tmp;
	/// On r�cup�re le nombre de s�quences de la partition
  	n = nb_seq[num_part];
  	
  	/// On parcourt les s�quences
  	for(j=0;j<n;++j)
  	{
  		/// On recopie la s�quence dans un tampon
  		k = size[j];
  		tmp = new char[k];
  		for(i=0;i<k;++i)
  		{
  			tmp[i] = data[seq[j]+i];
  		}
  		/// On �crit la s�quence compl�mentaire invers�e dans le tableau data
  		l = 0;
  		for(i=k-1;i>=0;--i)
  		{
  			data[seq[j]+l] = complNT(tmp[i]);
  			++l;
  		}
  		delete [] tmp;
  	}
}
