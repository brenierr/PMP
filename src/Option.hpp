#pragma once
#include <iostream>
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"
#include "parser.hpp"


/*! \class Option
 * \brief classe representant une option financière
 *
 *  La classe gere la lecture des données à partir
 * d'un fichier et le calcul du payoff de l'option
 */
class Option
{

public:
 
  double T_; /*! maturité */
  int nbTimeSteps_; /*! nombre de pas de temps de discrétisation */
  int size_; /*!  dimension du modèle, redondant avec BlackScholesModel::size_ */
  PnlVect *weights_; /*! les poids des sous-jacents */

public:
 
  /*!
   *  \brief Constructeur
   *
   *  Constructeur de la classe CPlayer
   *
   *  \param T : maturité 
   *  \param size : nombre d'actifs sous-jacents
   *  \param weights : vecteurs de poids des sous-jacents
   */
  Option(double T, int nbTimeSteps, int size, PnlVect* weights);
 
  /*!
   *  \brief Destructeur
   *
   *  Destructeur de la classe Option
   */
  ~Option();

  /*!
   *  \brief Constructeur
   *
   *  Constructeur de la classe CPlayer
   *
   *  \param infile : le fichier à parser
   */ 
  Option(char* infile);

  /*!
   * \brief Calcul la valeur du payoff sur la trajectoire
   *
   * @param[in] path est une matrice de taille N+1 x d
   * contenant une trajectoire du modèle telle que créée
   * par la fonction asset.
   * \return phi(trajectoire)
   */
  virtual double payoff(const PnlMat *path) = 0;

};


