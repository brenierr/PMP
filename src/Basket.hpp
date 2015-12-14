#ifndef BasketH
#define BasketH

#include "Option.hpp"
#include "parser.hpp"

/*! \class Basket
 * \brief classe representant une option type Basket
 *
 *  sous-classe de la classe Option
 *  La classe gere la lecture des données à partir
 * d'un fichier et le calcul du payoff de l'option
 */
class Basket : public Option
{

public :
  double strike_; /*! prix d'exercice */

  /*!
   *  \brief Constructeur
   *
   *  Constructeur de la classe Basket
   *
   *  \param T : maturité 
   *  \param size : nombre d'actifs sous-jacents
   *  \param stirke : prix d'exercice
   *  \param weights : vecteurs de poids des sous-jacents
   */
  Basket(double T, int nbTimeSteps, int size, double strike,  PnlVect* weights);
    
  /*!
   *  \brief Constructeur
   *
   *  Constructeur de la classe Basket
   *
   *  \param infile : le fichier à parser
   */ 
  Basket(char* infile);
  
  /*!
   * \brief Calcul la valeur du payoff de Basket sur la trajectoire
   *
   * @param[in] path est une matrice de taille N+1 x d
   * contenant une trajectoire du modèle telle que créée
   * par la fonction asset.
   * \return phi(trajectoire)
   */ 
  double payoff(const PnlMat *path);

};

#endif
