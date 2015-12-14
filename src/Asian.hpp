#ifndef Asian_H
#define Asian_H
 
#include "Option.hpp"
#include "parser.hpp"

/*! \class Asian
 * \brief classe representant une option type Asian
 *
 *  sous-classe de la classe Option
 *  La classe gere la lecture des données à partir
 * d'un fichier et le calcul du payoff de l'option
 */
class Asian : public Option
{

public:

  double strike_; /*! prix d'exercice */

  /*!
   *  \brief Constructeur
   *
   *  Constructeur de la classe Asian
   *
   *  \param T : maturité 
   *  \param size : nombre d'actifs sous-jacents
   *  \param stirke : prix d'exercice
   *  \param weights : vecteurs de poids des sous-jacents
   */
  Asian(double T, int nbTimeSteps, int size, double strike,  PnlVect* weights);

  /*!
   *  \brief Constructeur
   *
   *  Constructeur de la classe Asian
   *
   *  \param infile : le fichier à parser
   */ 
  Asian(char* infile);

  /*Destructeur*/
  // ~Asian();

  /*!
   * \brief Calcul la valeur du payoff de Asian sur la trajectoire
   *
   * @param[in] path est une matrice de taille N+1 x d
   * contenant une trajectoire du modèle telle que créée
   * par la fonction asset.
   * \return phi(trajectoire)
   */
  double payoff(const PnlMat *path);
 
};


#endif
