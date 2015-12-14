#ifndef Performance_H
#define Performance_H


#include "Option.hpp"
#include "parser.hpp"

/*! \class Performance
 * \brief classe representant une option type Performance
 *
 *  sous-classe de la classe Option
 *  La classe gere la lecture des données à partir
 * d'un fichier et le calcul du payoff de l'option
 */
class Performance : public Option
{

public:
  /*!
   *  \brief Constructeur
   *
   *  Constructeur de la classe Performance
   *
   *  \param T : maturité 
   *  \param size : nombre d'actifs sous-jacents
   *  \param weights : vecteurs de poids des sous-jacents
   */
  Performance(double T, int nbTimeSteps, int size, PnlVect* weights);
 
  /*!
   *  \brief Constructeur
   *
   *  Constructeur de la classe Performance
   *
   *  \param infile : le fichier à parser
   */ 
  Performance(char* infile);

  /*!
   * \brief Calcul la valeur du payoff de Performance sur la trajectoire
   *
   * @param[in] path est une matrice de taille N+1 x d
   * contenant une trajectoire du modèle telle que créée
   * par la fonction asset.
   * \return phi(trajectoire)
   */ 
  double payoff(const PnlMat *path);
};


#endif
