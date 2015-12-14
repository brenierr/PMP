#pragma once

#include "pnl/pnl_random.h"
#include "pnl/pnl_vector.h"
#include "pnl/pnl_matrix.h"

/// \brief Modèle de Black Scholes
class BlackScholesModel
{
public:
  /*! nombre d'actifs du modèle*/
  int size_; 
  /*! taux d'intérêt*/
  double r_; 
  PnlVect *rVect_;
  /*! paramètre de corrélation*/
  double rho_;
  /*! vecteur de volatilités*/
  PnlVect *sigma_; 
  /*! valeurs initiales du sous-jacent*/
  PnlVect *spot_;
  /*! proba historique*/
  PnlVect *trend_;
  /*! matrice de correlation */
  PnlMat *chol_;
 
  /**
   * Constructeur de la classe BlackSholesModel
   *
   * @param[in] size le nombre d'actifs sous-jacents
   * @param[in] r le taux d'intérêt
   * @param[in] rho le coefficient de corrélation
   * @param[in] sigma le vecteur des volatilités
   * @param[in] spot le vecteur des valeurs initiales du sous-jacent
   */   
  BlackScholesModel(int size, double r, double rho, PnlVect *sigma, PnlVect *spot);
    
  /**
   * Destructeur de la classe BlackSholes
   */
  ~BlackScholesModel();
  /**
   *Créé la matrice de cholesky a partir de rho_
   * @param[out] M la futur matrice de cholesky
   */
  void chol(PnlMat *M);

  /**
   * Génère une trajectoire du modèle et la stocke dans path
   *
   * @param[out] path contient une trajectoire du modèle.
   * C'est une matrice de taille (N+1) x d
   * @param[in] T  maturité
   * @param[in] nbTimeSteps nombre de dates de constatation
   */
  void asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng);

  /**
   * Calcule une trajectoire du sous-jacent connaissant le
   * passé jusqu' à la date t
   *
   * @param[out] path  contient une trajectoire du sous-jacent
   * donnée jusqu'à l'instant T par la matrice past
   * @param[in] t date jusqu'à laquelle on connait la trajectoire.
   * t n'est pas forcément une date de discrétisation
   * @param[in] nbTimeSteps nombre de pas de constatation
   * @param[in] T date jusqu'à laquelle on simule la trajectoire
   * @param[in] past trajectoire réalisée jusqu'a la date t
   */
  void asset(PnlMat *path, double t, int nbTimeSteps, double T,
	     PnlRng *rng, const PnlMat *past);

  /**
   * Shift d'une trajectoire du sous-jacent
   *
   * @param[in]  path contient en input la trajectoire
   * du sous-jacent
   * @param[out] shift_path contient la trajectoire path
   * dont la composante d a été shiftée par (1+h)
   * à partir de la date t.
   * @param[in] t date à partir de laquelle on shift
   * @param[in] h pas de différences finies
   * @param[in] d indice du sous-jacent à shifter
   * @param[in] timestep pas de constatation du sous-jacent
   */
  void shiftAsset(PnlMat *shift_path, const PnlMat *path,
		  int d, double h, double t, double timestep, PnlVect *V);

  /**
   * Génère une trajectoire du modèle avec la proba historique et la stocke dans path
   *
   * @param[out] path contient une trajectoire du modèle.
   * C'est une matrice de taille (N+1) x d
   * @param[in] T  maturité
   * @param[in] H nombre de dates de constatation
   */
  void simul_market(PnlMat *path, double T, int H, PnlRng *rng);


  /**
   * Met le vecteur des proba historiques à trend 
   */  
  void setTrend(PnlVect *trend){
    trend_ = trend;
  }

private:
    /**
     * Calcule une trajectoire du sous-jacent connaissant le
     * passé jusqu' à la date t
     *
     * @param[in/out] path  contient une trajectoire du sous-jacent
     * donnée jusqu'à l'instant T par la matrice past. En entrée, elle contient les trajectoires
     * du sous-jacent jusqu'à l'indice ind_last
     * @param[in] ind_last date à partir de laquelle il faut remplir la trajectoire
     * t n'est pas forcément une date de discrétisation
     * @param[in] nbTimeSteps nombre de pas de constatation
     * @param[in] T date jusqu'à laquelle on simule la trajectoire
     * @param[in] proba : probabilité à utiliser (historique ou r)
     */
    void asset(PnlMat* path, int ind_last, double T, int nbTimeSteps, PnlRng *rng, PnlVect *proba);


};


