#pragma once

#include "Option.hpp"
#include "BlackScholesModel.hpp"
#include "pnl/pnl_random.h"

/// \brief Modèle de Monte Carlo
class MonteCarlo
{
    public:
        /*! pointeur vers le modèle */
        BlackScholesModel *mod_;
        /*! pointeur sur l'option */
        Option *opt_;
        /*! pointeur sur le générateur */
        PnlRng *rng_; 
        /*! pas de différence finie */
        double fdStep_; 
        /*! nombre de tirages Monte Carlo */
        int nbSamples_; 

        /**
         * Constructeur de la classe MonteCarlo
         *
         * @param[in] mod le modèle de BlackSholees correspondant
         * @param[in] opt l'option à pricer
         * @param[in] rng le générateur de la libraire pnl
         * @param[in] fdStep le pas de différence finie
         * @param[in] nbSamples le nombres de tirage Monte Carlo
         */   
        MonteCarlo(BlackScholesModel *mod, Option *opt, PnlRng *rng, double fdStep, int nbSamples);

        /**
         * Destructeur de la classe MonteCarlo
         */
        ~MonteCarlo();

        /**
         * Calcule le prix de l'option à la date 0
         *
         * @param[out] prix valeur de l'estimateur Monte Carlo
         * @param[out] ic largeur de l'intervalle de confiance
         */
        void price(double &prix, double &ic);

        /**
         * Calcule le prix de l'option à la date t
         *
         * @param[in]  past contient la trajectoire du sous-jacent
         * jusqu'à l'instant t
         * @param[in] t date à laquelle le calcul est fait
         * @param[out] prix contient le prix
         * @param[out] ic contient la largeur de l'intervalle
         * de confiance sur le calcul du prix
         */
        void price(const PnlMat *past, double t, double &prix, double &ic);

        /**
         * Calcule le delta de l'option à la date t
         *
         * @param[in] past contient la trajectoire du sous-jacent
         * jusqu'à l'instant t
         * @param[in] t date à laquelle le calcul est fait
         * @param[out] delta contient le vecteur de delta
         * de confiance sur le calcul du delta
         */
        void delta(const PnlMat *past, double t, PnlVect *delta);

        /**
         * Calcule l'erreur du portefeuille de couverture
         * (Profit and Loss)
         * @param[in] H le nombre de dates de la simulation
         * \return l'erreur de couverture en double
         */
        double profitAndLoss(int H);


        void PriceSlave(int rank, int size, double* rbuf, int numberOfSamples); 
        bool PriceMaster(double&prix, double&ic, int size, double* rbuf, double precision, int numberOfSamplesFirstStep, int minSamples); 
        void PriceMaster(double&prix, double&var, int size, double* rbuf, int numberOfSamples);
};


