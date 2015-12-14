#include "Performance.hpp"

//#include <math.h>

/*!
 * \file Performance.cpp
 */

Performance::Performance(double T, int nbTimeSteps, int size, PnlVect* weights):Option(T, nbTimeSteps, size, weights){}

Performance::Performance(char* infile):Option(infile){}

double Performance::payoff(const PnlMat *path)
{
    PnlVect* Sti_d = pnl_vect_new();
    PnlVect* Sti_moins_1_d = pnl_vect_new();
    double numerateur, denominateur;
    double sum = 0;
    for(int i=1; i<nbTimeSteps_+1; i++)
    {
        pnl_mat_get_row(Sti_d, path, i); // On récupère la dernière ligne = Sti,d
        pnl_mat_get_row(Sti_moins_1_d, path, (i-1)); // On récupère l'avant dernière ligne = Sti-1,d

        numerateur = pnl_vect_scalar_prod(Sti_d, weights_);
        denominateur = pnl_vect_scalar_prod(Sti_moins_1_d, weights_);
        sum += std::max(0.0, numerateur/denominateur - 1);
    }

    // free memory
    pnl_vect_free(&Sti_d);
    pnl_vect_free(&Sti_moins_1_d);

    return (1 + sum);
}  

