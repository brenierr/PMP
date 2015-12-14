#include "Asian.hpp"

//#include <math.h>

/*!
 * \file Asian.cpp
 */

Asian::Asian(double T, int nbTimeSteps, int size, double strike,  PnlVect* weights):Option(T, nbTimeSteps, size, weights)
{
    strike_ = strike;
}

Asian::Asian(char* infile):Option(infile)
{
    Param *P = new Parser(infile);
    P->extract("strike", strike_);
    delete P;
}

double Asian::payoff(const PnlMat *path){
    double payoff = 0.0;
    PnlVect *prices = pnl_vect_create(path->m);
    for(int d=0; d<size_; d++) {
        pnl_mat_get_col(prices, path, d); // On récupère la dernière col = Sti,d
        payoff += pnl_vect_get(weights_, d) * (pnl_vect_sum(prices)) / (nbTimeSteps_ + 1); // lambda_d * sum(Sti,d)
    }
    payoff = payoff - strike_;

    // free memory
    pnl_vect_free(&prices);

    return (payoff>=0)?payoff:0;
}  

