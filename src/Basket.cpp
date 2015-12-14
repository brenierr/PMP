#include "Basket.hpp"

Basket::Basket(double T, int nbTimeSteps, int size, double strike,  PnlVect* weights):Option(T, nbTimeSteps, size, weights)
{
    strike_ = strike;
}

Basket::Basket(char* infile):Option(infile)
{
    Param *P = new Parser(infile);
    P->extract("strike", strike_);
    delete P;
}

double Basket::payoff(const PnlMat *path)
{
    double s=0;
    for ( int d=0; d<size_; d++)
        s += (pnl_mat_get(path,nbTimeSteps_,d)*pnl_vect_get(weights_,d));

    double payoff = s - strike_;

    return (payoff>=0)?payoff:0;
}
