#include "Option.hpp"


Option::Option( double T, int nbTimeSteps, int size, PnlVect* weights){

  T_ = T;
  nbTimeSteps_ = nbTimeSteps;
  size_ = size;
  weights_ = weights;
}



Option::Option(char* infile){

    Param *P = new Parser(infile);

    P->extract("maturity", T_);
    P->extract("option size", size_);
    P->extract("timestep number", nbTimeSteps_);
    P->extract("payoff coefficients", weights_, size_);  
    delete P;
}

Option::~Option()
{
    pnl_vect_free(&weights_);
}
