#include "MonteCarlo.hpp"
#include <cmath>
#include <ctime>
#include "pnl/pnl_finance.h"
#include <mpi.h>

MonteCarlo::MonteCarlo(BlackScholesModel *mod, Option *opt, PnlRng *rng, double fdStep, int nbSamples)
{
    mod_ = mod;
    opt_ = opt;
    rng_ = rng;
    fdStep_ = fdStep;
    nbSamples_ = nbSamples;
}

void MonteCarlo::price(double &prix, double &ic)
{
    double sum = 0;
    double var = 0;
    double phi;
    PnlMat *path = pnl_mat_create(opt_->nbTimeSteps_+1, opt_->size_);
    for(int M = 0; M < nbSamples_; M++){
        mod_->asset(path, opt_->T_, opt_->nbTimeSteps_, rng_);
        phi = opt_->payoff(path);
        sum += phi;
        var += phi*phi;
    }
    prix = std::exp(-mod_->r_*opt_->T_) * sum / nbSamples_;
    ic = std::exp(-2*mod_->r_*opt_->T_)*(var/nbSamples_ - std::pow(sum/nbSamples_, 2));
    ic = 2.0*1.96*std::sqrt(ic/nbSamples_);

    pnl_mat_free(&path);
}

void MonteCarlo::price(const PnlMat *past, double t, double &prix, double &ic)
{
    double sum = 0;
    double var = 0;
    double phi;
    PnlMat *path = pnl_mat_create_from_scalar(opt_->nbTimeSteps_+1, opt_->size_,0);
    for(int M = 0; M < nbSamples_; M++){
        mod_->asset(path, t, opt_->nbTimeSteps_, opt_->T_, rng_, past);
        phi = opt_->payoff(path);
        sum += phi;
        var += phi*phi;
    }
    prix = std::exp(-mod_->r_*opt_->T_) * sum / nbSamples_;
    ic = std::exp(-2*mod_->r_*opt_->T_)*(var/nbSamples_ - std::pow(sum/nbSamples_,2));
    ic = 2.0*1.96*std::sqrt(ic/nbSamples_);
    pnl_mat_free(&path);
}

void MonteCarlo::delta(const PnlMat *past, double t, PnlVect *delta)
{

    double timestep = (opt_->T_/opt_->nbTimeSteps_);

    PnlMat *shift_path_1 =  pnl_mat_create(opt_->nbTimeSteps_ +1, opt_->size_);// la matrice déjà transposée
    PnlMat *shift_path_2 =  pnl_mat_create(opt_->nbTimeSteps_ +1, opt_->size_);
    PnlMat *path2 = pnl_mat_create(opt_->nbTimeSteps_ +1, opt_->size_);

    double phi_1;
    double phi_2;

    pnl_vect_resize(delta, opt_->size_);
    pnl_vect_set_all(delta, 0.0);
    PnlVect *V = pnl_vect_new();

    for (int M = 0; M < nbSamples_; ++M){
        mod_->asset(path2, t, opt_->nbTimeSteps_, opt_->T_, rng_, past);
        for (int d = 0; d < opt_->size_; ++d){

            mod_->shiftAsset(shift_path_1, path2, d, fdStep_, t, timestep, V); // h = fdStep
            mod_->shiftAsset(shift_path_2, path2, d, -fdStep_, t, timestep, V);

            phi_1 = opt_->payoff(shift_path_1);
            phi_2 = opt_->payoff(shift_path_2);

            pnl_vect_set(delta, d, pnl_vect_get(delta,d) + phi_1 - phi_2);

        }
    }
    PnlVect prixt = pnl_vect_wrap_mat_row(past, past->m-1);
    pnl_vect_mult_scalar(delta,  std::exp(-mod_->r_*(opt_->T_-t))/(nbSamples_*2*fdStep_));
    pnl_vect_div_vect_term(delta, &prixt);

    pnl_vect_free(&V);
    pnl_mat_free(&path2);
    pnl_mat_free(&shift_path_1);
    pnl_mat_free(&shift_path_2);
}


double MonteCarlo::profitAndLoss(int H)
{
    //condition inclusion sur la discretisation
    if(std::floor(H/opt_->nbTimeSteps_) != ((double)H/opt_->nbTimeSteps_))
        throw "condition d'inclusion des temps de discretisation non respectée";

    PnlMat *path = pnl_mat_create(H+1, opt_->size_);
    mod_->simul_market(path, opt_->T_, H, rng_);

    PnlMat *deltas = pnl_mat_create(H+1, opt_->size_);
    PnlVect *investTauxSansRisque = pnl_vect_create(H+1);

    PnlMat *past = pnl_mat_new();
    PnlMat *pastTot = pnl_mat_create(opt_->nbTimeSteps_ + 1, opt_->size_);

    //on recupere tous les 'N'timesteps  N = H/nbTimeSteps
    PnlVect *current = pnl_vect_new();
    for(int i = 0; i<= opt_->nbTimeSteps_; i++){
        pnl_mat_get_row(current, path, i*(H/opt_->nbTimeSteps_));
        pnl_mat_set_row(pastTot, current, i);
    }

    //le delta en 0 n'est pas le même on le traite separement
    pnl_mat_extract_subblock(past, path, 0, 1, 0, opt_->size_);

    PnlVect *delt = pnl_vect_create(opt_->size_);

    delta(past, 0, delt);

    pnl_mat_set_row(deltas, delt, 0);

    double prixOpt, ic;
    price(prixOpt, ic);

    PnlVect *prixSimu = pnl_vect_new();
    pnl_mat_get_row(prixSimu, path,0);
    pnl_vect_set(investTauxSansRisque,0, prixOpt-pnl_vect_scalar_prod(delt, prixSimu));// calcul de V0

    double coeff;
    double diff;
    PnlVect *deltAvant = pnl_vect_create(opt_->size_);
    int mult = H/opt_->nbTimeSteps_;
    int index;
    for(int i = 0; i< opt_->nbTimeSteps_; i++){
        //on extrait jusqu'a i et on parcourt les differences en h
        pnl_mat_extract_subblock(past, pastTot, 0, i+2, 0, opt_->size_);

        for(int j=1; j<= mult; j++){
            double t = i*opt_->T_/opt_->nbTimeSteps_ + j*opt_->T_/H;
            index = mult*i+j;
            pnl_mat_get_row(current, path, index);
            pnl_mat_set_row(past, current, past->m-1);

            delta(past, t,delt);

            pnl_mat_set_row(deltas, delt, index);

            pnl_mat_get_row(prixSimu, path,index);
            coeff = std::exp(mod_->r_*opt_->T_/H);
            pnl_mat_get_row(deltAvant, deltas, index-1);
            pnl_vect_minus_vect(delt, deltAvant);
            diff = pnl_vect_scalar_prod(delt ,prixSimu); // diff = (delta(i) - delta(i-1))*Sti
            pnl_vect_set(investTauxSansRisque, index, pnl_vect_get(investTauxSansRisque,index-1)*coeff-diff);
        }
    }

    double VH = pnl_vect_get(investTauxSansRisque, H);
    pnl_mat_get_row(delt, deltas, H);
    pnl_mat_get_row(prixSimu, path, H);
    double actifRisque = pnl_vect_scalar_prod(delt ,prixSimu);

    double pl = VH + actifRisque - opt_->payoff(pastTot);

    pnl_mat_free(&path);
    pnl_mat_free(&deltas);
    pnl_vect_free(&investTauxSansRisque);
    pnl_mat_free(&past);
    pnl_mat_free(&pastTot);
    pnl_vect_free(&current);
    pnl_vect_free(&delt);
    pnl_vect_free(&prixSimu);
    pnl_vect_free(&deltAvant);

    return pl;
}

void MonteCarlo::PriceMaster(double&sum, double&var, int size, double* rbuf, int numberOfSamples)
{
    PriceSlave(0, size, rbuf, numberOfSamples);
    sum=0.;
    var=0.;
    // j=0 et j=1 sont pour le thread root (qui n'a pas travaillé)
    for(int j=0; j<2*size; j=j+2) {
        sum += rbuf[j];
        var += rbuf[j+1];
    }
}

void MonteCarlo::PriceSlave(int rank, int nb_threads, double* rbuf, int numberOfSamples) 
{
    double prix_t=0., var_t=0.;

    double payoff;
    PnlMat *pathMatrix = pnl_mat_create(opt_->nbTimeSteps_+1, opt_->size_);
    for (size_t i = 0 ; i < numberOfSamples/(nb_threads); i++)
    {
        mod_->asset(pathMatrix, opt_->T_, opt_->nbTimeSteps_, rng_);
        payoff = opt_->payoff(pathMatrix);
        prix_t += payoff;
        var_t += payoff * payoff;
    }
    double buffer[2];
    buffer[0] = prix_t;
    buffer[1] = var_t;

    MPI_Gather(buffer, 2, MPI_DOUBLE, rbuf, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    pnl_mat_free(&pathMatrix);
}


MonteCarlo::~MonteCarlo()
{
}
