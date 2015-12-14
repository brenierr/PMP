#include "BlackScholesModel.hpp"
#include <cmath>
#include <iostream>

BlackScholesModel::BlackScholesModel(int size, double r, double rho, PnlVect *sigma, PnlVect *spot)
{
    size_ = size;
    r_ = r;
    rho_ = rho;
    sigma_ = sigma;
    spot_ = spot;
    chol_ = pnl_mat_create_from_scalar(size, size, rho_);
    chol(chol_);
    trend_ = pnl_vect_create_from_scalar(size_, r_);
    rVect_ = pnl_vect_create_from_scalar(size_, r_);
}


void BlackScholesModel::chol(PnlMat * M)
{
    for(int i = 0; i < size_; i++){
      pnl_mat_set(M, i, i, 1);
    }
    pnl_mat_chol(M);
}

void BlackScholesModel::asset(PnlMat* path, int ind_last, double T, int nbTimeSteps, PnlRng *rng, PnlVect *proba)
{
    PnlVect *G = pnl_vect_new();
    PnlVect Ld;
    double temps = T/nbTimeSteps;
    double taux;
    double brown;
    double sigmaD;
    double prodScal;

    for(int i = ind_last; i <= nbTimeSteps; i++){
        //Wt
        pnl_vect_rng_normal(G, size_, rng);
        for(int d = 0; d < size_; d++){
            Ld = pnl_vect_wrap_mat_row(chol_, d);
            prodScal = pnl_vect_scalar_prod(&Ld, G);
            sigmaD = pnl_vect_get(sigma_,d);
            taux = pnl_vect_get(proba,d) - sigmaD*sigmaD/2;
            brown = sigmaD*std::sqrt(temps)*prodScal;
            pnl_mat_set(path, i, d, pnl_mat_get(path, i-1, d)*std::exp(taux*temps+brown));
        }
   }
   pnl_vect_free(&G);
}

void BlackScholesModel::asset(PnlMat *path, double T, int nbTimeSteps, PnlRng *rng)
{
    for(int d = 0; d < size_; d++){
        pnl_mat_set(path, 0, d, pnl_vect_get(spot_, d));
    }
    asset(path, 1, T, nbTimeSteps, rng, rVect_);
}

void BlackScholesModel::asset(PnlMat *path, double t, int nbTimeSteps, double T,
           PnlRng *rng, const PnlMat *past)
{
    PnlVect currentRow;
    // copie des lignes de past dans path
    for(int i=0; i<past->m; i++)
    {
        currentRow = pnl_vect_wrap_mat_row(past, i);
        pnl_mat_set_row(path, &currentRow, i);
    }

    // calcul du prochain pas de discretisation à simuler
    double t_suiv = (past->m-1)*(T/nbTimeSteps);
    double delta_t = t_suiv - t;

    // calcul St, marche pour t un pas de discrétisation ou non.
    PnlVect *G = pnl_vect_new();
    PnlVect Ld;

    double taux;
    double brown;
    double sigmaD;
    double prodScal;
    pnl_vect_rng_normal(G, size_, rng);
    // calcul de la ligne à mettre dans path pour t = t_suiv (St_suiv)
    for(int d = 0; d < size_; d++){
        //Ld*Wt
        Ld = pnl_vect_wrap_mat_row(chol_, d);
        prodScal = pnl_vect_scalar_prod(&Ld, G);

        sigmaD = pnl_vect_get(sigma_,d);

        taux = r_ - sigmaD*sigmaD/2;

        brown = sigmaD*std::sqrt(delta_t)*prodScal;

        pnl_mat_set(path, past->m-1, d, pnl_mat_get(path, past->m-1, d)*std::exp(taux*delta_t+brown));
    }
    // remplie les lignes de path de t_suiv à T
    asset(path, past->m, T, nbTimeSteps, rng, rVect_);
    pnl_vect_free(&G);
}


void BlackScholesModel::shiftAsset(PnlMat *shift_path, const PnlMat *path,
                int d, double h, double t, double timestep, PnlVect *V)
{
    // copie de path dans shift_path
  //PnlVect *V = pnl_vect_new();//coute cher
  pnl_mat_clone(shift_path,path);
  //pnl_mat_print(shift_path);
  pnl_mat_get_col(V,shift_path,d);
  int start = ceil(t/timestep); 
  for( int i = start; i< path->m; i++ ) {
    pnl_vect_set(V,i,GET(V,i)*(1+h));
  } 
  pnl_mat_set_col(shift_path,V,d); 
  

}

void BlackScholesModel::simul_market(PnlMat *path, double T, int H, PnlRng *rng)
{
    for(int d = 0; d < size_; d++)
    {
        pnl_mat_set(path, 0, d, pnl_vect_get(spot_, d));
    }
    //proba historique est mise dans trend avec un set
    asset(path, 1, T, H, rng, trend_);
}

BlackScholesModel::~BlackScholesModel()
{
    pnl_mat_free(&chol_);
    pnl_vect_free(&trend_);
    pnl_vect_free(&rVect_);
}
