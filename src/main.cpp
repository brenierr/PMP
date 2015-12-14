#include <iostream>
#include "BlackScholesModel.hpp"
#include "MonteCarlo.hpp"
#include <boost/program_options.hpp>
#include <mpi.h>
#include "parser.hpp"
#include "../src/Basket.hpp"
#include "../src/Asian.hpp"
#include "../src/Performance.hpp"
#include <string>
#include <omp.h>
#include <cmath>
#include <time.h>

namespace po = boost::program_options;
using namespace std;

// création d'un générateur aléatoire en parallèle en fonction du rang du processus courant 
PnlRng* create_pnl_rng_generator(int rank) 
{
    PnlRng* rng;
    rng = pnl_rng_dcmt_create_id(rank, PNL_RNG_MERSENNE);
    pnl_rng_sseed(rng, time(NULL));
    return rng;
}

/**
 * Le root lit le fichier infile et renvoie un pack dans buf (de taille bufsize) 
 */
int ReadAndPackInfileData(char* infile, char** buf, int *bufsize, BlackScholesModel** bs, MonteCarlo** mc, Option** option, int rank) {
    // Parsing du fichier infile
    Param *P = new Parser(infile);
    PnlVect *spot = pnl_vect_new(), *sigma = pnl_vect_new();
    double r;
    int size, nbIterations;
    double correlation;

    char* option_type;
    P->extract("option type", option_type);
    int num_option_type;

    if (std::strcmp (option_type, "basket") == 0)
    {
        num_option_type = 1;
        *option = new Basket(infile);
    }
    else if (std::strcmp (option_type, "asian") == 0)
    {
        num_option_type = 2;
        *option = new Asian(infile);
    }
    else if (std::strcmp (option_type, "performance") == 0)
    {
        num_option_type = 3;
        *option = new Performance(infile);
    }
    else
    {
        std::cout << "Type de l'option ("<<option_type<<") dans le fichier source invalide" << std::endl;
        std::cout << "Entrer un fichier correct " << std::endl;
        return 1;
    }
    size = (*option)->size_;
    P->extract("interest rate", r);
    P->extract("correlation", correlation);
    P->extract("sample number", nbIterations);
    P->extract("spot", spot, size);
    P->extract("volatility", sigma, size);

    // Pack de toutes les données à envoyer aux autres processus 
    *bufsize = 0;
    int info;
    int pos=0, count;
    // je les pack et je les envoie à mes fils
    // pack pour BlackScolesModel :  
    info = MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count); // taile de spot
    if (info) return info;
    *bufsize += count;
    info = MPI_Pack_size(spot->size, MPI_DOUBLE, MPI_COMM_WORLD, &count); // spot
    if (info) return info;
    *bufsize += count;
    info = MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count); // taille de sigma 
    if (info) return info;
    *bufsize += count;
    info = MPI_Pack_size(sigma->size, MPI_DOUBLE, MPI_COMM_WORLD, &count); // sigma
    if (info) return info;
    *bufsize += count;
    info = MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &count); // correlation
    if (info) return info;
    *bufsize += count;
    info = MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &count); // r 
    if (info) return info;
    *bufsize += count;
    info = MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count); // size
    if (info) return info;
    *bufsize += count;
    // pack pour Option : 
    info = MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count); // type de l'option
    if (info) return info;
    *bufsize += count;
    info = MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &count); // T
    if (info) return info;
    *bufsize += count;    
    info = MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count); // nbTimeSteps
    if (info) return info;
    *bufsize += count;
    info = MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count); // size
    if (info) return info;
    *bufsize += count;
    info = MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD, &count); // strike (on ne s'en sert pas dans le cas du perf)
    if (info) return info;
    *bufsize += count;
    info = MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count); // taile de weights
    if (info) return info;
    *bufsize += count;
    info = MPI_Pack_size((*option)->weights_->size, MPI_DOUBLE, MPI_COMM_WORLD, &count); // weights
    if (info) return info;
    *bufsize += count;
    // pack pour MonteCarlo : 
    info = MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &count); // nbIterations
    if (info) return info;
    *bufsize += count;

    // on alloue la place qu'il faut pour envoyer le pack
    *buf = (char*)malloc(*bufsize);

    // on met les informations dans le pack :
    // Pour BlackScholesModel : 
    info = MPI_Pack(&(spot->size), 1, MPI_INT, *buf, *bufsize, &pos, MPI_COMM_WORLD); // spot size
    if (info) return info;
    info = MPI_Pack((spot->array), spot->size, MPI_DOUBLE, *buf, *bufsize, &pos, MPI_COMM_WORLD); // spot array 
    if (info) return info;
    info = MPI_Pack(&(sigma->size), 1, MPI_INT, *buf, *bufsize, &pos, MPI_COMM_WORLD); // sigma size 
    if (info) return info;
    info = MPI_Pack((sigma->array), sigma->size, MPI_DOUBLE, *buf, *bufsize, &pos, MPI_COMM_WORLD); // sigma array 
    if (info) return info;
    info = MPI_Pack(&correlation, 1, MPI_DOUBLE, *buf, *bufsize, &pos, MPI_COMM_WORLD); // correlation
    if (info) return info;
    info = MPI_Pack(&r, 1, MPI_DOUBLE, *buf, *bufsize, &pos, MPI_COMM_WORLD); // r 
    if (info) return info;
    info = MPI_Pack(&size, 1, MPI_INT, *buf, *bufsize, &pos, MPI_COMM_WORLD); // size
    if (info) return info;
    // Pour Option :
    info = MPI_Pack(&(num_option_type), 1, MPI_INT, *buf, *bufsize, &pos, MPI_COMM_WORLD); // type de l'option
    if (info) return info;
    info = MPI_Pack(&((*option)->T_), 1, MPI_DOUBLE, *buf, *bufsize, &pos, MPI_COMM_WORLD); // T
    if (info) return info;
    info = MPI_Pack(&((*option)->nbTimeSteps_), 1, MPI_INT, *buf, *bufsize, &pos, MPI_COMM_WORLD); // nbTimeSteps
    if (info) return info;
    info = MPI_Pack(&((*option)->size_), 1, MPI_INT, *buf, *bufsize, &pos, MPI_COMM_WORLD); // size
    if (info) return info;
    if (std::strcmp(option_type, "performance") == 0) {
        double strike = 0.0;
        info = MPI_Pack(&strike, 1, MPI_DOUBLE, *buf, *bufsize, &pos, MPI_COMM_WORLD); // pour perf on met 0 
    } else if (std::strcmp(option_type, "basket") == 0)   {
        info = MPI_Pack(&(((Basket*)*option)->strike_), 1, MPI_DOUBLE, *buf, *bufsize, &pos, MPI_COMM_WORLD); // strike
    } else {
        info = MPI_Pack(&(((Asian*)*option)->strike_), 1, MPI_DOUBLE, *buf, *bufsize, &pos, MPI_COMM_WORLD); // strike
    }
    if (info) return info;
    info = MPI_Pack(&((*option)->weights_->size), 1, MPI_INT, *buf, *bufsize, &pos, MPI_COMM_WORLD); // taille de weights 
    if (info) return info;
    info = MPI_Pack((*option)->weights_->array, ((*option)->weights_->size), MPI_DOUBLE, *buf, *bufsize, &pos, MPI_COMM_WORLD); // weights 
    if (info) return info;
    // Pour MonteCarlo : 
    info = MPI_Pack(&nbIterations, 1, MPI_INT, *buf, *bufsize, &pos, MPI_COMM_WORLD); // nbIterations 
    if (info) return info;
   
    // création des objets pour le processeur de rank==0 
    *bs = new BlackScholesModel(size, r, correlation, sigma, spot);
    *mc = new MonteCarlo(*bs, *option, create_pnl_rng_generator(rank), 0.1, nbIterations);
    delete P;
}

/**
 * Les processus esclaves lisent le pack buf et crée les instances mc, bs et option
 */
int UnpackData(const char* buf,const int bufsize, BlackScholesModel** bs, MonteCarlo** mc, Option** option, int rank) {
    PnlVect *spot = pnl_vect_new(), *sigma = pnl_vect_new(), *weights = pnl_vect_new();
    int size, nbIterations, nbTimeSteps_option, size_option, num_type_option;
    double correlation, T_option, strike_option, r;
    // lire buf pour reconstruire BlackScholesModel :
    int pos = 0, info=0, n;
    info = MPI_Unpack(buf, bufsize, &pos, &n, 1, MPI_INT, MPI_COMM_WORLD); // taille de spot 
    if (info) return info;
    pnl_vect_resize(spot, n);
    info = MPI_Unpack(buf, bufsize, &pos, spot->array, n, MPI_DOUBLE, MPI_COMM_WORLD); // array de spot
    if (info) return info;
    info = MPI_Unpack(buf, bufsize, &pos, &n, 1, MPI_INT, MPI_COMM_WORLD); // taille de sigma 
    if (info) return info;
    pnl_vect_resize(sigma, n);
    info = MPI_Unpack(buf, bufsize, &pos, sigma->array, n, MPI_DOUBLE, MPI_COMM_WORLD); // array of sigma 
    if (info) return info;
    info = MPI_Unpack(buf, bufsize, &pos, &correlation, 1, MPI_DOUBLE, MPI_COMM_WORLD); // correlation 
    if (info) return info;
    info = MPI_Unpack(buf, bufsize, &pos, &r, 1, MPI_DOUBLE, MPI_COMM_WORLD); // r 
    if (info) return info;
    info = MPI_Unpack(buf, bufsize, &pos, &size, 1, MPI_INT, MPI_COMM_WORLD); // size 
    if (info) return info;
    // lire buf pour reconstruire Option :  
    info = MPI_Unpack(buf, bufsize, &pos, &num_type_option, 1, MPI_INT, MPI_COMM_WORLD); // type de l'option 
    if (info) return info;
    info = MPI_Unpack(buf, bufsize, &pos, &T_option, 1, MPI_DOUBLE, MPI_COMM_WORLD); // T
    if (info) return info;    
    info = MPI_Unpack(buf, bufsize, &pos, &nbTimeSteps_option, 1, MPI_INT, MPI_COMM_WORLD); // nbTimeSteps
    if (info) return info;
    info = MPI_Unpack(buf, bufsize, &pos, &size_option, 1, MPI_INT, MPI_COMM_WORLD); // size 
    if (info) return info;
    info = MPI_Unpack(buf, bufsize, &pos, &strike_option, 1, MPI_DOUBLE, MPI_COMM_WORLD); // strike 
    if (info) return info;
    info = MPI_Unpack(buf, bufsize, &pos, &n, 1, MPI_INT, MPI_COMM_WORLD); // taille de weights
    if (info) return info;
    pnl_vect_resize(weights, n);
    info = MPI_Unpack(buf, bufsize, &pos, weights->array, n, MPI_DOUBLE, MPI_COMM_WORLD); // array de weights 
    if (info) return info;
    // lire buf pour reconstruire MonteCarlo : 
    info = MPI_Unpack(buf, bufsize, &pos, &nbIterations, 1, MPI_INT, MPI_COMM_WORLD); // nbIterations 
    if (info) return info;
    // en fonction du type de l'option, on appelle le bon constructeur 
    if (num_type_option==1) {
        *option = new Basket(T_option, nbTimeSteps_option, size_option, strike_option, weights);
    } else if (num_type_option==2) {
        *option = new Asian(T_option, nbTimeSteps_option, size_option, strike_option, weights);
    } else {
        *option = new Performance(T_option, nbTimeSteps_option, size_option, weights);
    }

    // création des objets pour le processus courant 
    *bs = new BlackScholesModel(size, r, correlation, sigma, spot);
    *mc = new MonteCarlo(*bs, *option, create_pnl_rng_generator(rank), 0.1, nbIterations);
}

int main(int argc, char **argv)
{
    char* infile =  &(argv[1])[0u];
    Option *option;
    BlackScholesModel* bs;
    MonteCarlo *mc;

    int nb_threads, rank, bufsize;  

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nb_threads);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    char* buf;
    if (rank==0) { // je suis le maître, je lis les données du fichier et je les met dans buf 
        ReadAndPackInfileData(infile, &buf, &bufsize, &bs, &mc, &option, rank);
    }

    // Envoie du pack à tout le monde
    MPI_Bcast(&bufsize, 1, MPI_INT, 0, MPI_COMM_WORLD); // récupération de la taille du pack
    if (rank != 0) {
        buf = (char*)malloc(bufsize); // allocation de la taille qu'il faut pour récupérer le pack
    }
    MPI_Bcast(buf, bufsize, MPI_PACKED, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        UnpackData(buf, bufsize, &bs, &mc, &option, rank); // lit le pack buf et crée les objets bs et mc
    }
    // calcul du prix en 0
    double prix_MPI = 0.;
    double ic_MPI = 0.;   

    // les attributs liés à la précision 
    bool hasPrecision = (argc>2); // hasPrecision est faux si on garde le nombre de samples précisé dans le fichier
    double precision;
    if (hasPrecision) {
        precision = atof(argv[2]);
    }

    // rbuf = buffer de reception de tous les prix et var par le root
    double *rbuf;    
    if (rank==0) {
        rbuf = (double*)malloc(2*nb_threads*sizeof(double));
    }
    double sum=0., var=0., sum_total = 0., var_total = 0.;

    int minNbSamples = 10000*nb_threads; // nombre d'itérations minimales à séparer entre les threads
    int nbSamples = 10000*nb_threads; // pour commencer (premiere estimation)
    int nbSamplesTotal = nbSamples;

    double start=0., end=0.;
    
    start = MPI_Wtime();
    
    while(nbSamples > 0) { // tant qu'il reste du travail a faire pour atteindre la précision
        if(rank==0) {
            mc->PriceMaster(sum, var, nb_threads, rbuf, nbSamples); 
            sum_total += sum;
            var_total += var;
            prix_MPI = std::exp(-bs->r_*option->T_)*sum_total / nbSamplesTotal;
            ic_MPI = std::exp(-2*bs->r_*option->T_)*(var_total/nbSamplesTotal - std::pow(sum_total/nbSamplesTotal, 2));
            // estimation combien il en reste à faire pour atteindre la précision :
            nbSamples = std::max( (int) std::pow(2*1.96*sqrt(ic_MPI)/precision, 2), minNbSamples); // on ne fait pas moins de minSamples itérations (sinon le gather est trop long par rapport au travail fourni)
            ic_MPI = 2*1.96*std::sqrt(ic_MPI/nbSamplesTotal); 
            //std::cout << "nbSamples en plus : "<<nbSamples << std::endl;
            if (ic_MPI <= precision || !hasPrecision) { // on a fini si on est assez précis ou si on doit faire qu'un tour de boucle (nbSamples donné en paramètre) 
                nbSamples = 0;
            }
            nbSamplesTotal += nbSamples;
        } else {
            mc->PriceSlave(rank, nb_threads, rbuf, nbSamples);
        }
        // envoi du nombre de samples à réaliser pour le prochain tour
        MPI_Bcast(&nbSamples, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    end = MPI_Wtime();
    // on affiche le résultat pour un utilisateur : 
    /*
    if (rank==0) {
        cout << "*******************************" << endl;
        if (hasPrecision) {
            cout << "precision : " << precision << endl;
            cout << "nbSamplesTotal : " << nbSamplesTotal << endl;
        }
        else {
            cout << "nbSamples : " << nbSamplesTotal << endl;
        }
        cout << "Prix MPI : " << prix_MPI << endl;
        cout << "IC MPI : " <<  ic_MPI << endl;
        cout << "temps d'execution : " <<  end-start << endl;
    }*/
    // on affiche le résultat pour le script : 
    if (rank==0) {
        cout << nbSamplesTotal << " " << end-start << endl;    
    }

    if (rank==0) {
        free(rbuf); // seul le master l'a alloué 
    }

    pnl_rng_free(&(mc->rng_)); // (pas fait dans le destructeur de mc) 
    delete bs;
    delete mc;
    delete option;
    free(buf);

    MPI_Finalize();
    exit(0);
}


