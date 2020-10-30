/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
double walker;

// averages
double blk_av,blk_norm,accepted,attempted, contatore, contatore2;
double glob_av,glob_av2;
double stima_ham,err_ham;

// simulation
int nstep, nblk, step_T, step_opt;
double delta, mu_opt, sigma_opt, ene_opt, x0, x1, init_mu, init_sigma, delta_opt, delta_T, T;

//functions
void Input(void);
void InputOptimization(void);
void Optimization(void);
void Reset(int);
void Accumulate(void);
void OptimalParameter(int);
void Averages(int);
void Move(double, double);
void Measure(double, double);
double H_mean(double, double, double);
double phi(double, double, double);
double Energy(double, double, double);
double Error(double,double,int);

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
