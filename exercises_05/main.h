/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include "random.h"
#include <fstream>


int seed[4];
Random rnd;

//variable
double x_0;
double y_0;
double z_0;
double X_0;
double Y_0; //coordinate del punto di partenza
double Z_0;
double X;
double Y;
double Z;

double x;
double y;
double z;
double costGS;
double cost2p;

int Nsim;
int Nblocchi;






//function
void Input(void);
double MoveGS(double, double, double, int);
double MoveEX(double, double, double, int);

double GS(double, double);
double EXCITED(double, double, double, double);





















/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
