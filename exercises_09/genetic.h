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
#include <vector>
#include <fstream>

int seed[4];
Random rnd;

std::vector<std::vector<int>> sample;
std::vector<std::vector<int>> sons;
std::vector<int> v1;
std::vector<int> v2;
std::vector<double> dist;
std::vector<double> prob;


int N=30;
int n=10;
int N_generations=2500;
const int M=30;
double x[M], y[M];
int c[2];
int check[2];

//functions

void Check(int, std::vector<std::vector<int>>);
void Circle(int, double);
void Crossover(double, std::vector<int>, std::vector<int>);
double Fitness(int, int);
void Generatrice(int, int);
void Input(void);
void Extract(int);
void Selection(int, int);
void Square(int, double);
void WriteLenght(int, std::vector<double>);

void FinalConfXY(int);



std::vector<int> Change(int, int, int, std::vector<int>);
std::vector<int> Mutazione1(double, std::vector<int>);
std::vector<int> Mutazione2(double, std::vector<int>);
//std::vector<int> Mutazione3(double, std::vector<int>);
std::vector<int> Mutazione4(double, std::vector<int>);
std::vector<int> Mutazione5(double, std::vector<int>);
/****************************************************************
 *****************************************************************
 _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
 _/_/  _/ _/       _/       Physics Department
 _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
 _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/























