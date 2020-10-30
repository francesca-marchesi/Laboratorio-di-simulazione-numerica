//Random numbers
#include "random.h"
#include <vector>
#include <fstream>

int seed[4];
Random rnd;

int N=30;
const int M=30;
double x[M], y[M];

double T_init;
double delta_T;
double step;
double beta;
double T;
int iterations;
double probability;
std::vector <int> old_vector;
//MPI constant
int my_size, my_rank;


//functions
double Boltzmann(int, std::vector<int>);
void Check(int, std::vector<std::vector<int>>);
void Circle(int, double);
void Crossover(double, std::vector<int>, std::vector<int>);
double Fitness(int, std::vector<int>);
void FinalConf(std::vector<int>);
void Input(void);
void Move(std::vector<int>, int, double);
void Square(int, double);
void WriteLenght(int, std::vector<int>);

void FinalConfXY(int);







std::vector<int> Mutazione1(double, std::vector<int>);
std::vector<int> Mutazione2(double, std::vector<int>);
//std::vector<int> Mutazione3(double, std::vector<int>);
std::vector<int> Mutazione4(double, std::vector<int>);
std::vector<int> Mutazione5(double, std::vector<int>);


























