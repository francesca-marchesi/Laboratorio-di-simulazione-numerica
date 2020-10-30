#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
using namespace std;

int main(int argc, const char * argv[]) {
    
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();
    
    ifstream input("seed.in");
    string property;
    if (input.is_open()){
        while ( !input.eof() ){
            input >> property;
            if( property == "RANDOMSEED" ){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;
    rnd.SaveSeed();

    double gamma=1.;
    double lambda=1.;
    double mu=0.;
    int n=10000;
    int N=1;
    
    double y=0.;
    double z=0.;
    
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
//distribzione di lorentz
    ofstream Lorentz;
    double S;
    double r;
    double min=0.;
    double max=0.;
    int num_bin=1000;
    double bin=0.;
    int prova=0;
    double lor[n];
    int isto_lor[num_bin];

    Lorentz.open("HistoL.100");
    for(int k=0; k<num_bin; k++){
        isto_lor[k]=0;
    }
    for (int i=0; i<n; i++) {
        S=0.;
        for (int j=0; j<N; j++) {
            r=rnd.Rannyu();
            y=gamma*tan(M_PI*(r-(1./2.)));
            S=S+y;
        }
        lor[i]=S/N;
        if(lor[i]>max){max=lor[i];}
        if(lor[i]<min){min=lor[i];}
    }

    bin=(max-min)/num_bin;
    
    for(int k=0; k<M; k++){
        prova=(lor[k]-min)/bin;
        isto_lor[prova]=isto_lor[prova]+1;
    }
    

    for(int k=0; k<num_bin; k++){
        Lorentz<<isto_lor[k]<<endl;
    }

    for(int k=0; k<n; k++){
        Lorentz<<lor[k]<<endl;
    }
    Lorentz.close();
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
// distribuzione esponenziale
    ofstream Exp;
    double isto_exp[n];
    double S2;
    double s;
    Exp.open("HistoE.100");
    
    for (int i=0; i<n; i++) {
        S2=0.;
        
        for (int j=0; j<N; j++) {
            s=rnd.Rannyu();
            z=-1/lambda*log(1-s);
            S2=S2+z;
        }
        isto_exp[i]=S2/N;
    }
    for(int k=0; k<n; k++){
        Exp<<isto_exp[k]<<endl;
    }

    Exp.close();
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
// distribuzione uniforme
    ofstream Uni;

    double isto_uni[n];
    double S3;
    double t;
    
    Uni.open("HistoU.1");

    
    for (int i=0; i<n; i++) {
        S3=0.;
        
        for (int j=0; j<N; j++) {
            t=rnd.Rannyu();
            S3=S3+t;
        }
        isto_uni[i]=S3/N;
    }
    for(int k=0; k<n; k++){
        Uni<<isto_uni[k]<<endl;
    }

    Uni.close();
   return 0;
}
