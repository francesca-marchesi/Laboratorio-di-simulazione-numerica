//
//  main.cpp
//  01.3
//
//  Created by Francesca on 04/05/20.
//  Copyright © 2020 Francesca. All rights reserved.
//

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

    double L;
    double d;
    double pi=0.;
    double pi2;
    int N_tot=1000000;
    int N_block=100;
    int N_int=N_tot/N_block;
    double pi_blocco;
    double pi2_blocco;
    double sigma_pi=0.;
    double prova=0.;
    double prova_2=0;
    L=5.;
    d=15.;
    
    double r;
    double s;
    double x;
    double sum;
    int N_hit=0;
    
    ofstream Buffon;
    Buffon.open("Buffon2");
    
    for(int j=0; j<N_block; j++){
        pi_blocco=0.;
        pi2_blocco=0.;
        N_hit=0.;
        for(int i=0; i<N_int; i++){
            r=rnd.Rannyu(0, d);
            s=rnd.Rannyu(0, 360);
            x=L*cos(s);
            sum=r+x;
            if (sum<0. || sum>d){
                N_hit=N_hit+1;
            }
        }
        
        pi_blocco=(2*L*N_int)/(N_hit*d);
        
        pi2_blocco=pi_blocco*pi_blocco;
        
        prova= prova+ pi_blocco;
        prova_2=prova_2+pi2_blocco;
        pi=prova/(j+1);
        pi2=prova_2/(j+1);
        sigma_pi=sqrt((pi2-pi*pi)/(N_block));
        Buffon<<pi<<" "<<sigma_pi<<endl;
        
    }
    
    cout<<"pigreco "<<pi<<" sigma_pi "<<sigma_pi<<endl;
    Buffon.close();
    
    return 0;
}




































