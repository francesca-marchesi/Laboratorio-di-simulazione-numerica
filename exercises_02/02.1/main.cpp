#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <math.h>

using namespace std;

int main (int argc, char *argv[]){
    
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
    ofstream ES1;
    ofstream ES2;
    
    
    int M=10000000;
    int N=100;
    int N_block=M/N;
    double a=0.;
    double b=1.;
    double x;
    double y;
    double f=0.;
    double g=0.;
    double I_blocco= 0.;
    double stima_I =0.;
    double stima_I2 =0.;
    double error_I = 0.;
    
    //distribuzione uniforme
    ES1.open("1.1");
    for(int j=0; j<N; j++){
        I_blocco=0.;
        f=0.;
    //stima di I nel blocco
        for(int i=0; i<N_block; i++){
            x=rnd.Rannyu(a, b);
            f=f+M_PI/2*cos(M_PI*x/2);
        }
        I_blocco=f*(b-a)/N_block;
        cout<<"l'integrale di blocco vale "<<I_blocco<<endl;
   
        stima_I += I_blocco;
        stima_I2 += I_blocco*I_blocco;
        error_I = sqrt((stima_I2/(j+1) - pow(stima_I/(j+1),2))/(j+1));
        ES1<<stima_I/(j+1)<<" "<<error_I<<endl;
        
    }
    ES1.close();
    cout<< "l'integrale vale "<<stima_I/(N)<<endl;
    
    //importance sampling
    I_blocco=0.;
    stima_I = 0.;
    stima_I2 = 0.;
    error_I = 0.;
    ES2.open("1.2");
    for(int j=0; j<N; j++){
        I_blocco=0.;
        g=0.;
        
        for(int i=0; i<N_block; i++){
            y=rnd.Rannyu();
            x=1-sqrt(1-y);
            g=g+(M_PI/2*cos(M_PI*x/2)/(2-2*x));
        }
        I_blocco=g/(double)N_block;
        cout<<"l'integrale di blocco vale "<<I_blocco<<endl;
        
        stima_I += I_blocco;
        stima_I2 += I_blocco*I_blocco;
        error_I = sqrt((stima_I2/(j+1) - pow(stima_I/(j+1),2))/(j+1));
        ES2<<stima_I/(j+1)<<" "<<error_I<<endl;
        

    }
    
    cout<< "l'integrale vale "<<stima_I/(N)<<endl;
    return 0;
}

























