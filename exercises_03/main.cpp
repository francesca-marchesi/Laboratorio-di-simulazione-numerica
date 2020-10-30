#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <algorithm>
#include <cmath>

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
    
    
    int N=100000;
    int N_blocchi=100;
    int M=N/N_blocchi;
    
    double S_0=100.;
    double T=1.;
    double K=100.;
    double r=0.1;
    double sigma=0.25;
    
    
//_______________Black-Scholes__________________
    
   /* double PUT=0.;
    double CALL=0.;
    double d1;
    double d2;
    double N1;
    double N2;
    
    d1= (log(S_0/K)+(r+sigma*sigma*T/2))/(sigma*sqrt(T));
    d2= d1-sigma*sqrt(T);
    
    N1=0.5*(1+erf(d1/sqrt(2)));
    N2=0.5*(1+erf(d2/sqrt(2)));
    
    CALL=S_0*N1-K*exp(-r*T)*N2;
    PUT=S_0*(N1-1)-K*exp(-r*T)*(N2-1);*/
    
    
    
//_______________parte uno: simulo solo S(T)_________________________
    /*double w=0.; //variabile distribuita come N(0,t)
    double S_T=0.;
    double Pi_put=0.;
    double Pi_call=0.;
    double call=0.;
    double put=0.;
    double Prezzo_put[N_blocchi];
    double Prezzo_call[N_blocchi];
    double Incertezza_put[N_blocchi];
    double Incertezza_call[N_blocchi];
    double prezzo_put=0.;
    double prezzo_call=0.;
    double prezzo_put2=0.;
    double prezzo_call2=0.;
    
    for(int s=0; s<N_blocchi;s++){
        Incertezza_put[s]=0.;
        Incertezza_call[s]=0.;
    }
    
    ofstream Put1, Call1;
    Put1.open("Put1");
    Call1.open("Call1");
    
    for (int j=0; j<N_blocchi; j++){
        Pi_put=0.;
        Pi_call=0.;

        for(int i=0; i<M; i++){
            w=rnd.Gauss(0, T);
            S_T=S_0*exp((r-pow(sigma,2)/2)*T+sigma*w);
        
            call=exp(-r*T)*max(S_T-K,0.);
            put=exp(-r*T)*max(K-S_T,0.);

            Pi_call=Pi_call+call; //dopo aver ciclato sul singolo blocco questa variabile contiene la somma di tutti i prezzi delle call nel blocco
            Pi_put=Pi_put+put;
            //cout<<Pi_call<<endl;
        }
    
    
    
        prezzo_put=prezzo_put+Pi_put/M;
        Prezzo_put[j]=prezzo_put/(j+1);
        prezzo_put2=prezzo_put2+pow(Pi_put/M,2);
        //cout<<prezzo_put<<endl;
        cout<<Prezzo_put[j]<<endl;
        
        prezzo_call=prezzo_call+Pi_call/M;
        Prezzo_call[j]=prezzo_call/(j+1);
        prezzo_call2=prezzo_call2+pow(Pi_call/M,2);
        
        if(j!=0){
            Incertezza_put[j]=sqrt((prezzo_put2/(j+1)-pow(prezzo_put/(j+1), 2))/j);
            Incertezza_call[j]=sqrt((prezzo_call2/(j+1)-pow(prezzo_call/(j+1), 2))/j);
        }
        cout<<Incertezza_put[j]<<endl;

        
    }
    
     for (int k=0; k<N_blocchi; k++){
         Put1<<Prezzo_put[k]<<" "<<Incertezza_put[k]<<endl;
         Call1<<Prezzo_call[k]<<" "<<Incertezza_call[k]<<endl;
     }
    Put1.close();
    Call1.close();

    */
    
//_______________parte due: simulo tempi discreti_______________________________
    int N_passi=100;
    double z=0.;
    double S_ti;
    double price_call=0.;
    double price_put=0.;
    double call_pr=0.;
    double put_pr=0.;
    double call_pr2=0.;
    double put_pr2=0.;
    double final_put[N_blocchi];
    double final_call[N_blocchi];
    double Inc_put[N_blocchi];
    double Inc_call[N_blocchi];
 
 ofstream Put2, Call2;
 Put2.open("Put2");
 Call2.open("Call2");

    for(int s=0; s<N_blocchi;s++){
        final_put[s]=0.;
        final_call[s]=0.;
        Inc_put[s]=0.;
        Inc_call[s]=0.;
        
    }
    
    for(int l=0; l<N_blocchi; l++){
    
        price_call=0.;
        price_put=0.;
    
        for(int s=0; s<M; s++){
            S_ti=S_0;
    
            for(int i=0; i<N_passi; i++){
                z=rnd.Gauss(0, 1);
                S_ti=S_ti*exp((r-pow(sigma,2)/2)/100+sigma*z*sqrt(1./100.));
                }
        
            price_call=price_call+exp(-r*T)*max(S_ti-K,0.);
            price_put= price_put+exp(-r*T)*max(K-S_ti,0.);
        }
        
        call_pr=call_pr+price_call/M;
        call_pr2=call_pr2+pow(price_call/M, 2);
        final_call[l]=call_pr/(l+1);
        
        put_pr=put_pr+price_put/M;
        put_pr2=put_pr2+pow(price_put/M, 2);
        final_put[l]=put_pr/(l+1);
        
        
        if(l!=0){
            Inc_put[l]=sqrt((put_pr2/(l+1)-pow(put_pr/(l+1), 2))/l);
            Inc_call[l]=sqrt((call_pr2/(l+1)-pow(call_pr/(l+1), 2))/l);
        }
    
    }
 for (int k=0; k<N_blocchi; k++){
     Put2<<final_put[k]<<" "<<Inc_put[k]<<endl;
     Call2<<final_call[k]<<" "<<Inc_call[k]<<endl;
 }
    Put2.close();
    Call2.close();

    return 0;
}
