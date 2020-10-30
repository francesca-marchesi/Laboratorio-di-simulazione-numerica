#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
using namespace std;

int main(int argc, const char * argv[]) {
 
    int M=1000000000;
    int N=100;
    int lanci=M/N;
    
    
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
    
    
    ofstream output;
    ofstream output_prova;
    output.open("ese_1.dat");
    output_prova.open("prova");
    
    //esercizio 1.1.1
    double incertezza_stat[N];
    incertezza_stat[0]=0;
    double vettore[N];
    double r2_medio[N];
    double mu_r[N];
    double sum_mu=0.;
    double sum_mu2=0.;
    
    for (int j=0; j<N; j++){
        double stima_mu=0.;
        double mu_i=0.;

        
        for(int i=0; i<lanci; i++){
            stima_mu +=rnd.Rannyu();
        }
        
        mu_i=stima_mu/(double)lanci;
        sum_mu += mu_i ;
        sum_mu2 += mu_i*mu_i;
        
        mu_r[j]=sum_mu/(double)(j+1);
        r2_medio[j]=sum_mu2/(double)(j+1);
        vettore[j]=(sum_mu/(double)(j+1))-0.5;
        
        if(j!=0){
            incertezza_stat[j]=sqrt((r2_medio[j]-pow(mu_r[j],2))/(double)(j+1));
        }
    }

    for (int j=0; j<N; j++){
        output<<mu_r[j]<< "\t"<<incertezza_stat[j]<<endl;
        output_prova<<vettore[j]<<endl;
    }
    
    output.close();
    output_prova.close();
    
    
    
    
    //esercizio 1.1.2
    ofstream output2;
    output2.open("ese_2.dat");

   
    double sum_sigma=0.;
    double sum_sigma2=0.;
    double sigma_medio[N];
    double sigma2_medio[N];
    double vettore2[N];
    double incertezza_sigma[N];
    double r;
    incertezza_sigma[0]=0.;
    
    for (int j=0; j<N; j++){
        double stima_sigma=0.;
        double sigma_i=0.;
        
        
        for(int i=0; i<lanci; i++){
            r=rnd.Rannyu();
            stima_sigma+=pow(r-0.5,2);
        }
        
        sigma_i =stima_sigma/(double)lanci;
        sum_sigma += sigma_i;
        sum_sigma2+= sigma_i*sigma_i;
        
        sigma_medio[j]=sum_sigma/(j+1);
        sigma2_medio[j]=sum_sigma2/(j+1);
        vettore2[j]=(sum_sigma/(j+1))-(1./12.);
        
        if(j!=0){
            incertezza_sigma[j]=sqrt((sigma2_medio[j]-pow(sigma_medio[j],2))/(j+1));
        }
    }
    
    
    for (int j=0; j<N; j++){
        output2<<vettore2[j]<<" "<<incertezza_sigma[j]<<endl;
    }
    
    output2.close();
    
    //esercizio 1.1.3
    ofstream output3;
    output3.open("ese_3.dat");

    M=100;
    int n=10000;
    
    int misure[M];
    double contatore=0.;
    int campione=100;
    double chi2[campione];
    
    for (int k=0; k<campione; k++){
        
        for(int j=0; j<M; j++){
            misure[j]=0;
        }
        
        for(int i=0; i<n; i++){
            double s=rnd.Rannyu();
            int indice=0;
            contatore=1./M;
            while (contatore<=s) {
                contatore=contatore+(1./M);
                indice=indice+1;
            }
            misure[indice]=misure[indice]+1;
        }
        double chi2_k=0.;
        for(int l=0; l<M; l++){
            chi2_k=chi2_k+pow(misure[l]-(n/M),2);
            
        }
        chi2[k]=(chi2_k/n)*M;
        
        cout<<"chi quadro per la ripetizione "<<k+1<<"-esima è: "<<chi2[k]<<endl;
    }
    
    for(int s=0; s< campione; s++){
        output3<<chi2[s]<<endl;
    }
    output3.close();
    return 0;
}
