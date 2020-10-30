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
    
    for(int i=0; i<20; i++){
    }
    
    rnd.SaveSeed();
    
    //random walk discreto su un reticolo cubico 3D-----------------------------
    ofstream output;
    output.open("02.2.1.txt");
    int N= 10000;
    int N_blocchi=100;
    int n=N/N_blocchi;
    int M=100;
    double a=1.;
    double lunghezza;
    double b=1/3.;
    
    double x=0.;
    double y=0.;
    double z=0.;
    double xi=0.;
    double yi=0.;
    double zi=0.;
    double dist_quad=0.;
    double dist_quarta=0.;
    double v[M];
    double w[M];
    double mu_dist[M];
    double incertezza_stat[M];
    double A[M];
    double A2[M];
    double mu_rad[M];
    v[0]=0.;
    w[0]=0.;
    for(int k=0; k<M; k++){
        A[k]=0.;
        A2[k]=0.;
        mu_dist[k]=0.;
        incertezza_stat[k]=0.;
        mu_rad[k]=0.;
    }
    double r;
    double s;
    for (int h=0; h<N_blocchi; h++){
        for(int k=0; k<M; k++){
            v[k]=0.;
            w[k]=0.;
        }
        
        for(int i=0; i< n; i++){
            x=0.;
            y=0.;
            z=0.;
            for(int j=0; j<M; j++){
                r=rnd.Rannyu();
                s=rnd.Rannyu();
                if(r<0.5){
                    lunghezza=-a;
                }
                else{
                    lunghezza=a;
                }
                if(s<b){
                    xi=lunghezza;
                    yi=0.;
                    zi=0.;
                }
                else if (s>=b && s<2*b){
                    xi=0.;
                    yi=lunghezza;
                    zi=0.;
                }
                else{
                    xi=0.;
                    yi=0.;
                    zi=lunghezza;
                }
                x=x+xi;
                y=y+yi;
                z=z+zi;
                dist_quad=x*x+y*y+z*z;
                dist_quarta=dist_quad*dist_quad;
            
                v[j+1]=v[j+1]+dist_quad;
                w[j+1]=w[j+1]+dist_quarta;
            }
        }
        
        for(int l=1; l<M; l++){
            mu_dist[l]=mu_dist[l]+(v[l]/n);
            A[l]=A[l]+v[l]/n;
            A2[l]=A2[l]+w[l]/n;
        }
    }
        
        for(int q=0; q<M; q++){
            mu_rad[q]=sqrt(mu_dist[q]/N_blocchi);
            incertezza_stat[q]=sqrt((A2[q]/N_blocchi-pow((A[q]/N_blocchi), 2))/(N_blocchi-1));
            output<<mu_rad[q]<< ","<<incertezza_stat[q]<<endl;
        }
    
    
    output.close();
    
    
    //continuo-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
    
    ofstream output2;
    output2.open("02.2.2.txt");
    //N= 10000 numero di simulazioni
    //N_blocchi=100 numero di blocchi
    //n=N/N_blocchi

    double theta;
    double phi;
    double xpos=0.;
    double ypos=0.;
    double zpos=0.;
    
    double xpos_i=0.;
    double ypos_i=0.;
    double zpos_i=0.;
    
    double distance_quad=0.;
    double distance_quarta=0.;
    
    double vpos[M];
    double vpos_quarta[M];
    
    double distanza_media[M];
    double A_continuo[M];
    double A2_continuo[M];
    double radice_dist[M];
    double incertezza_stat_cont[M];
    
    for(int k=0; k<M; k++){
        A_continuo[k]=0.;
        A2_continuo[k]=0.;
        distanza_media[k]=0.;
        incertezza_stat_cont[k]=0.;
        radice_dist[k]=0.;
    }
    
    
    
    for (int h=0; h<N_blocchi; h++){
        for(int k=0; k<M; k++){
            vpos[k]=0.;
            vpos_quarta[k]=0.;
        }

        for(int K=0; K<n; K++){
            xpos=0.;
            ypos=0.;
            zpos=0.;
        
            for(int I=0; I<M; I++){
                theta=rnd.Rannyu(0, M_PI);
                phi=rnd.Rannyu(0, 2*M_PI);
                xpos_i=a*sin(theta)*cos(phi);
                ypos_i=a*sin(theta)*sin(phi);
                zpos_i=a*cos(theta);
            
                xpos=xpos+xpos_i;
                ypos=ypos+ypos_i;
                zpos=zpos+zpos_i;
                distance_quad=xpos*xpos+ypos*ypos+zpos*zpos;
                distance_quarta=distance_quad*distance_quad;
            
                vpos[I+1]=distance_quad+vpos[I+1];
                vpos_quarta[I+1]=vpos_quarta[I+1]+distance_quarta;
            }
        }
        
        for(int l=1; l<M; l++){
            distanza_media[l]=distanza_media[l]+(vpos[l]/n);
            A_continuo[l]=A_continuo[l]+vpos[l]/n;
            A2_continuo[l]=A2_continuo[l]+vpos_quarta[l]/n;
        }
        
    }

    for(int q=0; q<M; q++){
        radice_dist[q]=sqrt(distanza_media[q]/N_blocchi);
        incertezza_stat_cont[q]=sqrt((A2_continuo[q]/N_blocchi-pow((A_continuo[q]/N_blocchi), 2))/(N_blocchi-1));
        output2<<radice_dist[q]<< ","<<incertezza_stat_cont[q]<<endl;
    }

    output2.close();

    return 0;
}

























