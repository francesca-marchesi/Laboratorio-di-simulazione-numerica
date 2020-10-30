/****************************************************************
 *****************************************************************
 _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
 _/_/  _/ _/       _/       Physics Department
 _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
 _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "main.h"

using namespace std;

int main (int argc, char *argv[]){
    Input();
    int Ntrial=Nsim/Nblocchi;
    ofstream WritePosition, WriteR;
   
    
    
    //GROUND STATE
    double r0;
    double r_proposto;
    double ri=0.;
    double ri2=0.;
    double rblocco=0.;
    double r2blocco=0.;
    double sigma=0.;
    double rapporto;
    int contatore=0;
    
    WriteR.open("R_gs");
    WritePosition.open("PositionsXYZ_gs");
    WritePosition<<x_0<<" "<<y_0<<" "<<z_0<<" "<<sqrt(x_0*x_0+y_0*y_0+z_0*z_0)<<endl;
    
    r0=sqrt(x_0*x_0+y_0*y_0+z_0*z_0);
    
    for(int i=0; i<100; i++){
        x=rnd.Rannyu(x_0-costGS,x_0+costGS);
        y=rnd.Rannyu(y_0-costGS,y_0+costGS);
        z=rnd.Rannyu(z_0-costGS,z_0+costGS);
        r_proposto=sqrt(x*x+y*y+z*z);
        rapporto=GS(r0, r_proposto);
        r0=MoveGS(r0, r_proposto, rapporto, contatore);
    }
    
    
    for(int blocks=0; blocks<Nblocchi; blocks++){
        ri=0.;
        ri2=0.;
        for(int trials=0; trials<Ntrial; trials++){ //blocco-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
            
            x=rnd.Rannyu(x_0-costGS,x_0+costGS); //IN COORDINATE CARTESIANE
            y=rnd.Rannyu(y_0-costGS,y_0+costGS);
            z=rnd.Rannyu(z_0-costGS,z_0+costGS);
    
            r_proposto=sqrt(x*x+y*y+z*z);
            rapporto=GS(r0, r_proposto);
            
           // cout<<"rapporto: "<<rapporto<<endl;
            r0=MoveGS(r0, r_proposto, rapporto, contatore);
            
            ri=ri+r0;
            ri2=ri2+r0*r0;
            if(blocks!=0 && blocks!=1){
                WritePosition<<x_0<<" "<<y_0<<" "<<z_0<<" "<<sqrt(x_0*x_0+y_0*y_0+z_0*z_0)<<endl;
            }
            } //-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
        
        rblocco=rblocco+ri/(double)Ntrial;
        r2blocco=r2blocco+ri2/(double)Ntrial;
        if(blocks==0){sigma=0.;}
        if(blocks!=0){
            sigma=sqrt((r2blocco/(double)(blocks+1)-pow(rblocco/(double)(blocks+1), 2))/(blocks));
        }
            WriteR<<rblocco/(double)(blocks+1)<<" "<<sigma<<endl;
        
    }
    cout<<"stima GS "<<rblocco/(double)(Nblocchi)<<" valore vero "<<3./2.<<endl;
    
    WritePosition.close();
    WriteR.close();
    
    
     //2P EXCITED STATE
    ofstream WritePosition2, WriteR2;
    double R0;
    double R_proposto;
    double Ri=0.;
    double Ri2=0.;
    double Rblocco=0.;
    double R2blocco=0.;
    double Sigma=0.;
    double Rapporto;
    int contatore2=0;
    
    WritePosition2.open("PositionsXYZ_2p");
    WriteR2.open("R_2p");
    
    R0=sqrt(X_0*X_0+Y_0*Y_0+Z_0*Z_0);
    
    
    for(int i=0; i<100; i++){
        X=rnd.Rannyu(X_0-cost2p,X_0+cost2p);
        Y=rnd.Rannyu(Y_0-cost2p,Y_0+cost2p);
        Z=rnd.Rannyu(Z_0-cost2p,Z_0+cost2p);
        R_proposto=sqrt(X*X+Y*Y+Z*Z);
        Rapporto=EXCITED(R0, R_proposto, Z, Z_0);
        R0=MoveEX(R0, R_proposto, Rapporto, contatore2);
    }
    for(int blocks=0; blocks<Nblocchi; blocks++){//--------
        Ri=0.;
        Ri2=0.;
        for(int trials=0; trials<Ntrial; trials++){
            
            X=rnd.Rannyu(X_0-cost2p,X_0+cost2p);
            Y=rnd.Rannyu(Y_0-cost2p,Y_0+cost2p);
            Z=rnd.Rannyu(Z_0-cost2p,Z_0+cost2p);
            
            R_proposto=sqrt(X*X+Y*Y+Z*Z);
            Rapporto=EXCITED(R0, R_proposto, Z, Z_0);
            
            R0=MoveEX(R0, R_proposto, Rapporto, contatore2);
            
            Ri=Ri+R0;
            Ri2=Ri2+R0*R0;
            if(blocks!=0 && blocks!=1){
                WritePosition2<<X_0<<" "<<Y_0<<" "<<Z_0<<" "<<sqrt(X_0*X_0+Y_0*Y_0+Z_0*Z_0)<<endl;
            }
        }
        Rblocco=Rblocco+Ri/(double)Ntrial;
        R2blocco=R2blocco+Ri2/(double)Ntrial;
        if(blocks==0){sigma=0.;}
        if(blocks!=0){
            Sigma=sqrt((R2blocco/(double)(blocks+1)-pow(Rblocco/(double)(blocks+1), 2))/(blocks));
        }
            WriteR2<<Rblocco/(double)(blocks+1)<<" "<<Sigma<<endl;
    }//-------------------
    cout<<"stima 2p "<<Rblocco/(double)(Nblocchi)<<" valore vero "<<5<<endl;
    WriteR2.close();
    WritePosition2.close();
    return 0;
    
}

//-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

void Input(void){
    ifstream ReadInput;
    
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    
    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
    
    ReadInput.open("input.dat");
    
    ReadInput >> x_0;
    ReadInput >> y_0;
    ReadInput >> z_0;
    cout << "La posizione iniziale è per il GS = (" <<x_0 <<" , " <<y_0<<" , "<<z_0<<")"<<endl;
    ReadInput >> X_0;
    ReadInput >> Y_0;
    ReadInput >> Z_0;
    cout << "La posizione iniziale è per lo stato 2p = (" <<X_0 <<" , " <<Y_0<<" , "<<Z_0<<")"<<endl;
    
    
    ReadInput >> Nsim;
    cout << "Il numero di simulazioni totali è = " << Nsim << endl;
    
    ReadInput >> Nblocchi;
    cout << "Il numero di blocchi fatti è = " << Nblocchi << endl;
    
    ReadInput >> costGS;
    cout << "La larghezza del passo di metropolis è = " << 2*costGS << endl;
    
    ReadInput >> cost2p;
    cout << "La larghezza del passo di metropolis è = " << 2*cost2p << endl;
    cout<<endl;
    cout<<endl;
    
    
    ReadInput.close();
    
    
    
}


//____________________________________________________________________________
double MoveGS(double r0, double r_proposto, double rapporto, int contatore){
    double A;
    double s;
    A=min(1., rapporto);
    s=rnd.Rannyu();
    if(A>=s){
        x_0=x;
        y_0=y;
        z_0=z;
        contatore=contatore+1.;
        r0=r_proposto;
        
    }
    return r0;
}

//____________________________________________________________________________
double MoveEX(double R0, double R_proposto, double rapporto, int contatore){
    double A;
    double S;
    A=min(1., rapporto);
    S=rnd.Rannyu();
    if(A>=S){
        X_0=X;
        Y_0=Y;
        Z_0=Z;
        contatore=contatore+1.;
        R0=R_proposto;
    }
    return R0;
}
//____________________________________________________________________________
double GS(double r0, double r_proposto){
    double rapporto=exp(-2.0*(r_proposto-r0));

    return rapporto;
}

//____________________________________________________________________________

double EXCITED(double raggio, double raggio_proposto, double zeta, double zeta_0){
    double rapporto=pow(zeta/zeta_0,2)*exp(raggio-raggio_proposto);

    return rapporto;
}



/****************************************************************
 *****************************************************************
 _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
 _/_/  _/ _/       _/       Physics Department
 _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
 _/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
 *****************************************************************
 *****************************************************************/
