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
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include "main.h"

using namespace std;

int main(){
    
    ofstream Histo;
    Histo.open("histogram_16ott");
    InputOptimization();
    Optimization();
    OptimalParameter(step_T);
    x1 = x0;
    walker = 0;
    Input();
    
    for(int iblk=1; iblk <= nblk; ++iblk){
        Reset(iblk);
        for(int istep=1; istep <= nstep; ++istep){
            Move(mu_opt, sigma_opt);
            Measure(mu_opt, sigma_opt);
            Accumulate();
            Histo<<x1<<endl;
        }
        Averages(iblk);
    }
    Histo.close();
    return 0;
}

//-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void Input(void){
    
        ifstream ReadInput;
    
        ReadInput.open("input.dat");
        
        ReadInput >> nblk;
        cout << "Number of blocks = " << nblk << endl;
        
        ReadInput >> nstep;
        cout << "Number of steps in each block = " << nstep << endl;
        
    
        cout << "The program performs Metropolis moves with uniform transaction probability" << endl;
        cout << "Moves parameter = " << delta << endl;
        cout << "Number of blocks = " << nblk << endl;
        cout << "Number of steps in one block = " << nstep << endl << endl;
        ReadInput.close();
        
    }


//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
void InputOptimization(void){
    
    ifstream ReadInput;
    
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    
    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
    
    //Read input informations
    ReadInput.open("opt.dat");
    
    ReadInput >> T;
    cout << "Initial temperature = "<< T <<endl;
    ReadInput >> init_mu;
    cout << "initial Mu = " << init_mu << endl;
    
    ReadInput >> init_sigma;
    cout << "initia Sigma = " << init_sigma << endl;
    

    ReadInput >> step_T;
    ReadInput >> delta_T;
    ReadInput >> step_opt;
    ReadInput >> delta_opt;
    ReadInput >> x0;
    ReadInput >> delta;
    
    
    
    cout<<" "<<endl;
    cout << "The program performs Simulated Annealing in order to calculate optimal parameters mu and sigma" << endl;
    cout << "Moves parameter = " << delta_T << endl;
    cout << "Number of steps per each values of T = " << step_opt << endl << endl;
    ReadInput.close();
    
}


//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
void Optimization(void){
    int pippo=0;
    
    
    ofstream OPTIMIZER;
    double beta;
    double old_mu, old_sigma;
    double new_mu, new_sigma, p;
    double h_old, h_new;
    contatore=0;
    contatore2=0;
    OPTIMIZER.open("parameters_16ott");
    
    old_mu=init_mu;
    old_sigma=init_sigma;
    x1=x0;
    
    h_old=H_mean(old_mu, old_sigma, step_opt);
    
    for(int i=0; i < step_T; i++){
        x1=x0;
        beta=1./(T - i*delta_T);
        contatore=0;
        contatore2=0;
        
        for(int s=0; s < 100; s++){
            x1 = x0;
            h_old=H_mean(old_mu, old_sigma, step_opt);
            
            new_sigma = rnd.Rannyu(old_sigma - delta_opt, old_sigma + delta_opt);
            new_mu = rnd.Rannyu(old_mu - delta_opt, old_mu + delta_opt);
            h_new = H_mean(new_mu, new_sigma, step_opt);
            
            if(h_new > h_old){
                p = rnd.Rannyu();
                if(p <= exp(beta*(h_old - h_new))){
                    old_mu = new_mu;
                    old_sigma = new_sigma;
                    h_old = h_new;
                    contatore += 1;
                }
            }else{
                    old_mu = new_mu;
                    old_sigma = new_sigma;
                    h_old = h_new;
                contatore +=1;
                }
            contatore2+=1;
        }
        OPTIMIZER << old_mu <<" "<< old_sigma <<" "<< h_old << endl;
       
    }
    OPTIMIZER.close();
    
}
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
void OptimalParameter(int dimension){
    ifstream optMuSigma;
    optMuSigma.open("parameters_16ott");
    double mu, sigma, ene;
    optMuSigma>>mu_opt;
    optMuSigma>>sigma_opt;
    optMuSigma>>ene;
    for(int i=1; i<dimension; i++){
        optMuSigma>>mu;
        optMuSigma>>sigma;
        optMuSigma>>ene;
        if(ene< ene_opt){
            mu_opt=mu;
            sigma_opt=sigma;
            ene_opt=ene;
        }
    }
    optMuSigma.close();
    cout<<"_______________________"<<endl<<endl;
    cout<<"mu opt "<<mu_opt<<" sigma opt "<<sigma_opt<<" ene minima "<<ene_opt;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
double H_mean(double mu, double sigma, double step){
    walker=0.;
    for(int j=0; j<step; j++){
        Move(mu, sigma);
        Measure(mu, sigma);
    }
    return walker/(double)step;
    
}
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

void Move(double mu, double sigma){
    double p, p_old, p_new;
    double xold, xnew;
    xold=x1;
    p_old=phi(xold, mu, sigma)*phi(xold, mu, sigma);
    xnew=rnd.Rannyu(xold-delta, xold+delta);
    p_new=phi(xnew, mu, sigma)*phi(xnew, mu, sigma);
    
    p=rnd.Rannyu();
    if(p<=p_new/p_old){
        xold=xnew;
        accepted = accepted + 1.0;
    }
    attempted = attempted + 1.0;
    x1=xold;
    
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

void Measure(double mu, double sigma){
    walker=Energy(x1,mu,sigma);
}
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
double phi(double x, double mu, double sigma){
    double phi;
    
    phi= exp(-pow((x-mu), 2)/(2*sigma*sigma))+exp(-pow((x+mu), 2)/(2*sigma*sigma));
    return phi;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
double Energy(double x1, double mu, double sigma){
    double ene_kin;
    double ene_pot;
    
    ene_kin= 0.5*(exp(-pow((x1-mu),2)/(2*sigma))*(pow(sigma,2)-pow(mu,2)-pow(x1,2)+2*x1*mu)/(pow(sigma,4))
                  + exp(-pow((x1+mu),2)/(2*sigma))*(pow(sigma,2)-pow(mu,2)-pow(x1,2)-2*x1*mu)/(pow(sigma,4)));
    ene_pot= pow(x1,4) - 5./2.*pow(x1,2);
    
    return ene_kin+ene_pot;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

void Reset(int iblk){ //Reset block averages
    
    if(iblk == 1){
        glob_av = 0;
        glob_av2 = 0;
    }
        
    blk_av= 0;
        
    blk_norm = 0;
    attempted = 0;
    accepted = 0;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
void Accumulate(void){
    
    blk_av = blk_av + walker;
    blk_norm = blk_norm + 1.0;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

void Averages(int iblk){
    
    ofstream Hamiltonian;
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Hamiltonian.open("hamiltonian_opt_16ott",ios::app);
    
    stima_ham = blk_av/blk_norm;
    glob_av += stima_ham;
    glob_av2 += stima_ham*stima_ham;
    err_ham=Error(glob_av,glob_av2,iblk);
    
    Hamiltonian << iblk << "," << stima_ham << "," << glob_av/(double)iblk << "," << err_ham << endl;
    
    Hamiltonian.close();
}

//-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_

double Error(double sum, double sum2, int iblk){
    double errore;
    if(iblk==0) {errore=0.;}
    else {errore=sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);}
    return errore;
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
