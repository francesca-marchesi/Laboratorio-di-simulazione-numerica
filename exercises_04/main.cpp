/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
    Input();             //Inizialization
    int nconf = 1;
    //int Nblocchi=50;
    //int Nmisure=nstep/10;
   // int N=Nmisure/Nblocchi;
    int contatore=0;
    
    ofstream WriteE, WriteK, WriteU, WriteT, WriteP, WriteGdR, WriteGave;
    
    double temp_blocco=0.;
    double u_blocco=0.;
    double e_blocco=0.;
    double k_blocco=0.;
    double p_blocco=0.;
    double temp_i=0.;
    double u_i=0.;
    double e_i=0.;
    double k_i=0.;
    double p_i=0.;
    double temp_i2=0.;
    double u_i2=0.;
    double e_i2=0.;
    double k_i2=0.;
    
    double p_i2=0.;
    
    double sigma_temp=0.;
    double sigma_u=0.;
    double sigma_e=0.;
    double sigma_k=0.;
    double sigma_p=0.;
    double gdir=0.;
    double sigma_gdir=0.;
    double r, deltaV;
    
   // WriteT.open("temp_equil_5",ios::app);
    WriteT.open("ave_temp",ios::app);
    WriteE.open("ave_etot",ios::app);
    WriteK.open("ave_ekin",ios::app);
    WriteU.open("ave_epot",ios::app);
    WriteP.open("ave_pressure",ios::app);
    WriteGdR.open("gof", ios::app);
    WriteGave.open("g_ave", ios::app);
    
    
 for(int istep=1; istep <= nstep; ++istep){
        Move();           //Move particles with Verlet algorithm
        if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
        if(istep%10 == 0){
        Measure();     //Properties measurement
       // ConfXYZ(nconf);//Write actual configuration in XYZ format
        
        nconf += 1;
         
        temp_blocco=temp_blocco+stima_temp;
        e_blocco=e_blocco+stima_etot;
        k_blocco=k_blocco+stima_kin;
        p_blocco=p_blocco+stima_press;
        u_blocco=u_blocco+stima_pot;
        for (int k=0; k<nbins; ++k){
            blk_av[k]=blk_av[k] + walker[k];
        }
     
        if(istep%200 == 0){ //in questo ciclo entra 50 volte: cioè il numero i blocchi che vengono fatti
            contatore =contatore+1;
             
            temp_i=temp_i+temp_blocco/20.;
            e_i=e_i+e_blocco/20.;
            k_i=k_i+k_blocco/20.;
            p_i=p_i+p_blocco/20.;
            u_i=u_i+u_blocco/20.;
             
            temp_i2=temp_i2+pow(temp_blocco/20., 2);
            e_i2=e_i2+pow(e_blocco/20., 2);
            k_i2=k_i2+pow(k_blocco/20., 2);
            p_i2=p_i2+pow(p_blocco/20., 2);
            u_i2=u_i2+pow(u_blocco/20., 2);
            
            for(int i=0; i<nbins; ++i){
                r=i*bin_size;
                deltaV=(4./3.)*M_PI*(pow(r+bin_size,3)-pow(r,3));
                gdir=blk_av[i]/(rho*npart*deltaV)/20.;
                glob_av[i] += gdir;
                glob_av2[i] += gdir*gdir;
                sigma_gdir=sqrt((glob_av2[i]/contatore - pow(glob_av[i]/contatore,2))/(double)contatore);
                WriteGdR<< glob_av[i]/(double)contatore<<" "<<sigma_gdir<<endl;
                
            }
            if (contatore!=1){
                sigma_temp=sqrt((temp_i2/contatore-pow(temp_i/contatore,2))/(contatore-1));
                sigma_e=sqrt((e_i2/contatore-pow(e_i/contatore,2))/(contatore-1));
                sigma_k=sqrt((k_i2/contatore-pow(k_i/contatore,2))/(contatore-1));
                sigma_p=sqrt((p_i2/contatore-pow(p_i/contatore,2))/(contatore-1));
                sigma_u=sqrt((u_i2/contatore-pow(u_i/contatore,2))/(contatore-1));
             }
             
             WriteT<<temp_i/contatore<<" "<<sigma_temp<<endl;
             WriteE<<e_i/contatore<<" "<<sigma_e<<endl;
             WriteK<<k_i/contatore<<" "<<sigma_k<<endl;
             WriteP<<p_i/contatore<<" "<<sigma_p<<endl;
             WriteU<<u_i/contatore<<" "<<sigma_u<<endl;
             temp_blocco=0.;
             e_blocco=0.0;
             k_blocco=0.0;
             p_blocco=0.0;
             u_blocco=0.0;
            for (int k=0; k<nbins; ++k){
                blk_av[k]=0.0;
            }
         }
     }
     
  }
    for(int i=0; i<nbins; ++i){
        r=i*bin_size;
        sigma_gdir=sqrt((glob_av2[i]/contatore - pow(glob_av[i]/contatore,2))/(double)contatore);
        WriteGave << r <<  "," << glob_av[i]/(double)contatore << "," << sigma_gdir <<endl;
    }
    WriteGdR.close();
    WriteGave.close();
    WriteT.close();
    WriteE.close();
    WriteK.close();
    WriteU.close();
    WriteP.close();
    ConfFinal();         //Write final configuration to restart
  return 0;
}

//_________________________________________________________________________________________________





void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  ifstream ReadOld;//aggiunto io
  //double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator -> li uso per le velocità
    
  ReadInput.open("input_gas.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> restart;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  if(restart==1){                                         //aggiunto
      cout<< "The program restarts the simulation "<<endl;//aggiunto
  }                                                       //aggiunto
  ReadInput.close(); //----------------------------------------------------

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables
    
    //measurement of g(r)
  nbins = 100;
  bin_size = (box/2.0)/(double)nbins;


//Read initial configuration
 /* cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (int i=0; i<npart; ++i){
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = x[i] * box;
        y[i] = y[i] * box;
        z[i] = z[i] * box;
    }
    ReadConf.close();*/

  
    
//esercizio 4.1 aggiunto

    if(restart==1){
        cout << "Read last configuration from file config.final " << endl << endl;
        ReadConf.open("config.final");
        for (int i=0; i<npart; ++i){
            ReadConf >> x[i] >> y[i] >> z[i];
            x[i] = x[i] * box;
            y[i] = y[i] * box;
            z[i] = z[i] * box;
        }
        ReadConf.close();
        
        cout << "Read old configuration from file old.final " << endl << endl;
        ReadOld.open("old.final");
        for (int i=0; i<npart; ++i){
            ReadOld >> xold[i] >> yold[i] >> zold[i];
            xold[i] = xold[i] * box;
            yold[i] = yold[i] * box;
            zold[i] = zold[i] * box;
        }
        ReadOld.close();
        
        Move();
        
        double t=0.;
        for (int k=0; k<npart; ++k){
            t += 0.5 * (vx[k]*vx[k] + vy[k]*vy[k] + vz[k]*vz[k]);
        }
        
        
        double temperatura = (2.0 / 3.0) * t/(double)npart;
        double alpha=sqrt(temp/temperatura);
        
        cout<<"temperatura impostata: "<<temp<<" temperatura iniziale: "<<temperatura<<endl;
        cout<<"Rescale velocity factor: "<<alpha<<endl;
            
        for(int i=0; i<npart; ++i){
            vx[i] = vx[i]*alpha;
            vy[i] = vy[i]*alpha;
            vz[i] = vz[i]*alpha;
                
            xold[i] = Pbc(x[i]-delta*vx[i]);
            yold[i] = Pbc(y[i]-delta*vy[i]);
            zold[i] = Pbc(z[i]-delta*vz[i]);
            
        } //fine del for
        t=0.0;
        for (int k=0; k<npart; ++k){
            t += 0.5 * (vx[k]*vx[k] + vy[k]*vy[k] + vz[k]*vz[k]);
        }
        cout<< "la temperatura iniziale è: "<<(2.0 / 3.0) * t/(double)npart<< " la temperatura impostata era: "<<temp<<endl;
    } //fine if(restart==1)

    
    if(restart==0){
//Prepare initial velocities
        cout << "Read initial configuration from file config.0 " << endl << endl;
        ReadConf.open("config.0");
        for (int i=0; i<npart; ++i){
            ReadConf >> x[i] >> y[i] >> z[i];
            x[i] = x[i] * box;
            y[i] = y[i] * box;
            z[i] = z[i] * box;
        }
        ReadConf.close();
        
        cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
        double sumv[3] = {0.0, 0.0, 0.0};
        for (int i=0; i<npart; ++i){
            vx[i] = rand() - 0.5;
            vy[i] = rand() - 0.5;
            vz[i] = rand() - 0.5;

            sumv[0] += vx[i];
            sumv[1] += vy[i];
            sumv[2] += vz[i];
        }
        for (int idim=0; idim<3; ++idim){
            sumv[idim] /= (double)npart;
        }
    
        double sumv2 = 0.0, fs;
    
        for (int i=0; i<npart; ++i){
            vx[i] = vx[i] - sumv[0];
            vy[i] = vy[i] - sumv[1];
            vz[i] = vz[i] - sumv[2];

            sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        }
        sumv2 /= (double)npart;

        fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
        for (int i=0; i<npart; ++i){
            vx[i] *= fs;
            vy[i] *= fs;
            vz[i] *= fs;

            xold[i] = x[i] - vx[i] * delta;
            yold[i] = y[i] - vy[i] * delta;
            zold[i] = z[i] - vz[i] * delta;
        }
    } //fine di if(restart==0)
   return;
}
//-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_ fine di Input

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}//-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_fine di Move

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}//-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_fine di Force

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij, pij;
  double p=0.;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Press;

  Epot.open("output_epot",ios::app);
  Ekin.open("output_ekin",ios::app);
  Temp.open("output_temp",ios::app);
  Etot.open("output_etot",ios::app);
  Press.open("output_press",ios::app);

  v = 0.0; //reset observables
  t = 0.0;

  for (int k=0; k<nbins; ++k) walker[k]=0.0; //reset histogram


//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
    //update of the histogram of g(r)
     if(dr<=box/2.){
            bin=dr/bin_size;
            walker[bin] = walker[bin]+1;
     }

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       pij=48.0/pow(dr,12)-24.0/pow(dr,6);

//Potential energy
       v += vij;
//pressure
       p += pij;
         
     }
    }          
  }

//Kinetic energy
    for (int i=0; i<npart; ++i){
        t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    }
   
    stima_pot = v/(double)npart; //Potential energy tutte per unità di particelle
    stima_kin = t/(double)npart; //Kinetic energy
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total enery
    stima_press= rho*stima_temp+(p*rho)/(3.0*(double)npart);

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Press << stima_press<<endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Press.close();
    return;
}//-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_fine di Measure

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
    
  double xvecchio, yvecchio, zvecchio;
    
  ofstream WriteOld; //parte nuova
      cout << "salvo le penultime posizioni per ripartire "<<endl;
    WriteOld.open("old.final");
    for (int i=0; i<npart; ++i){
        xvecchio=x[i]-vx[i]*(2*delta);
        yvecchio= y[i]-vy[i]*(2*delta);
        zvecchio= z[i]-vz[i]*(2*delta);
        WriteOld << xvecchio/box << "   " <<  yvecchio/box << "   " << zvecchio/box << endl;
    }
  WriteOld.close();

  return;
} //-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_fine di ConfFinal

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;
  
  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
    
}//-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_fine di ConfXYZ

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}//-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_fine di PBC



/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
