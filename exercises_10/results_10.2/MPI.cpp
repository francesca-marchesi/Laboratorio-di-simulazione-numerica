#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"
#include "MPI_simulatedannealing.h"
#include <vector>

#include <algorithm>    // std::random_shuffle
#include <ctime>        // std::time
#include <cstdlib>
#include <numeric>      // std::accumulate
#include <stdlib.h>
#include "mpi.h"

using namespace std;

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &my_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    //Circle(N, 2.0);
    Square(N, 4.0);
    
    Input();
    vector <int> new_vector;
    double r;
    double lenght;
    for(int i=0; i<N; i++){
        old_vector.push_back(i+1);
    }
    
    random_shuffle(old_vector.begin(),old_vector.end());
    
    for(int it=0; it<iterations; it++){
        T=T_init-it*delta_T;
        beta=1./T;
        for(int j=0; j<step; j++){
            new_vector=old_vector;
            r=rnd.Rannyu();
            if (r<=1./4.){
                new_vector=Mutazione1(probability, new_vector);
            }
            if (r<=2./4. && r>1./4.){
                new_vector=Mutazione2(probability, new_vector);
            }
            if (r<=3./4. && r>2./4.){
                new_vector=Mutazione4(probability, new_vector);
            }
            if (r<=4./4. && r>3./4.){
                new_vector=Mutazione5(probability, new_vector);
            }
            
            Move(new_vector, N, beta);
            WriteLenght(N,old_vector);
            
        }
    }
    lenght=Fitness(N, old_vector);
    double *lenghts = NULL;
    if (my_rank==0){
        lenghts = (double *)malloc(sizeof(double)*my_size);
    }
    
    MPI_Gather(&lenght, my_size, MPI_FLOAT, lenghts, my_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    int best_rank=0;
    if(my_rank==0){
        cout<<"pippo"<<endl;
        double index=lenghts[0];
        for(int i=1; i<my_size; i++){
            if(index>lenghts[i]){
                index=lenghts[i];
                best_rank=i;
            }
        }
    }
    if(my_rank==best_rank){
        FinalConf(old_vector);
    }
    
    MPI_Finalize();
    return 0;
}


//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
double Fitness(int N, vector<int> vec){
    double L=0.;
    
    for(int i=0; i<N-1; i++){
        L=L+ sqrt(pow(x[vec[i]-1]-x[vec[i+1]-1],2)+pow(y[vec[i]-1]-y[vec[i+1]-1],2));
    }
    
    L=L+sqrt(pow(x[vec[0]-1]-x[vec[N-1]-1],2)+pow(y[vec[0]-1]-y[vec[N-1]-1],2));
    
    return L;
}
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-


void Move(vector<int> vec2, int N, double beta){
    double p, energy_old, energy_new, o;
    energy_old=Boltzmann(N, old_vector); //lunghezza della sequenza vecchia
    energy_new=Boltzmann(N, vec2); //lunghezza della sequenza nuova
    if(energy_old>=energy_new){
        for(int i=0; i<N; i++){
            old_vector[i]=vec2[i];
        }
    }else{
        p=exp(-beta*(energy_new-energy_old));
        o=rnd.Rannyu();
        if(o<=p){
            old_vector=vec2;
        }
    }
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

double Boltzmann(int N, vector <int> vec){
    double ene = Fitness(N, vec);
    return ene;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

 std::vector<int> Mutazione1(double probability, vector<int> vec){
     double mut1=rnd.Rannyu();
     if(mut1<=probability){
         int k;
         int t;
         k=ceil(rnd.Rannyu(-1,vec.size()-1));
         t=ceil(rnd.Rannyu(-1, vec.size()-1));
         swap (vec[k], vec[t]);
     }
     return vec;
 }

//------------------------
std::vector<int> Mutazione2(double probability, vector<int> vec ){
    double mut1=rnd.Rannyu();
    int shift;
    if (mut1<=probability){ //shift fisso di 3
        shift=ceil(rnd.Rannyu(-1, vec.size()-1));
        rotate(vec.begin(),vec.begin()+shift,vec.end());
    }
    return vec;
}

 //-------------------------
 std::vector<int> Mutazione4(double probability, vector<int> vec){
 double mut4=rnd.Rannyu();
 if(mut4<=probability){
 
    int dimension=ceil(rnd.Rannyu(0,vec.size()/2.));
    int first=ceil(rnd.Rannyu(-1, vec.size()-2*dimension));
    int second=ceil(rnd.Rannyu(first+dimension, vec.size()-dimension));
    vector <int> vettore;
    for(int i=0; i<dimension; i++){
        vettore.push_back(vec[first+i]);
    }
    for(int i=0; i<dimension; i++){
        vettore.push_back(vec[second+i]);
    }
    random_shuffle (vettore.begin(), vettore.end());
    
    for(int i=0; i<dimension; i++){
        vec[first+i]=vettore[i];
    }
 
    for(int i=0; i<dimension; i++){
        vec[second+i]=vettore[dimension+i];
    }
 }
 return vec;
 }
//------------------------

std::vector<int> Mutazione5(double probability, vector<int> vec){
    double mut5=rnd.Rannyu();
   
    if (mut5<=probability){
        int end=ceil(rnd.Rannyu(-1, vec.size()));
        int start=ceil(rnd.Rannyu(-1, end));
 
        reverse(vec.begin()+start,vec.begin()+end);
    }
    return vec;
}
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
void WriteLenght(int N, vector<int> vec){
    ofstream WriteL;
    ofstream WriteSequence;
    WriteL.open("2lunghezza_C_"+to_string(my_rank), ios::app);
    WriteL<<T<<"    "<<Fitness(N, vec)<<endl;
    WriteL.close();
    
    for(int i=0; i<N; i++){
        WriteSequence.open("2sequenza_C_"+to_string(my_rank), ios::app);
        WriteSequence<<old_vector[i]<<",";
        WriteSequence.close();
    }
    WriteSequence.open("2sequenza_C_"+to_string(my_rank), ios::app);
    WriteSequence<<endl;
    WriteSequence.close();
    
    
}
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
void FinalConf(vector<int> vec){
    ofstream best_configuration;
    ofstream best_path;
    int indice=0;
    for(int s=0; s<vec.size(); s++){
        best_configuration.open("2SequenzaFinale_C"+to_string(my_rank), ios::app);
        indice=vec[s]-1;
        best_configuration<<x[indice]<<" "<<y[indice]<<endl;
        best_configuration.close();
    }
    
    best_path.open("2best_path_C");
    for(int i=0; i<N; i++){
        best_path<<old_vector[i]<<",";
    }
    best_path.close();
}
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
void Check(int N, vector<vector<int>> vec){
    int somma=N*(N+1)/2;
    int somma2;
    for(int d=0; d<vec.size(); d++){
        
        somma2= accumulate(vec[d].begin(),vec[d].end(),0);
        
        if(somma2!=somma){
            cout<<"la sequenza "<<d<<"-esima non rispetta i vincoli"<<endl;
            exit(0);
        }
    }
    
}
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
void Crossover(double probability, vector<int> vec1, vector<int> vec2){
    double cross=rnd.Rannyu();
    int a[2], b[2];
    int X, Y, Z, W;
    if (cross<=probability){
        X=vec1[vec1.size()-2];
        Y=vec1[vec1.size()-1];
        Z=vec2[vec2.size()-2];
        W=vec2[vec2.size()-1];
        
        for(int j=0; j<vec2.size(); j++){
            if(vec2[j]==X){
                a[0]=j;
            }
            if(vec2[j]==Y){
                a[1]=j;
            }
            if(vec1[j]==Z){
                b[0]=j;
            }
            if(vec1[j]==W){
                b[1]=j;
            }
        }
        
        if(a[0]>a[1]){
            swap(vec1[vec1.size()-2], vec1[vec1.size()-1]);
        }
        if(b[0]>b[1]){
            swap(vec2[vec2.size()-2], vec2[vec2.size()-1]);
        }
        
    }
}


//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
void Input(void){
    ifstream ReadInput;
    if (my_rank==0){
        cout << "Simulated annealing for Traveling Salesman Problem" <<endl;
    }
    
    int p1[my_size], p2[my_size];
    ifstream Primes("Primes");
    
    for (int i=0; i<my_size; i++){
        Primes >> p1[i] >> p2[i];       //ho una coppia di numeri per nodo
    }
    Primes.close();
    
    
    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1[my_rank],p2[my_rank]);
    input.close();
    
    
    ReadInput.open("input.dat");
    
    ReadInput >> N;
    //cout << "Il numero di città considerate è = " <<  N << endl;
    
    ReadInput >> T_init;
    //cout << "Temperatura iniziale = " <<  T_init << endl;
    //cout <<"Beta = "<<1./T_init <<endl;
    
    ReadInput >> delta_T;
    //cout << "Delta  = " << delta_T << endl;
    
    ReadInput >> step;
    //cout << "Numero di step in ogni iterazione = " << step << endl;
    
    ReadInput >> iterations;
    //cout << "Numero di iterazioni totali = " << iterations << endl;
   
    ReadInput >> probability;
    //cout << "Probabilità delle singole iterazioni = " << probability << endl;
    
    ReadInput.close();
    
}
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

void Square(int n_cities, double l){
    ifstream CoordS;
    CoordS.open("Coordinate_Quadrato");
    for(int i=0; i<n_cities; i++){
        CoordS>>x[i];
        CoordS>>y[i];
    }
    CoordS.close();
}


//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

void Circle(int n_cities, double R){
    ifstream CoordC;
    CoordC.open("Coordinate_Cerchio");
    for(int i=0; i<n_cities; i++){
        CoordC>>x[i];
        CoordC>>y[i];
    }
    CoordC.close();
}


