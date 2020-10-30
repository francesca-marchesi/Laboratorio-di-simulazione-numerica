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
#include "genetic.h"
#include <vector>

#include <algorithm>    // std::random_shuffle
#include <ctime>        // std::time
#include <cstdlib>
#include <numeric>      // std::accumulate
#include <stdlib.h>



using namespace std;

int main(int argc, const char * argv[]) {
    
    ofstream Seq;
    double l=0.;
    int seq=0;
    
    Input();
    
    //Circle(N, 2.0);
    Square(N, 4.0);
    Generatrice(n, N);

    for(int generation=0; generation<N_generations; generation++){
        Check(N, sample);
        Selection(n, N);
        WriteLenght(N_generations, dist);
        dist.clear();
        
        for(int o=0; o<n/2; o++){
            Extract(n);
            v1=sample[c[0]];
            v2=sample[c[1]];
            
            v1=Mutazione1(.1, v1);
            v2=Mutazione1(.1, v2);
            v1=Mutazione2(.1, v1);
            v2=Mutazione2(.1, v2);
            v1=Mutazione4(.1, v1);
            v2=Mutazione4(.1, v2);
            v1=Mutazione5(.1, v1);
            v2=Mutazione5(.1, v2);
            Crossover(.7, v1, v2);
            sons.push_back(v1);
            sons.push_back(v2);
        }
        sample=sons;
        sons.clear();
        prob.clear();
        
    }
    l=Fitness(N, 0);
    seq=0;
    for(int f=1; f<n; f++){
        if(Fitness(N, f)<l){
            l=Fitness(N, f);
            seq=f;
        }
    }
    FinalConfXY(seq);
    return 0;
}


//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
double Fitness(int N, int o){
    double L=0.;
    
    for(int i=0; i<N-1; i++){
        L=L+ sqrt(pow(x[sample[o][i]-1]-x[sample[o][i+1]-1],2)+pow(y[sample[o][i]-1]-y[sample[o][i+1]-1],2));
    }
   
    L=L+sqrt(pow(x[sample[o][0]-1]-x[sample[o][N-1]-1],2)+pow(y[sample[o][0]-1]-y[sample[o][N-1]-1],2));
    return L;
}
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
void Selection(int n, int N){
    double lunghezza=0.;
    double sum=0.;
    for(int j=0; j<n; j++){
        lunghezza=Fitness(N, j);
        sum=sum+1./lunghezza;
        dist.push_back(lunghezza);
    }
    for(int k=0; k<n; k++){
        prob.push_back(1./dist[k]/sum);
    }
}
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
void Extract(int n){
    double contatore1=0.;
    double contatore2=0.;
    double s1=rnd.Rannyu();
    double s2=rnd.Rannyu();
    
    for(int l=0; l<n; l++){
        contatore1=contatore1+prob[l];
        if(s1<=contatore1){
            c[0]=l;
            break;}
    }
    for(int l2=0; l2<n; l2++){
        contatore2=contatore2+prob[l2];
        if(s2<=contatore2 && l2!=c[0]){
            c[1]=l2;
            break;}
    }

}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

 void Crossover(double probabilità, vector<int> vec1, vector<int> vec2){
    double cross=rnd.Rannyu();
    int a[2], b[2];
    int X, Y, Z, W;
    if (cross<=probabilità){
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
 //cout<<"crossover avvenuto con successo"<<endl;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

std::vector<int> Change(int a, int b, int check, vector<int> vec){
    if (check==1){
        int x;
        x=vec[a];
        vec[a]=vec[b];
        vec[b]=x;
    }
    
    return vec;
}
 
 
 //_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
 
 std::vector<int> Mutazione1(double probabilità, vector<int> vec){
     double mut1=rnd.Rannyu();
     if(mut1<=probabilità){
         int k;
         int t;
         k=ceil(rnd.Rannyu(-1,vec.size()-1));
         t=ceil(rnd.Rannyu(-1, vec.size()-1));
         swap (vec[k], vec[t]);
 
     }
 
 
 //cout<<"prima mutazione avvenuta con successo"<<endl;
 return vec;
 }
 
 

 
 
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

std::vector<int> Mutazione2(double probabilità, vector<int> vec ){
    double mut1=rnd.Rannyu();
    int shift;
    if (mut1<=probabilità){ //shift fisso di 3
    shift=ceil(rnd.Rannyu(-1, vec.size()-1));
    rotate(vec.begin(),vec.begin()+shift,vec.end());

 
 }
 
// cout<<"second mutazione avvenuta con successo"<<endl;
    return vec;
}
 //_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
 std::vector<int> Mutazione4(double probabilità, vector<int> vec){
 double mut4=rnd.Rannyu();
 if(mut4<=probabilità){
 
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
 
// cout<<"quarta mutazione avvenuta con successo"<<endl;
 
 return vec;
 }
 
 
 

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-


std::vector<int> Mutazione5(double probabilità, vector<int> vec){
    double mut5=rnd.Rannyu();
   
    if (mut5<=probabilità){
 int end=ceil(rnd.Rannyu(-1, vec.size()));
 int start=ceil(rnd.Rannyu(-1, end));
 
 reverse(vec.begin()+start,vec.begin()+end);
 
    }
 
// cout<<"quinta mutazione avvenuta con successo"<<endl;
          return vec;
}
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
void Generatrice(int n, int N){
    vector<int> vector1;
    for (int j = 0; j < N; j++) {
        vector1.push_back(j+1);
    }
 
    for (int i=0; i<n; i++){ //n è il numero di righe= numero di elementi nel sample che genero
        random_shuffle (vector1.begin(), vector1.end());
        sample.push_back(vector1);
    }
    /*cout<<"rows "<<sample.size()<<endl;
    cout<<"columns "<<sample[0].size()<<endl;
    for (int i = 0; i < sample.size(); i++)       // loops through each row of vy
    {  for (int j = 0; j < sample[i].size(); j++) // loops through each element of each row
        cout << " " << sample[i][j];           // prints the jth element of the ith row
        cout << endl;
    }*/
}


//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

void Square(int n_cities, double l){
    ofstream CoordS;
    CoordS.open("results/CoordinateQuadrato");
    
    for(int i=0; i<n_cities; i++){
        x[i]=rnd.Rannyu(-l/2., l/2.);
        y[i]=rnd.Rannyu(-l/2., l/2.);
        CoordS<<x[i]<<" "<<y[i]<<endl;
    }
    CoordS.close();
}


//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

void Circle(int n_cities, double R){
    ofstream CoordC;
    CoordC.open("results/CoordinateCerchio");
    
    double tetha;
    for(int i=0; i<N; i++){
        tetha=rnd.Rannyu(0, 2*M_PI);
        x[i]=R*cos(tetha);
        y[i]=R*sin(tetha);
        
        CoordC<<x[i]<<" "<<y[i]<<endl;
    }
    CoordC.close();
}
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
void WriteLenght(int generation, vector<double> dist){
    ofstream WriteMean, WriteL;
    WriteMean.open("results/Lunghezza_media_quadrato", ios::app);
    WriteL.open("results/Lunghezza_minima_quadrato", ios::app);
    double l_mean=0.;
    int max=dist.size()/2;
    
    sort(dist.begin(), dist.end());
    for(int m=0; m<max; m++){
        l_mean=l_mean+dist[m];
    }
    WriteMean<<generation<<" "<<l_mean/double(max)<<endl;
    WriteL<<generation<<" "<<dist[0]<<endl;
    WriteMean.close();
    WriteL.close();
    
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
void Input(void)
{
    ifstream ReadInput;
    
    cout << "Genetic Algotithm for Traveling Salesman Problem" << endl;
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    
    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
    
    
    ReadInput.open("input.dat");
    
    ReadInput >> N;
    cout << "Numero di città = " << N << endl;
    
    ReadInput >> N_generations;
    cout << "Numero di generazioni = " << N_generations << endl;
    
    ReadInput >> n;
    cout << "Numero di elementi in ogni generazione = " << n << endl;
    ReadInput.close();
    
}
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
void FinalConfXY(int best){
    ofstream WriteConf;
    WriteConf.open("results/configurazione_finale_quadrato");
    int indice=0;
    
    for(int s=0; s<sample[best].size(); s++){
        indice=sample[best][s];
        WriteConf<<x[indice-1]<<" "<<y[indice-1]<<endl;
    }
    WriteConf.close();
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














