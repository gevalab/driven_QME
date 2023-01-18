#ifndef CONSTANT_H
#define CONSTANT_H
#include <vector>
#include <complex>
#include <armadillo>
using namespace std;
using namespace arma;

const double h_ = 1.0;         // reduced plank constant, in atomic unit
const double pi = 3.14159265358979323846; //pi

class Hamiltonian
{
    public:
    static double  w;          //frequency of the field
    static double  w01;        //frequency of the 2-level system
    static double  epsilon;    //system-field coupling
    static double  Lambda;     //system-bath coupling
    static int     N;          //number of bath modes
    static double  wc;         //cutoff frequency of Ohmic spectral density
    static double  eta;        //coupling strength of Ohmic spectral density 
    vec w_;         //frequency of the bath modes
    vec c;          //coupling constant of the bath modes   
    vec Req;        //equilibrium position of the bath modes
    void print() const;
    int Ohmic_mode_generator();
};
double  Hamiltonian::w       = 0.0;
double  Hamiltonian::w01     = 0.0;
double  Hamiltonian::epsilon = 0.0;
double  Hamiltonian::Lambda  = 0.0;
int     Hamiltonian::N       = 1;
double  Hamiltonian::wc      = 1.0;
double  Hamiltonian::eta     = 1.0;
void Hamiltonian::print() const{
    cout << "w          " << w          << endl;
    cout << "w01        " << w01        << endl;
    cout << "epsilon    " << epsilon    << endl;
    cout << "Lambda     " << Lambda     << endl;
    cout << "N          " << N          << endl;
    cout << "wc         " << wc         << endl;
    cout << "eta        " << eta        << endl;
    //cout << "w_         " << endl;
    //cout << w_ << endl;
    //cout << "Req        " << endl;
    //cout << Req << endl;
}
int Hamiltonian::Ohmic_mode_generator(){
    w_  = vec(N,fill::zeros);
    c   = vec(N,fill::zeros);
    Req = vec(N,fill::zeros);
    for(int i=0;i<N;i++){
        w_(i)  = wc*log(N/(N-0.5-i));
        c(i)   = sqrt(2*eta*wc/(pi*N))*w_(i);
        Req(i) = 2*c(i)/pow(w_(i),2);
    }
    w_.print("w_:");
    c.print("c:");
    Req.print("Req:");
    return 0;
}

class timer
{
    public:
    double T;            //total simulation time
    double dt;           //timestep
    double t;            //current time, default=0
    timer(double=1,double=0.1,double=0);
    void print() const;
};
timer::timer(double T_set, double dt_set, double t_set){
    T = T_set;
    dt = dt_set;
    t = t_set;
}
void timer::print() const{
    cout << "t/T =  " << t << "/" << T << endl;
}

class state
{
    public:
    cx_mat sigma = cx_mat(2,2,fill::zeros);  //density matrix
};

#endif