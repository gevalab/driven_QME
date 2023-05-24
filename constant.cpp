#include "constant.h"

// method
char method::projector[] = "pop";           // default = pop
char method::time_convolution[] = "tcl";    // default = tcl

// Hamiltonian 
double  Hamiltonian::w       = 0.0;
double  Hamiltonian::w01     = 0.0;
double  Hamiltonian::epsilon = 0.0;
double  Hamiltonian::Lambda  = 0.0;
int     Hamiltonian::N       = 1;
double  Hamiltonian::Omega   = 1.0;
double  Hamiltonian::y0      = 0.0;
double  Hamiltonian::wc      = 1.0;
double  Hamiltonian::eta     = 1.0;
double  Hamiltonian::beta    = 1.0e16;      // default = infinity
int     Hamiltonian::N_point = 100;
void Hamiltonian::print() const{
    cout << "w          " << w          << endl;
    cout << "w01        " << w01        << endl;
    cout << "epsilon    " << epsilon    << endl;
    cout << "Lambda     " << Lambda     << endl;
    cout << "N          " << N          << endl;
    cout << "Omega      " << Omega      << endl;
    cout << "y0         " << y0         << endl;
    cout << "wc         " << wc         << endl;
    cout << "eta        " << eta        << endl;
    cout << "beta       " << beta       << endl;
}
int Hamiltonian::Ohmic_mode_generator(){
    int Ns = N-1;                                           // number of the secondary mode
    arma::vec w = arma::vec(Ns,arma::fill::zeros);          // frequency of the 2nd modes
    arma::vec c = arma::vec(Ns,arma::fill::zeros);          // couple constant of the 2nd modes
    arma::mat Hessian = arma::mat(N,N,arma::fill::zeros);   // Hessian matrix for 1st+2nd modes
    w_  = arma::vec(N,arma::fill::zeros);                   // frequency of the normalized modes
    Req = arma::vec(N,arma::fill::zeros);                   // Equilibrium position of the normalized modes
    // generate the Ohmic spectral density
    for(int i=0;i<Ns;i++){
        w(i) = wc*log(Ns/(Ns-0.5-i));
        c(i) = sqrt(2*eta*wc/(pi*Ns))*w(i);
    }
    w.print("w:");
    c.print("c:");
    // generate the Hessian matrix
    arma::vec A = arma::pow(c,2)/arma::pow(w,2);
    Hessian(0,0) = pow(Omega,2) + arma::sum(A);
    for(int i=1;i<N;i++){
        Hessian(0,i) = c(i-1);
        Hessian(i,0) = c(i-1);
        Hessian(i,i) = pow(w(i-1),2);
    }
    //Hessian.print("Hessian:");
    // get the eigenvalue and eigenvector
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval,eigvec,Hessian);
    //eigval.print("eigval:");
    //eigvec.print("eigvec:");
    // check eigen
    arma::mat D = eigvec * arma::diagmat(eigval) * eigvec.t();
    arma::mat diff = Hessian - D;
    //D.print("D:");
    //diff.print("diff:");    // should be zeros
    // get w_ and Req
    arma::vec one(1,arma::fill::ones);
    arma::vec req = join_cols(one,-c/arma::pow(w,2));
    //req.print("req:");

    w_ = arma::pow(eigval,0.5);
    Req = 2*y0*eigvec.t()*req;
    w_.print("w_:");
    Req.print("Req:");

    return 0;
}
std::complex<double> Hamiltonian::V01_V10(double t, double tau) const{
    arma::vec vec2 = 0.5/h_ * w_%arma::pow(Req,2) % arma::pow(arma::tanh(0.5*beta*h_*w_),-1) % (arma::cos(w_*(t-tau))-1); // where is beta come from ???
    arma::vec vec3 = 0.5/h_ * w_%arma::pow(Req,2) % arma::sin(w_*(t-tau));
    std::complex<double> part1 = (exp( 1i*w*(t-tau))*pow(Lambda,2) + pow(epsilon,2) + epsilon*exp(1i*w*t)*Lambda + epsilon*exp(-1i*w*tau)*Lambda) * exp( 1i*w01*(t-tau));
    std::complex<double> part2 = exp(    arma::sum(vec2));
    std::complex<double> part3 = exp(-1i*arma::sum(vec3));
    return part1*part2*part3;
}
std::complex<double> Hamiltonian::V10_V01(double t, double tau) const{
    arma::vec vec2 = 0.5/h_ * w_%arma::pow(Req,2) % arma::pow(arma::tanh(0.5*beta*h_*w_),-1) % (arma::cos(w_*(t-tau))-1);
    arma::vec vec3 = 1.0/h_ * w_%arma::pow(Req,2) % (arma::sin(w_*t)-arma::sin(w_*tau)-0.5*arma::sin(w_*(t-tau)));
    std::complex<double> part1 = (exp(-1i*w*(t-tau))*pow(Lambda,2) + pow(epsilon,2) + epsilon*exp(1i*w*t)*Lambda + epsilon*exp(-1i*w*tau)*Lambda) * exp(-1i*w01*(t-tau));
    std::complex<double> part2 = exp(    arma::sum(vec2));
    std::complex<double> part3 = exp( 1i*arma::sum(vec3));
    return part1*part2*part3;
}

// timer
std::vector<double> timer::t_history;
timer::timer(double T_set, double dt_set, double t_set){
    T = T_set;
    dt = dt_set;
    t = t_set;
}
void timer::print() const{
    cout << "t/T =  " << t << "/" << T << endl;
}

// state
std::vector<arma::cx_mat> state::sigma_history;
void state::print() const{
    sigma.print("sigma:");
}
int state::check_normal(){
    double tr = real(arma::trace(sigma));
    double error_limit = 1e-10;
    double error_limit_sq = pow(error_limit,2);
    if((tr-1)*(tr-1) < error_limit_sq){
        return 1;
    }else{return 0;}
}
class state state::evolution(double dt, arma::cx_mat dsigma_dt) const{
    class state state_new;
    state_new.sigma = sigma + dsigma_dt * dt;
    return state_new;
}
class state & state::operator= (const class state & state_old){
    if(this != &state_old){
        this->sigma = state_old.sigma;
    }
    return *this;
}
double state::sigma_x(){
    std::complex<double> s_x = sigma(0,1)+sigma(1,0);
    return s_x.real();
}
double state::sigma_y(){
    std::complex<double> s_y = -1i*sigma(0,1)+1i*sigma(1,0);
    return s_y.real();
}
double state::sigma_z(){
    std::complex<double> s_z = sigma(0,0)-sigma(1,1);
    return s_z.real();
}
double state::sigma_tr(){
    std::complex<double> s_tr= arma::trace(sigma);
    return s_tr.real();
}
