#include "evolution.h"

int evolution(class Hamiltonian & Hami_sim, class timer & time_sim, class state & stat_sim){
    //RK4
    ofstream outfile;
    outfile.open("sigma_z.txt",ios::app);
    arma::cx_mat k1,k2,k3,k4,k_sum;
    for(time_sim.t=0; time_sim.t<time_sim.T; time_sim.t+=time_sim.dt){
        k1 = dsigma_dt(Hami_sim,time_sim.t              ,stat_sim                            );
        k2 = dsigma_dt(Hami_sim,time_sim.t+time_sim.dt/2,stat_sim.evolution(time_sim.dt/2,k1));
        k3 = dsigma_dt(Hami_sim,time_sim.t+time_sim.dt/2,stat_sim.evolution(time_sim.dt/2,k2));
        k4 = dsigma_dt(Hami_sim,time_sim.t+time_sim.dt  ,stat_sim.evolution(time_sim.dt  ,k3));
        k_sum = (k1+2*k2+2*k3+k4)/6;
        stat_sim = stat_sim.evolution(time_sim.dt,k_sum);
        if(int(time_sim.t/time_sim.dt+0.5)%10==0){
            std::cout << time_sim.t << '\t' << stat_sim.sigma_z() << '\t' << stat_sim.sigma_tr() << std::endl;
            outfile   << time_sim.t << '\t' << stat_sim.sigma_z() << endl;
        }
    }
    return 0;
}

arma::cx_mat dsigma_dt(const class Hamiltonian & Hami_sim, double t, const class state & stat_sim){
    std::complex<double> k01;   // rate constant
    std::complex<double> k10;
    // calculate k01 and k10
    k01 = k_calculation(Hami_sim,0,1,t,Hamiltonian::N_point);
    k10 = k_calculation(Hami_sim,1,0,t,Hamiltonian::N_point);
    // calculate dsigma_dt
    arma::cx_mat ds_dt = arma::cx_mat(2,2,arma::fill::zeros);
    ds_dt(0,0) = -k10*stat_sim.sigma(0,0) + k01*stat_sim.sigma(1,1);
    ds_dt(1,1) =  k10*stat_sim.sigma(0,0) - k01*stat_sim.sigma(1,1);
    return ds_dt;
}

std::complex<double> k_calculation(const class Hamiltonian & Hami_sim, int j, int k, double t, int N){
    double h = (t-0)/N;
    std::complex<double> sum = 0;
    std::complex<double> kjk = 0;
    // simpson integral
    if(j==0 && k==1){
        // k01 ~ <V10_V01>
        // trapzoid integral
        /*
        sum = 0.5*Hami_sim.V10_V01(t,0) + 0.5*Hami_sim.V10_V01(t,t);
        for(int i=1;i<N;i++){
            sum += Hami_sim.V10_V01(t,0+i*h);
        }
        sum *= h;
        */
        // simpson integral
        sum = Hami_sim.V10_V01(t,0) + Hami_sim.V10_V01(t,t);
        for(int i=1;i<N;i++){
            if(i%2==1){
                sum += 4.*Hami_sim.V10_V01(t,0+i*h);
            }else{
                sum += 2.*Hami_sim.V10_V01(t,0+i*h);
            }
        }
        sum *= (h/3);
        kjk = 2*pow(h_,-2)*sum.real();
        return kjk;
    }else if(j==1 && k==0){
        // k10 ~ <V01_V10>
        // trapzoid integral
        /*
        sum = 0.5*Hami_sim.V01_V10(t,0) + 0.5*Hami_sim.V01_V10(t,t);
        for(int i=1;i<N;i++){
            sum += Hami_sim.V01_V10(t,0+i*h);
        }
        sum *= h;
        */
        // simpson intnegral
        sum = Hami_sim.V01_V10(t,0) + Hami_sim.V01_V10(t,t);
        for(int i=1;i<N;i++){
            if(i%2==1){
                sum += 4.*Hami_sim.V01_V10(t,0+i*h);
            }else{
                sum += 2.*Hami_sim.V01_V10(t,0+i*h);
            }
        }
        sum *= (h/3);
        kjk = 2*pow(h_,-2)*sum.real();
        return kjk;
    }else{return 0;}
}