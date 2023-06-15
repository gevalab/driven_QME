#include "evolution.h"

int evolution(class Hamiltonian & Hami_sim, class timer & time_sim, class state & stat_sim){
    std::cout << "evolution begin" << std::endl;
    ofstream outfile;
    outfile.open("sigma_z.txt",ios::app);
    arma::cx_mat k1,k2,k3,k4,k_sum;
    for(time_sim.t=0; time_sim.t<time_sim.T; time_sim.t+=time_sim.dt){
        // save the sigma history
        timer::t_history.push_back(time_sim.t);
        state::sigma_history.push_back(stat_sim.sigma);
        state::sigma_vec_history.push_back(stat_sim.mat_to_vec(stat_sim.sigma));
        // RK4
        k1 = dsigma_dt(Hami_sim,time_sim.t              ,stat_sim                            );
        k2 = dsigma_dt(Hami_sim,time_sim.t+time_sim.dt/2,stat_sim.evolution(time_sim.dt/2,k1));
        k3 = dsigma_dt(Hami_sim,time_sim.t+time_sim.dt/2,stat_sim.evolution(time_sim.dt/2,k2));
        k4 = dsigma_dt(Hami_sim,time_sim.t+time_sim.dt  ,stat_sim.evolution(time_sim.dt  ,k3));
        k_sum = (k1+2*k2+2*k3+k4)/6;
        stat_sim = stat_sim.evolution(time_sim.dt,k_sum);
        // print the sigma
        if(int(time_sim.t/time_sim.dt+0.5)%10==0){
            std::cout << time_sim.t << '\t' << stat_sim.sigma_z() << '\t' << stat_sim.sigma_tr() << std::endl;
            outfile   << time_sim.t << '\t' << stat_sim.sigma_z() << endl;
        }
    }
    return 0;
}

arma::cx_mat dsigma_dt(const class Hamiltonian & Hami_sim, double t, const class state & stat_sim){
    if(strcmp(method::projector,"pop")==0){
        if(strcmp(method::time_convolution,"tcl")==0){
            // pop-tcl
            std::complex<double> k01;   // rate constant
            std::complex<double> k10;
            // calculate k01 and k10
            k01 = k_calculation(Hami_sim,0,1,t,Hamiltonian::N_point);
            k10 = k_calculation(Hami_sim,1,0,t,Hamiltonian::N_point);
            // calculate dsigma_dt
            arma::cx_mat ds_dt = arma::cx_mat(2,2,arma::fill::zeros);
            ds_dt(0,0) = - k10*stat_sim.sigma(0,0) + k01*stat_sim.sigma(1,1);
            ds_dt(1,1) =   k10*stat_sim.sigma(0,0) - k01*stat_sim.sigma(1,1);
            return ds_dt;
        }else if(strcmp(method::time_convolution,"tc")==0){
            // pop-tc
            std::complex<double> I01;
            std::complex<double> I10;
            // calculate I01 and I10
            I01 = I_calculation(Hami_sim,0,1,t);
            I10 = I_calculation(Hami_sim,1,0,t);
            // calculate dsigma_dt
            arma::cx_mat ds_dt = arma::cx_mat(2,2,arma::fill::zeros);
            ds_dt(0,0) = - I10 + I01;
            ds_dt(1,1) =   I10 - I01;
            return ds_dt;
        }else{return 0;}
    }else if(strcmp(method::projector,"full")==0){
        if(strcmp(method::time_convolution,"tcl")==0){
            // full-tcl
            arma::cx_mat L = arma::cx_mat(4,4,arma::fill::zeros);
            arma::cx_mat k = arma::cx_mat(4,4,arma::fill::zeros);
            // calculate L and k
            L = Hami_sim.L(t);
            k = k_calculation_full(Hami_sim,t,Hamiltonian::N_point);
            // calculate dsigma_dt
            arma::cx_mat ds_dt_mat = arma::cx_mat(2,2,arma::fill::zeros);
            arma::cx_mat ds_dt_vec = arma::cx_mat(4,1,arma::fill::zeros);
            arma::cx_mat sigma_vec = stat_sim.mat_to_vec(stat_sim.sigma);
            ds_dt_vec = -1i/h_ * L * sigma_vec - k * sigma_vec;
            ds_dt_mat = stat_sim.vec_to_mat(ds_dt_vec);
            return ds_dt_mat;
        }else if(strcmp(method::time_convolution,"tc")==0){
            // full-tc
            arma::cx_mat L = arma::cx_mat(4,4,arma::fill::zeros);
            arma::cx_mat I = arma::cx_mat(4,4,arma::fill::zeros);
            // calculate L and I
            L = Hami_sim.L(t);
            I = I_calculation_full(Hami_sim,t);
            // calculate dsigma_dt
            arma::cx_mat ds_dt_mat = arma::cx_mat(2,2,arma::fill::zeros);
            arma::cx_mat ds_dt_vec = arma::cx_mat(4,1,arma::fill::zeros);
            arma::cx_mat sigma_vec = stat_sim.mat_to_vec(stat_sim.sigma);
            ds_dt_vec = -1i/h_ * L * sigma_vec - I;
            ds_dt_mat = stat_sim.vec_to_mat(ds_dt_vec);
            return ds_dt_mat;
        }else{return 0;}
    }else{return 0;}
}
arma::cx_mat k_calculation_full(const class Hamiltonian & Hami_sim, double t, int N){
    double h = (t-0)/N;
    arma::cx_mat sum = arma::cx_mat(4,4,arma::fill::zeros);
    // simpson integral
    sum = Hami_sim.K(t,0) + Hami_sim.K(t,t);
    for(int i=1;i<N;i++){
        if(i%2==1){
            sum += 4.*Hami_sim.K(t,0+i*h);
        }else{
            sum += 2.*Hami_sim.K(t,0+i*h);
        }
    }
    sum *= (h/3);
    sum *= pow(h_,-2);
    return sum;
}
arma::cx_mat I_calculation_full(const class Hamiltonian & Hami_sim, double t){
    int N = timer::t_history.size();
    double h = (timer::t_history[N-1]-timer::t_history[0])/N;
    double h_tail = t - timer::t_history[N-1];
    arma::cx_mat sum = arma::cx_mat(4,4,arma::fill::zeros);
    // trapzoid integral
    sum = 0.5*Hami_sim.K(t,timer::t_history[0])*state::sigma_vec_history[0] + 0.5*Hami_sim.K(t,timer::t_history[N-1])*state::sigma_vec_history[N-1];
    for(int i=1;i<N-1;i++){
        sum += Hami_sim.K(t,timer::t_history[i])*state::sigma_vec_history[i];
    }
    sum *= h;
    sum += h_tail * Hami_sim.K(t,timer::t_history[N-1])*state::sigma_vec_history[N-1];
    sum *= pow(h_,-2);
    return sum;
}
std::complex<double> I_calculation(const class Hamiltonian & Hami_sim, int j, int k, double t){
    int N = timer::t_history.size();
    double h = (timer::t_history[N-1]-timer::t_history[0])/N;
    double h_tail = t - timer::t_history[N-1];
    std::complex<double> sum = 0;
    std::complex<double> Ijk = 0;
    if(j==0 && k==1){
        // trapzoid integral
        sum = 0.5*Hami_sim.V10_V01(t,timer::t_history[0])*state::sigma_history[0](1,1) + 0.5*Hami_sim.V10_V01(t,timer::t_history[N-1])*state::sigma_history[N-1](1,1);
        for(int i=1;i<N-1;i++){
            sum += Hami_sim.V10_V01(t,timer::t_history[i])*state::sigma_history[i](1,1);
        }
        sum *= h;
        sum += h_tail * Hami_sim.V10_V01(t,timer::t_history[N-1])*state::sigma_history[N-1](1,1);
        // simpson integral
        /*
        if(N%2==1){
            sum = Hami_sim.V10_V01(t,timer::t_history[0])*state::sigma_history[0](1,1) + Hami_sim.V10_V01(t,timer::t_history[N-1])*state::sigma_history[N-1](1,1);
            for(int i=1;i<N-1;i++){
                if(i%2==1){
                    sum += 4.*Hami_sim.V10_V01(t,timer::t_history[i])*state::sigma_history[i](1,1);
                }else{
                    sum += 2.*Hami_sim.V10_V01(t,timer::t_history[i])*state::sigma_history[i](1,1);
                }
            }
            sum *= (h/3);
        }else{
            sum = Hami_sim.V10_V01(t,timer::t_history[0])*state::sigma_history[0](1,1) + Hami_sim.V10_V01(t,timer::t_history[N-2])*state::sigma_history[N-2](1,1);
            for(int i=1;i<N-2;i++){
                if(i%2==1){
                    sum += 4.*Hami_sim.V10_V01(t,timer::t_history[i])*state::sigma_history[i](1,1);
                }else{
                    sum += 2.*Hami_sim.V10_V01(t,timer::t_history[i])*state::sigma_history[i](1,1);
                }
            }
            sum *= (h/3);
            sum += h*0.5*(Hami_sim.V10_V01(t,timer::t_history[N-2])*state::sigma_history[N-2](1,1)+Hami_sim.V10_V01(t,timer::t_history[N-1])*state::sigma_history[N-1](1,1));
        }
        */
        Ijk = 2*pow(h_,-2)*sum.real();
        return Ijk;
    }else if(j==1 && k==0){
        sum = 0.5*Hami_sim.V01_V10(t,timer::t_history[0])*state::sigma_history[0](0,0) + 0.5*Hami_sim.V01_V10(t,timer::t_history[N-1])*state::sigma_history[N-1](0,0);
        for(int i=1;i<N-1;i++){
            sum += Hami_sim.V01_V10(t,timer::t_history[i])*state::sigma_history[i](0,0);
        }
        sum *= h;
        sum += h_tail * Hami_sim.V01_V10(t,timer::t_history[N-1])*state::sigma_history[N-1](0,0);
        Ijk = 2*pow(h_,-2)*sum.real();
        return Ijk;
    }else{return 0;}
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