/////////////////////////////
//    driven_QME 3.0       //
// edit by Zongwei Huang   //
//       2023.5.30         //
/////////////////////////////
// This is a population-only non-condon version for the driven system QME
// version 1.1 update: checked with yifan and pouya's result, known error fixed
// version 1.2 update: add correlation function check
// version 2.0 update: rename the projecto to 'driven_QME'
//                     add time-convolution version for population only case
// version 3.0 update: add full version for time-convolution and time-convolutionless
#include <iostream>
#include "constant.h"
#include "ioput.h"
#include "initialize.h"
#include "evolution.h"

int main(){
    class Hamiltonian Hami_sim;
    class timer time_sim;
    class state stat_sim;
    int error_flag = 0;
    input(Hami_sim,time_sim,stat_sim);
    error_flag = initialize(Hami_sim,time_sim,stat_sim);
    if(error_flag == 1){return 0;}
    check_correlation_function(Hami_sim);
    evolution(Hami_sim,time_sim,stat_sim);
    //analysis();
    //output();
    
    return 0;
}
