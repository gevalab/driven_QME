/////////////////////////////
//    driven_QME 4.0       //
// edit by Zongwei Huang   //
//       2023.7.25         //
/////////////////////////////
// This is a population-only non-condon version for the driven system QME
// version 1.1 update: checked with yifan and pouya's result, known error fixed
// version 1.2 update: add correlation function check
// version 2.0 update: rename the projecto to 'driven_QME'
//                     add time-convolution version for population only case
// version 3.0 update: add full version for time-convolution and time-convolutionless
// version 3.1 update: add new ohmic spectrum generator from (J. Phys. Chem. Lett. 2022, 13, 2330−2337)
//                     add steady state calculator (pop-tcl)
//                     add run time function
// version 4.0 update: add non-rotating frame version
// version 4.0.1 update: correct the V01_V01 and V10_V10 function in Hamiltonian

#include <iostream>
#include <chrono>
#include "constant.h"
#include "ioput.h"
#include "initialize.h"
#include "evolution.h"

int main(){
    class Hamiltonian Hami_sim;
    class timer time_sim;
    class state stat_sim;
    int error_flag = 0;
    auto start_time = std::chrono::high_resolution_clock::now();

    input(Hami_sim,time_sim,stat_sim);
    error_flag = initialize(Hami_sim,time_sim,stat_sim);
    if(error_flag == 1){return 0;}
    check_correlation_function(Hami_sim);
    if(strcmp(method::steady_state,"yes")==0){steady_state(Hami_sim,time_sim);}
    evolution(Hami_sim,time_sim,stat_sim);
    //analysis();
    //output();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::minutes>(end_time-start_time);
    std::cout << "total run time: " << duration.count() << " mins" << std::endl;
    
    return 0;
}
