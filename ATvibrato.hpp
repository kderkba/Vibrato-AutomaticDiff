/**************************************************************
    Code for 1st order antithetic Vibrato
    @file 
    @author kderkba
    @version
***************************************************************/

#ifndef ATVIBRATO_HPP
#define ATVIBRATO_HPP

#include "black_scholes.hpp"
#include "monte_carlo.hpp"
#include "diffusions.hpp"
#include <vector>
#include <random>
#include <ctime>
#include <cmath>
#include <iostream>
#include <functional>
#include <fstream>
#include <string>


/**************************************************************
    Function computes antithetic vibrato for the Black-Scholes Model
    @param M number of Monte Carlo paths
    @param Mz number of replication of random var Z for last 
    time step
    @param n time steps for Euler
    @param BS SDE of Black-scholes eq2.3, same for tangent
    @param C Black-Scholes option
    @param V Payoff function
    @param greek handles the special cases of theta and rho
    @return MeanVar object where mean is greek estimator
***************************************************************/
template <typename Payoff>
MeanVar AT_vibrato(unsigned M,unsigned Mz, unsigned n,SDE<> BS,
SDE<> tangent,BS_Option C,Payoff V,string greek=""){

    MeanVar expec(M*Mz,0.0,0.0);
    double expe_var = 0.0;
    double r = C.rate();
    
    double h = C.mat()/(double)n; //monte carlo time step
    std::normal_distribution<double> Z; //noise
    mt19937_64 G(time(NULL)); //RNG
    
    //Monte Carlo
    for(int i=0;i<M;i++){
        
        euler<> X(BS,h,0);
        euler<> Y(tangent,h,1);
        int k = 0;
        //computation euler Xn-1,Yn-1 eq2.4 eq2.12
        while(k<n-1){
            double rea = Z(G);
            Y.run(rea);
            X.run(rea);
            double tmp1 = X.State();
            double tmp2 = Y.State();
            Y.vecState() = {tmp1,tmp2};
            k++;
        }

        vector<double> state_X = X.vecState();
        vector<double> state_Y = Y.vecState();
        
        //equation 2.6, 2.13
        double mu_ = X.State() + BS.b(state_X)*h;
        double sig_ = BS.sigma(state_X)*sqrt(h);
        
        double dmu_ = Y.State() + tangent.b(state_Y)*h;
        double dsig_ = tangent.sigma(state_Y)*sqrt(h);


        //last time step, inner expectation
        // double inner_expec = 0.0;
        MeanVar inner_expec(Mz,0.0,0.0);
        // double inner_expec2 = 0.0;
        for(int j =0;j<Mz;j++){
            double rea = Z(G);

            double Xp = mu_ + sig_*rea;
            double Xm = mu_ - sig_*rea;
            double Xu = mu_;


            double Vp =V(Xp,C.strike());
            double Vm =V(Xm,C.strike());
            double Vu = V(Xu,C.strike());

            double E_mu = 0.5*(Vp - Vm)*(rea/sig_);
            double E_sig = 0.5*(Vp -2*Vu + Vm)*((rea*rea - 1)/sig_);
            
            if (greek == "rho"){
                double res = dmu_*E_mu + dsig_*E_sig - C.mat()*exp(-C.rate()*C.mat())*Vp;
                inner_expec.SumX() += res;
                inner_expec.SumXX() += res*res;
            } else if(greek == "theta"){
                double res = dmu_*E_mu + dsig_*E_sig - C.rate()*exp(-C.rate()*C.mat())*Vp;
                inner_expec.SumX() += res;
                inner_expec.SumXX() += res*res;
                //theta does not work. Seems like a term is missing
            } else{
                double res = dmu_*E_mu + dsig_*E_sig;
                inner_expec.SumX() += res;
                inner_expec.SumXX() += res*res;
            }
        }
        expe_var += inner_expec.var();
        // inner_expec/= (double)Mz;
        double sens = inner_expec.SumX()*exp(-r*C.mat()); //discounting sensibility
        expec.SumX() += sens;
        expec.SumXX() += sens*sens;
    }
    double var = (1/(double)M)*expec.var() + (1/(double)M*Mz)*(1/(double)M)*expe_var;
    expec.vrc() = var;
    return expec;
}

#endif