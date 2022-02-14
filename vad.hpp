/**************************************************************
    Code for Vibrato + Automatic Differentiation
    Using Ddoule2 by Guillaume Sall
    @file 
    @author kderkba
    @version
***************************************************************/


#ifndef VAD_HPP
#define VAD_HPP

#include "black_scholes.hpp"
#include "monte_carlo.hpp"
#include "diffusions.hpp"
#include "./ad_library/Ddouble2.hpp"
#include <vector>
#include <random>
#include <ctime>
#include <cmath>
#include <iostream>



using namespace std;


/**************************************************************
    Function computes VAD for the Black-Scholes Model
    @param M number of Monte Carlo paths
    @param Mz number of replication of random var Z for last 
    time step
    @param n time steps for Euler
    @param BS SDE of Black-scholes eq2.3, same for tangent
    @param C Black-Scholes option
    @param V Payoff function
    @param expec MeanVar struct to compute mean,var, confidence
    for Monte Carlo
    @return nothing. It updates MeanVar objects by reference
***************************************************************/
template<typename Payoff>
void VAD(int M,int Mz,int n,SDE<> BS,SDE<> tangent, BS_Option C,Payoff V,MeanVar& expec,MeanVar& dexpec){
    
    
    ddouble K(C.strike()); //to be able to pass it into payoff function
    
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

        double inner_V = 0.0;// //sum of first derivatives
        double inner_dVdx = 0.0; //sum of second derivatives
        for(int j=0;j<Mz;j++){
            double rea = Z(G);
            ddouble asset(state_X[0],state_Y[1]); //storing X and Y = dX/dtheta
        
            // ddouble mu(state_X[0]*(1 + C.rate()*h),state_Y[1]*(1+C.rate()*h));
            // ddouble sig(sqrt(h)*C.sigma()*state_X[0],state_Y[1]*sqrt(h)*C.sigma());

            ddouble mu(state_X[0] + BS.b(state_X)*h,state_Y[1] + tangent.b(state_Y)*h);
            ddouble sig(sqrt(h)*BS.sigma(state_X),tangent.sigma(state_Y)*sqrt(h));

            ddouble Vtp = V(mu + sig*rea,K);
            ddouble Vtm = V(mu - sig*rea,K);
            ddouble Vtd = V(mu,K);

            double dmu_ = mu[1]; //derivative of mu at n-1 wrt theta
            double dsig_ = sig[1]; ////derivative of sigma at n-1 wrt theta

            // computation of dV/dtheta and storing its derivative, eq 4.3
            ddouble dV = exp(-C.rate()*C.mat())*(dmu_*0.5*(Vtp - Vtm)*(rea/(asset*C.sigma()*sqrt(h))) +
                            dsig_*0.5*(Vtp - 2*Vtd + Vtm)*((rea*rea - 1)/(asset*C.sigma()*sqrt(h))));
        
            inner_V += dV[0];//updating sum of 1st derivatives
            inner_dVdx += dV[1]; //updating sum of 2nd derivatives
        }
        
        inner_V/=(double)Mz; //mean
        inner_dVdx /=(double)Mz; //mean
        
        expec.SumX() += inner_V;
        expec.SumXX() += (inner_V*inner_V); //updating variance of first derivative

        dexpec.SumX() += inner_dVdx;
        dexpec.SumXX() += (inner_dVdx*inner_dVdx); //updating variance of second derivative
    }
    // cout << "1st derivative  "<<expec <<endl;
    // cout << "2nd derivative "<<dexpec <<endl;
    // return expec;
}
#endif