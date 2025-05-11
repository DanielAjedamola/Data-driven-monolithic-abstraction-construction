/*
 * vehicle.cc
 *
 *  created on: 2023
 *      author: Daniel A.
 */

/*
 * information about this example is given in the readme file
 *
 */

#include <array>
#include <iostream>
#include <stdio.h>      /* printf */
#include <stdlib.h>     /* system, NULL, EXIT_FAILURE */
#include <unistd.h>     /* fork */
#include <fstream>     /* write to file*/
#include <sstream>
#include <string>
#include <cmath>
#include <list>

#include <random>
#include <chrono>
// #include <new>
// #include <cstdlib>

#include "gurobi_c++.h"

//#include "MatlabDataArray.hpp"
//#include "MatlabEngine.hpp"

/* memory profiling */
#include <sys/time.h>
#include <sys/resource.h>
struct rusage usage;

#include "cuddObj.hh"
#include "SymbolicSet.hh"
#include "cuddSymbolicset.hh"
#include "SymbolicModelGrowthBound.hh"

#include "TicToc.hh"
#include "RungeKutta4.hh"
#include "RungeKutta4 _distb.hh"
#include "FixedPoint.hh"

#ifndef M_PI
#define M_PI 3.14159265359
#endif

/* state space dim */
#define sDIM 3
#define iDIM 2

using namespace std;

/* data types for the ode solver */
typedef std::array<double,3> state_type;
typedef std::array<double,2> input_type;

/*Adaptive Method*/
typedef std::array<double,18> mat;
int global_flag = 0;
list<mat> map;

// double et = 0.05;
// double eta_x[sDIM]={et, et, et};
double eta_x[sDIM] = {.3, .3, .3};//{.0697, .0697, .2};
// double eta_x[sDIM]={10/6,10/6,2*(M_PI+0.4)/6};
const double samples = 1331; //pow(10,4);//1331;//2197;
// const double samples = 50000;
/* state space interval*/
double h = cbrt(samples);
double h_eta2[sDIM] = {eta_x[0]/h, eta_x[1]/h, eta_x[2]/h};
double sub_grid_size = *std::max_element(h_eta2, h_eta2 + sDIM);
double Lx = 1.1688 * 4;
double Lw = 1.21 * 4;
double w_max[sDIM] = {0.15, 0.15, 0.015}; // disturbance bound
double bias[sDIM] = {Lx*h_eta2[0]+Lw*w_max[0], Lx*h_eta2[1]+Lw*w_max[1], Lx*h_eta2[2]+Lw*w_max[2]}; //0.053 4*(Lx*h_eta2 + Lw*w_max)

/* sampling time */
const double tau = 0.25;
/* number of intermediate steps in the ode solver */
const int nint=5;
OdeSolver ode_solver(sDIM,nint,tau);
OdeSolver1 ode_solver1(sDIM,nint,tau);

/* we integrate the vehicle ode by tau sec (the result is stored in x)  */
auto  vehicle_post = [](state_type &x, input_type &u) -> void {

  /* the ode describing the vehicle */
  auto rhs =[](state_type& xx,  const state_type &x, input_type &u) {
      double alpha=std::atan(std::tan(u[1])/2.0);
      xx[0] = u[0]*std::cos(alpha+x[2])/std::cos(alpha);
      xx[1] = u[0]*std::sin(alpha+x[2])/std::cos(alpha);
      xx[2] = u[0]*std::tan(u[1]);
  };
  ode_solver(rhs,x,u);
};

/* sample with disturbance  */
auto  sample_post = [](state_type &x, input_type &u, input_type &w) -> void {
    
    /* the ode describing the vehicle */
    auto rhs =[](state_type& xx,   state_type &x, input_type &u, input_type &w) {
        double alpha=std::atan(std::tan(u[1])/2.0);
        xx[0] = u[0]*std::cos(alpha+x[2])/std::cos(alpha) + w[0];
        xx[1] = u[0]*std::sin(alpha+x[2])/std::cos(alpha) + w[1];
        xx[2] = u[0]*std::tan(u[1]) + w[2];
    };
    ode_solver1(rhs,x,u,w);
};

std::default_random_engine rng;
std::uniform_real_distribution<double> uniform_dist(-1.0, 1.0);

/* computation of the growth bound (the result is stored in r)  */
auto radius_post = [](state_type &r, input_type &u) {
    // TicToc t;
    
    // Define the range for the random values
    double min = 0.0;
    // if(global_flag == 0){
        try {
            // t.tic();
            // Create an environment
            GRBEnv env = GRBEnv(true);
            //remove all logs
            env.set(GRB_IntParam_OutputFlag, 0);
            env.set(GRB_IntParam_LogToConsole,0);
            //start environment
            env.start();
            
            // Create an empty model
            GRBModel model = GRBModel(env);
            // Create variables
            
            /* theta_1            
                M =             
                [m11 m12 m13;       
                m21 m22 m23;        
                m31 m32 m33];       
            */
            
            // diagonal values (can be anything) theta_1
            GRBVar m11 = model.addVar(-100, 100, 0, GRB_CONTINUOUS, "m11");
            GRBVar m22 = model.addVar(-100, 100, 0, GRB_CONTINUOUS, "m22");
            GRBVar m33 = model.addVar(-100, 100, 0, GRB_CONTINUOUS, "m33");
            
            // off-diagonal values (only positive) theta_1
            GRBVar m12 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "m12");
            GRBVar m13 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "m13");
            GRBVar m21 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "m21");
            GRBVar m23 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "m23");
            GRBVar m31 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "m31");
            GRBVar m32 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "m32");
            
            // boundary values for diagonals (for optimisation)
            GRBVar u11 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "u11");
            GRBVar u22 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "u22");
            GRBVar u33 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "u33");

            // // theta_2
            // GRBVar o1 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "o1");
            // GRBVar o2 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "o2");
            // GRBVar o3 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "o3");

            
            // Set objective: minimize c^T * theta(concatenated theta_1 and theta_2 by columns)
            // model.setObjective(u11 + m21 + m31 + m12 + u22 + m32 + m13 + m23 + u33 + o1 + o2 + o3, GRB_MINIMIZE);
            model.setObjective(u11 + m21 + m31 + m12 + u22 + m32 + m13 + m23 + u33, GRB_MINIMIZE);
            
            // Transitions modelled as from pre to post
            // rather than from first to next
            
            //setup other values
            state_type init_pre;
            state_type init_post;
            state_type other_pre;
            state_type other_post;
            state_type lhs;
            state_type rhs;
            // double random;
            input_type randI;
            input_type randO;
            // Add constraints by for loop
            for (int i = 0; i < h; i++)
            {
                // for (int j = 0; j < sDIM; j++)
                // {
                //     random = ((double) rand()) / (double) RAND_MAX;
                //     random = (2 * random - 1) * h_eta2[j];
                //     //random value +/- h of r for each dimension
                //     other_pre[j] = r[j] + random;
                // }
                // for (int j = 0; j < sDIM; j++)
                // {
                //     random = ((double) rand()) / (double) RAND_MAX;
                //     random = (2 * random - 1) * h_eta2[j];
                //     //random value +/- h of r for each dimension
                //     init_pre[j] = r[j] + random;
                // }
                // Access precomputed random numbers
                // int baseIndex = i * sDIM * 4;
                
                // for (int j = 0; j < sDIM; ++j) {
                //     other_pre[j] = r[j] + randNumbers[baseIndex++];
                //     init_pre[j] = r[j] + randNumbers[baseIndex++];
                // }
                for (int j = 0; j < sDIM; j++) {
                    double random = uniform_dist(rng);
                    random = (2 * random - 1) * h_eta2[j];
                    other_pre[j] = r[j] + random;
                }
                for (int j = 0; j < sDIM; j++) {
                    double random = uniform_dist(rng);
                    random = (2 * random - 1) * h_eta2[j];
                    init_pre[j] = r[j] + random;
                }
                init_post = init_pre;
                other_post = other_pre;
                for (int j = 0; j < sDIM; j++) 
                {
                randI[j] = min + (double)rand() / RAND_MAX * w_max[j];
                randO[j] = min + (double)rand() / RAND_MAX * w_max[j];
                }
                        
                sample_post(other_post,u,randI); //updates value other_post with actual post value
                sample_post(init_post, u,randO); //updates value init_post with actual post value
                
                for (int j = 0; j < sDIM; j++)
                {
                    rhs[j] = abs(other_post[j] - init_post[j]);
                    // lhs[j] = abs(other_pre[j] - init_pre[j]);
                    lhs[j] = abs(other_pre[j] - init_pre[j]) + tau*w_max[j];
                }
                
                // model.addConstr(-(lhs[0]*m11 + lhs[1]*m12 + lhs[2]*m13 + o1) + bias[0] <= -rhs[0]);
                // model.addConstr(-(lhs[0]*m21 + lhs[1]*m22 + lhs[2]*m23 + o2) + bias[1] <= -rhs[1]);
                // model.addConstr(-(lhs[0]*m31 + lhs[1]*m32 + lhs[2]*m33 + o3) + bias[2] <= -rhs[2]);
            
                model.addConstr(-(lhs[0]*m11 + lhs[1]*m12 + lhs[2]*m13) + bias[0] <= -rhs[0]);
                model.addConstr(-(lhs[0]*m21 + lhs[1]*m22 + lhs[2]*m23) + bias[1] <= -rhs[1]);
                model.addConstr(-(lhs[0]*m31 + lhs[1]*m32 + lhs[2]*m33) + bias[2] <= -rhs[2]);
            
            }
            
            //absolute values for optimising diagonals of matrix
            model.addConstr( m11 - u11 <= 0);
            model.addConstr(-m11 - u11 <= 0);
            model.addConstr( m22 - u22 <= 0);
            model.addConstr(-m22 - u22 <= 0);
            model.addConstr( m33 - u33 <= 0);
            model.addConstr(-m33 - u33 <= 0);
            
            // Optimize model
            model.optimize();
            
            // map.push_front({1,r[0],r[1],r[2],u[0],u[1],m11.get(GRB_DoubleAttr_X), m12.get(GRB_DoubleAttr_X),m13.get(GRB_DoubleAttr_X),m21.get(GRB_DoubleAttr_X), m22.get(GRB_DoubleAttr_X),m23.get(GRB_DoubleAttr_X),m31.get(GRB_DoubleAttr_X),m32.get(GRB_DoubleAttr_X),m33.get(GRB_DoubleAttr_X),o1.get(GRB_DoubleAttr_X),o2.get(GRB_DoubleAttr_X),o3.get(GRB_DoubleAttr_X)});
            
            
            
    
            // return r values based on optimised growth bound
            // r[0] = m11.get(GRB_DoubleAttr_X)*h_eta2[0]/2   + m12.get(GRB_DoubleAttr_X)*h_eta2[0]/2   + m13.get(GRB_DoubleAttr_X)*h_eta2[0]/2 + o1.get(GRB_DoubleAttr_X);
            // r[1] = m21.get(GRB_DoubleAttr_X)*h_eta2[1]/2   + m22.get(GRB_DoubleAttr_X)*h_eta2[1]/2   + m23.get(GRB_DoubleAttr_X)*h_eta2[1]/2 + o2.get(GRB_DoubleAttr_X);
            // r[2] = m31.get(GRB_DoubleAttr_X)*h_eta2[2]/2   + m32.get(GRB_DoubleAttr_X)*h_eta2[2]/2  + m33.get(GRB_DoubleAttr_X)*h_eta2[2]/2 + o3.get(GRB_DoubleAttr_X);            
            double h1 = 2;
            r[0] = m11.get(GRB_DoubleAttr_X)*h_eta2[0]/h1   + m12.get(GRB_DoubleAttr_X)*h_eta2[0]/h1   + m13.get(GRB_DoubleAttr_X)*h_eta2[0]/h1;
            r[1] = m21.get(GRB_DoubleAttr_X)*h_eta2[1]/h1   + m22.get(GRB_DoubleAttr_X)*h_eta2[1]/h1   + m23.get(GRB_DoubleAttr_X)*h_eta2[1]/h1;
            r[2] = m31.get(GRB_DoubleAttr_X)*h_eta2[2]/h1   + m32.get(GRB_DoubleAttr_X)*h_eta2[2]/h1  + m33.get(GRB_DoubleAttr_X)*h_eta2[2]/h1;            
            
            // double c = std::abs(u[0]*std::sqrt(std::tan(u[1])*std::tan(u[1])/4.0+1));
            // r[0] = r[0] + c*r[2]*0.3 + w_max[0]*tau + c*w_max[2]*tau*tau/2.0;
            // r[1] = r[1] + c*r[2]*0.3 + w_max[1]*tau + c*w_max[2]*tau*tau/2.0;
    
    // t.toc();
        }
        catch(GRBException e) {
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch(...) {
            cout << "Exception during optimization" << endl;
        }
    // }else{
    //     for (auto const& i : map) {
    //         if ((r[0] >= i[1] - h_eta2[0]/2) && (r[0] <= i[1] + h_eta2[0]/2) && (r[1] >= i[2] - h_eta2[1]/2) && (r[1] <= i[2] + h_eta2[1]/2) && (r[2] >= i[3] - h_eta2[2]/2) && (r[2] <= i[3] + h_eta2[2]/2) && (i[4] == u[0]) && (i[5] == u[1])){
    //             // r[0] = i[6]*h_eta2[0]/2    + i[7]*h_eta2[0]/2  +i[8]*h_eta2[0]/2 + i[15];
    //             // r[1] = i[9]*h_eta2[1]/2    + i[10]*h_eta2[0]/2 +i[11]*h_eta2[0]/2 +i[16];
    //             // r[2] = i[12]*h_eta2[2]/2   + i[13]*h_eta2[0]/2 +i[14]*h_eta2[0]/2 +i[17];
    //             double c = std::abs(u[0]*std::sqrt(std::tan(u[1])*std::tan(u[1])/4.0+1));
    //             r[0] = r[0] + c*r[2]*0.3 + randI[0]*tau + c*randI[2]*tau*tau/2.0;
    //             r[1] = r[1] + c*r[2]*0.3 + randI[1]*tau + c*randI[2]*tau*tau/2.0;
    //         }
    //     }
    // }
    // // compare with this 
    // for(int i=0; i<3; i++){
    //   cout << i << r[i] << "\n";
    // } 

};

/* forward declaration of the functions to setup the state space
 * input space and obstacles of the vehicle example */
scots::SymbolicSet vehicleCreateStateSpace(Cudd &mgr);
scots::SymbolicSet vehicleCreateInputSpace(Cudd &mgr);

void vehicleCreateObstacles(scots::SymbolicSet &obs);


// int run(Cudd &mgr){
int main() {
    /* to measure time */
    TicToc tt;
    Cudd mgr;
    /****************************************************************************/
    /* construct SymbolicSet for the state space */
    /****************************************************************************/
    scots::SymbolicSet ss=vehicleCreateStateSpace(mgr);
    ss.writeToFile("vehicle_ss.bdd");
    
    /****************************************************************************/
    /* construct SymbolicSet for the obstacles */
    /****************************************************************************/
    /* first make a copy of the state space so that we obtain the grid
     * information in the new symbolic set */
    scots::SymbolicSet obs(ss);
    vehicleCreateObstacles(obs);
    obs.writeToFile("vehicle_obst.bdd");
    
    /****************************************************************************/
    /* we define the target set */
    /****************************************************************************/
    /* first make a copy of the state space so that we obtain the grid
     * information in the new symbolic set */
    scots::SymbolicSet ts(ss);
    /* define the target set as a symbolic set */
    double H[4*sDIM]={-1, 0, 0,
                        1, 0, 0,
                        0,-1, 0,
                        0, 1, 0};
    /* compute inner approximation of P={ x | H x<= h1 }  */
    double h[4] = {-9,9.51,-0, .51};
    ts.addPolytope(4,H,h, scots::INNER);
    ts.writeToFile("vehicle_target.bdd");

    /****************************************************************************/
    /* construct SymbolicSet for the input space */
    /****************************************************************************/
    scots::SymbolicSet is=vehicleCreateInputSpace(mgr);
    
    /****************************************************************************/
    /* setup class for symbolic model computation */
    /****************************************************************************/
    /* first create SymbolicSet of post variables
     * by copying the SymbolicSet of the state space and assigning new BDD IDs */
    scots::SymbolicSet sspost(ss,1);
    /* instantiate the SymbolicModel */
    scots::SymbolicModelGrowthBound<state_type,input_type> abstraction(&ss, &is, &sspost);
    /* compute the transition relation */
    tt.tic();
    abstraction.computeTransitionRelation(vehicle_post, radius_post);
    std::cout << std::endl;
    tt.toc();

    if (!getrusage(RUSAGE_SELF, &usage)) {
        double memoryUsageGB = static_cast<double>(usage.ru_maxrss) / (1024 * 1024 * 1024); // Convert KB to GB
        std::cout << "Memory usage for the entire abstraction: " << memoryUsageGB << " GB" << std::endl;
    }
    /* get the number of elements in the transition relation */
    std::cout << std::endl << "Number of elements in the transition relation: " << abstraction.getSize() << std::endl;
    std::cout << std::endl << "Number of data sampled from the subgrid:" << samples << std::endl;
    std::cout << std::endl << "Norm of subgrid size:" << sub_grid_size << std::endl;

    /****************************************************************************/
    /* we continue with the controller synthesis */
    /****************************************************************************/
    /* we setup a fixed point object to compute reachabilty controller */
    scots::FixedPoint fp(&abstraction);
    /* the fixed point algorithm operates on the BDD directly */
    BDD T = ts.getSymbolicSet();
    BDD O = obs.getSymbolicSet();
    tt.tic();
    /* compute controller */
    BDD C=fp.reachAvoid(T,O,1);
    tt.toc();
    /****************************************************************************/
    /* last we store the controller as a SymbolicSet
     * the underlying uniform grid is given by the Cartesian product of
     * the uniform gird of the space and uniform gird of the input space */
    /****************************************************************************/
    scots::SymbolicSet controller(ss,is);
    controller.setSymbolicSet(C);
    cout << "Controller size: " << controller.getSize() << endl;
    if(controller.getSize() > 0){
        controller.writeToFile("vehicle_controller.bdd");
        return 1;}
        else{
        return 0;
    }
}

// Function to create the non-uniform grid
scots::SymbolicSet vehicleCreateStateSpace(Cudd &mgr) {

    /* lower bounds of the hyper rectangle */
    double lb[sDIM]={0,0,-M_PI-0.4};  
    /* upper bounds of the hyper rectangle */
    double ub[sDIM]={10,10,M_PI+0.4};  

    scots::SymbolicSet ss(mgr,sDIM,lb,ub,eta_x);

    /* add the grid points to the SymbolicSet ss */
    ss.addGridPoints();

    return ss;
}

void vehicleCreateObstacles(scots::SymbolicSet &obs) {
    
    /* add the obstacles to the symbolic set */
    /* the obstacles are defined as polytopes */
    /* define H* x <= h */
    double H[4*sDIM]={-1, 0, 0,
        1, 0, 0,
        0,-1, 0,
        0, 1, 0};
    /* add outer approximation of P={ x | H x<= h1 } form state space */
    double h1[4] = {-1,1.2,-0, 9};
    obs.addPolytope(4,H,h1, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h2 } form state space */
    double h2[4] = {-2.2,2.4,-0,5};
    obs.addPolytope(4,H,h2, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h3 } form state space */
    double h3[4] = {-2.2,2.4,-6,10};
    obs.addPolytope(4,H,h3, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h4 } form state space */
    double h4[4] = {-3.4,3.6,-0,9};
    obs.addPolytope(4,H,h4, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h5 } form state space */
    double h5[4] = {-4.6 ,4.8,-1,10};
    obs.addPolytope(4,H,h5, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h6 } form state space */
    double h6[4] = {-5.8,6,-0,6};
    obs.addPolytope(4,H,h6, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h7 } form state space */
    double h7[4] = {-5.8,6,-7,10};
    obs.addPolytope(4,H,h7, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h8 } form state space */
    double h8[4] = {-7,7.2,-1,10};
    obs.addPolytope(4,H,h8, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h9 } form state space */
    double h9[4] = {-8.2,8.4,-0,8.5};
    obs.addPolytope(4,H,h9, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h10 } form state space */
    double h10[4] = {-8.4,9.3,-8.3,8.5};
    obs.addPolytope(4,H,h10, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h11 } form state space */
    double h11[4] = {-9.3,10,-7.1,7.3};
    obs.addPolytope(4,H,h11, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h12 } form state space */
    double h12[4] = {-8.4,9.3,-5.9,6.1};
    obs.addPolytope(4,H,h12, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h13 } form state space */
    double h13[4] = {-9.3,10 ,-4.7,4.9};
    obs.addPolytope(4,H,h13, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h14 } form state space */
    double h14[4] = {-8.4,9.3,-3.5,3.7};
    obs.addPolytope(4,H,h14, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h15 } form state space */
    double h15[4] = {-9.3,10 ,-2.3,2.5};
    obs.addPolytope(4,H,h15, scots::OUTER);
    
}

scots::SymbolicSet vehicleCreateInputSpace(Cudd &mgr) {
    
    /* lower bounds of the hyper rectangle */
    double lb[sDIM]={-1,-1};
    /* upper bounds of the hyper rectangle */
    double ub[sDIM]={1,1};
    /* grid node distance diameter */
    double eta[sDIM]={.3,.3};
    
    scots::SymbolicSet is(mgr,iDIM,lb,ub,eta);
    is.addGridPoints();
    
    return is;
}
