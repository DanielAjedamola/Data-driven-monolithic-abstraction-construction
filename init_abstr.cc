/*
* Initializer of Abstraction, which generates the optimized grid parameter
* 
* created on: Aug 2023
*     author: Daniel A.
*/

#include <array>
#include <cmath>
#include <vector>
#include <random>
#include <iostream>

#include "gurobi_c++.h"
#include "scots.hh"
#include "cuddObj.hh"

#include "SymbolicSet.hh"
#include "SymbolicSets.hh"

#include "TicToc.hh"
#include "RungeKutta4.hh"
#include "RungeKutta4 _distb.hh"

#ifndef M_PI
#define M_PI 3.14159265359
#endif

/* state space dim */
#define sDIM 3
#define iDIM 2

/* data types for the ode solver */
typedef std::array<double,3> state_type;
typedef std::array<double,2> input_type;

// set up state-space initial parameters
/* lower bounds of the hyper rectangle */
double lb[sDIM]={0,0,-M_PI-0.4};  
/* upper bounds of the hyper rectangle */
double ub[sDIM]={10,10,M_PI+0.4}; 
/* grid node distance diameter */
double N = 28;
double n1 = 221;
// coarse grid to start with
// double eta_x[sDIM]={10/N,10/N,2*(M_PI+0.4)/n1};
double eta_x[sDIM]={1.25,1.25,.0015};
double gamma = std::log(eta_x[0]*eta_x[1]*eta_x[2]); // e^gamma = volume of cells

/* state space interval*/
int const samples = 1331;//2197;//2500;
double h = cbrt(samples);
double h_eta2[sDIM] = {eta_x[0]/h, eta_x[1]/h, eta_x[2]/h};
double Lx = 1.2 * 4;
double Lw = 1.21 * 4;
double w_max[sDIM] = {0.15, 0.15, 0.015}; // disturbance bound
double bias[sDIM] = {Lx*h_eta2[0]+Lw*w_max[0], Lx*h_eta2[1]+Lw*w_max[1], Lx*h_eta2[2]+Lw*w_max[2]}; //0.053 4*(Lx*h_eta2 + Lw*w_max)

/* sampling time */
const double tau = 0.25;
/* number of intermediate steps in the ode solver */
const int nint=5;
OdeSolver1 ode_solver1(sDIM,nint,tau);

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


// Define a struct to hold the output matrices
struct OutputMatrices {
    std::vector<std::vector<double>> M1;
    std::vector<double> M2;
};

OutputMatrices radius_post(const std::vector<double>& r, const std::vector<double>& u) {
    // Define the range for the random values
    double min = 0.0;
    // Convert u to input_type
    input_type u_array = {u[0], u[1]};
    try {
        // Create an environment
        GRBEnv env1 = GRBEnv(true);
        // Remove all logs
        env1.set(GRB_IntParam_OutputFlag, 0);
        env1.set(GRB_IntParam_LogToConsole, 0);
        // Start the environment
        env1.start();

        // Create an empty model
        GRBModel model1 = GRBModel(env1);

        // Define variables for M1 (3x3 matrix)
        std::vector<std::vector<GRBVar>> M1(sDIM, std::vector<GRBVar>(sDIM));
        for (int i = 0; i < sDIM; i++) {
            for (int j = 0; j < sDIM; j++) {
                if(i==j){M1[i][j] = model1.addVar(-200, 200, 0, GRB_CONTINUOUS, "m" + std::to_string(i + 1) + std::to_string(j + 1));}
                else{M1[i][j] = model1.addVar(0, 200, 0, GRB_CONTINUOUS, "m" + std::to_string(i + 1) + std::to_string(j + 1));}
            }
        }

        // // Define variables for M2 (vector)
        // std::vector<GRBVar> M2(sDIM);
        // for (int i = 0; i < sDIM; i++) {
        //     M2[i] = model1.addVar(0, 100, 0, GRB_CONTINUOUS, "o" + std::to_string(i + 1));
        // }

        // Define variables for the diagonal of M1 (theta_1)
        std::vector<GRBVar> theta1(sDIM);
        for (int i = 0; i < sDIM; i++) {
            theta1[i] = model1.addVar(0, 200, 0, GRB_CONTINUOUS, "u" + std::to_string(i + 1) + std::to_string(i + 1));
        }

        // Set objective: minimize c^T * theta (concatenated theta_1 and theta_2 by columns)
        GRBLinExpr objExpr;
        for (int i = 0; i < sDIM; i++) {
            for (int j = 0; j < sDIM; j++) {
                if(i!=j){objExpr += M1[i][j];}
            }
            // objExpr += M2[i];
            objExpr += theta1[i];
        }
        model1.setObjective(objExpr, GRB_MINIMIZE);

        // Setup other values
        state_type init_pre;
        state_type init_post;
        state_type other_pre;
        state_type other_post;
        state_type lhs;
        state_type rhs;
        double random;
        input_type randI;
        input_type randO;

        // Add constraints by a for loop
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < sDIM; j++) {
                random = ((double)rand()) / (double)RAND_MAX;
                random = (2 * random - 1) * h_eta2[j];
                // Random value +/- h of r for each dimension
                other_pre[j] = r[j] + random;
            }
            for (int j = 0; j < sDIM; j++) {
                random = ((double)rand()) / (double)RAND_MAX;
                random = (2 * random - 1) * h_eta2[j];
                // Random value +/- h of r for each dimension
                init_pre[j] = r[j] + random;
            }
            init_post = init_pre;
            other_post = other_pre;
            for (int j = 0; j < sDIM; j++) {
                randI[j] = min + (double)rand() / RAND_MAX * w_max[j];
                randO[j] = min + (double)rand() / RAND_MAX * w_max[j];
            }
            
            sample_post(other_post, u_array, randI); // Updates value other_post with actual post value
            sample_post(init_post, u_array, randO); // Updates value init_post with actual post value

            for (int j = 0; j < sDIM; j++) {
                rhs[j] = abs(other_post[j] - init_post[j]);
                // lhs[j] = abs(other_pre[j] - init_pre[j]);
                lhs[j] = abs(other_pre[j] - init_pre[j]) + tau*w_max[j];
            }

            for (int j = 0; j < sDIM; j++) {
                // model1.addConstr(-(lhs[0] * M1[j][0] + lhs[1] * M1[j][1] + lhs[2] * M1[j][2] + M2[j]) + bias[j] <= -rhs[j]);
                model1.addConstr(-(lhs[0] * M1[j][0] + lhs[1] * M1[j][1] + lhs[2] * M1[j][2]) + bias[j] <= -rhs[j]);
            }
        }

        // Absolute values for optimizing diagonals of the matrix
        for (int i = 0; i < sDIM; i++) {
            model1.addConstr(M1[i][i] - theta1[i] <= 0);
            model1.addConstr(-M1[i][i] - theta1[i] <= 0);
        }

        // Optimize the model
        model1.optimize();

        // Retrieve the optimized values for M1 (theta_1) and M2
        OutputMatrices output;
        output.M1.resize(sDIM, std::vector<double>(sDIM));
        for (int i = 0; i < sDIM; i++) {
            for (int j = 0; j < sDIM; j++) {
                if (i==j) {output.M1[i][j] = 1+M1[i][j].get(GRB_DoubleAttr_X);} // A=id+theta1
                else {output.M1[i][j] = M1[i][j].get(GRB_DoubleAttr_X);}
            }
        }
        output.M2.resize(sDIM);
        for (int i = 0; i < sDIM; i++) {
            // output.M2[i] = 2*M2[i].get(GRB_DoubleAttr_X); // p=2*theta_2
            output.M2[i] = 2*(M1[0][i].get(GRB_DoubleAttr_X)*w_max[0]*tau
                                + M1[1][i].get(GRB_DoubleAttr_X)*w_max[1]*tau
                                + M1[2][i].get(GRB_DoubleAttr_X)*w_max[2]*tau);
        }

        return output;
    }
    catch (GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
        return {}; // Return empty matrices on error
    }
    catch (...) {
        std::cout << "Exception during optimization" << std::endl;
        return {}; // Return empty matrices on error
    }
};

// OutputMatrices radius_post(const std::vector<double>& r, const std::vector<double>& u) {
    
//     // Convert u to input_type
//     input_type u_array = {u[0], u[1]};
//     std::vector<std::vector<double>> M1(sDIM, std::vector<double>(sDIM));
//     double c = std::abs(u[0]*std::sqrt(std::tan(u[1])*std::tan(u[1])/4.0+1));
//     M1[0][0] = 2;
//     M1[0][1] = 0;
//     M1[0][2] = c*tau;
//     M1[1][0] = 0;
//     M1[1][1] = 2;
//     M1[1][2] = c*tau;
//     M1[2][0] = 0;
//     M1[2][1] = 0;
//     M1[2][2] = 2;

//     // Define variables for M2 (vector)
//     std::vector<double> M2(sDIM);
//     M2[0] = w_max[0]*tau + w_max[2]*c*0.5*tau*tau;
//     M2[1] = w_max[1]*tau + w_max[1]*c*0.5*tau*tau;
//     M2[2] = w_max[2]*tau;

//     // Retrieve the optimized values for M1 and M2
//     OutputMatrices output;
//     output.M1.resize(sDIM, std::vector<double>(sDIM));
//     for (int i = 0; i < sDIM; i++) {
//         for (int j = 0; j < sDIM; j++) {
//             output.M1[i][j] = M1[i][j];
//         }
//     }
//     output.M2.resize(sDIM);
//     for (int i = 0; i < sDIM; i++) {
//         output.M2[i] = 2*M2[i]; // p=2*theta_2
//     }

//     return output;
// };

/* forward declaration of the functions to setup the state space
 * input space of the vehicle example */
scots::SymbolicSet vehicleCreateStateSpace(Cudd &mgr);
scots::SymbolicSet vehicleCreateInputSpace(Cudd &mgr);

// converting the grid encoded as BDDs to a vectors of points 
std::vector<double> bdd_to_grid_point(Cudd& manager, scots::SymbolicSet& symbolicSet) {
    BDD bdd = symbolicSet.getSymbolicSet();

    size_t dim = symbolicSet.getDimension();
    const double* eta = symbolicSet.getEta();
    const double* firstGridPoint = symbolicSet.getFirstGridPoint();

    // Find the number of grid points in the BDD
    size_t no_gp = symbolicSet.getSize();

    // Initialize the vector of grid points to be returned
    std::vector<double> gp(no_gp * dim);

    // Set up iteration to iterate over BDD cubes
    DdManager* dd = manager.getManager();
    int* cube;
    CUDD_VALUE_TYPE value;
    DdGen* gen;
    size_t counter = 0;

    // Iterate over BDD cubes
    Cudd_ForeachCube(dd, bdd.getNode(), gen, cube, value) {
        size_t offset = 1;
        for (size_t i = 0; i < dim; i++) {
            size_t no_vars = symbolicSet.getNofBddVars()[i];
            size_t* indBddVars = symbolicSet.getIndBddVars()[i];
            for (size_t j = 0; j < no_vars; j++) {
                size_t id = indBddVars[j];
                if (cube[id] == 1) {
                    for (size_t k = 0; k < offset; k++) {
                        gp[(counter + k) * dim + i] = firstGridPoint[i] + (size_t{1} << (no_vars - 1 - j)) * eta[i];
                    }
                }
                /* take care of don't care */
                if (cube[id] == 2) {
                    for (size_t k = 0; k < offset; k++) {
                        for (size_t l = 0; l <= i; l++) {
                            gp[(counter + k + offset) * dim + l] = gp[(counter + k) * dim + l];
                        }
                    }
                    for (size_t k = 0; k < offset; k++) {
                        gp[(counter + k + offset) * dim + i] = firstGridPoint[i] + (size_t{1} << (no_vars - 1 - j)) * eta[i];
                    }
                    offset = (offset << 1);
                }
            }
        }
        counter += offset;
    }
    // The vector 'gp' should contain the grid points, check if it's empty
    if (gp.empty()) {
        std::cout << "Error: Grid points vector is empty!" << std::endl;
    }
    return gp;
}

int main() {
    /* to measure time */
    TicToc tt;
    /* there is one unique manager to organize the bdd variables */
    Cudd mgr;

    tt.tic();
    /****************************************************************************/
    /* construct SymbolicSet for the state space */
    /****************************************************************************/
    scots::SymbolicSet ss = vehicleCreateStateSpace(mgr);
    std::vector<double> ss_vector = bdd_to_grid_point(mgr, ss);
    
    /****************************************************************************/
    /* construct SymbolicSet for the input space */
    /****************************************************************************/
    scots::SymbolicSet is = vehicleCreateInputSpace(mgr);
    std::vector<double> is_vector = bdd_to_grid_point(mgr, is);
    
    // // Sample ss_vector and is_vector
    // std::vector<double> ss_vector1 = { 1, 2, 3, 4, 5, 6 ,2,3,5};
    // std::vector<double> is_vector1 = { 7, 8, 9, 10,3,5};

    // Vectors to store the 3x3 matrices M1 and 3x1 vectors M2
    std::vector<std::vector<std::vector<double>>> M1_vector;
    std::vector<std::vector<double>> M2_vector;

    
    // Nested loop to iterate over each pair (s, i) in ss_vector x is_vector
    for (size_t s = 0; s < ss_vector.size(); s += 3) {
        for (size_t i = 0; i < is_vector.size(); i += 2) {
            // Extract the current (r, u) pair from ss_vector and is_vector
            std::vector<double> r(ss_vector.begin() + s, ss_vector.begin() + s + 3);
            std::vector<double> u(is_vector.begin() + i, is_vector.begin() + i + 2);

            // Calculate M1 and M2 for the current (r, u) pair
            OutputMatrices result = radius_post(r, u);

            // Store the calculated M1 and M2 in the vectors
            M1_vector.push_back(result.M1);
            M2_vector.push_back(result.M2);
            
            // // to print the values of pair (r,u)
            // std::cout << "r: ";
            // for (const auto& elem : r) {
            //     std::cout << elem << " ";
            // }
            // std::cout << std::endl;

            // // Print the elements of u
            // std::cout << "u: ";
            // for (const auto& elem : u) {
            //     std::cout << elem << " ";
            // }
            // std::cout << std::endl;
        }
    }

    // // Print M1_vector
    // std::cout << "M1_vector:" << std::endl;
    // for (const auto& matrix : M1_vector) {
    //     for (const auto& row : matrix) {
    //         for (double value : row) {
    //             std::cout << value << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << std::endl;
    // }

    // // Print M2_vector
    // std::cout << "M2_vector:" << std::endl;
    // for (const auto& row : M2_vector) {
    //     for (double value : row) {
    //         std::cout << value << " ";
    //     }
    //     std::cout << std::endl;
    // }
    tt.toc();
    std::cout << "Now optimising" << std::endl;
    try {
        // Set up Gurobi environment
        GRBEnv env = GRBEnv(true);
        env.set(GRB_IntParam_OutputFlag, 0);
        env.set(GRB_IntParam_LogToConsole, 0);
        // Start the environment
        env.start();

        // optimising the grid parameter for the least abstraction size.
        size_t J = M1_vector.size(); // Set the size of M1 and M2 here
        GRBModel model = GRBModel(env);
        // Decision variables bound
        double x_bd = 20.5;//0.50; 
        double n = -1;
        double y_bd = std::exp(x_bd);
        double z_bd = std::exp(2.0*x_bd);
        double s_bd = std::exp(3.0*x_bd);

        GRBVar x1 = model.addVar(n*x_bd, x_bd, 0, GRB_CONTINUOUS, "x1");
        GRBVar x2 = model.addVar(n*x_bd, x_bd, 0, GRB_CONTINUOUS, "x2");
        GRBVar x3 = model.addVar(n*x_bd, x_bd, 0, GRB_CONTINUOUS, "x3");
        //auxilliary variable y=exp(x)
        GRBVar y1 = model.addVar(0, y_bd, 0, GRB_CONTINUOUS, "y1");
        GRBVar y2 = model.addVar(0, y_bd, 0, GRB_CONTINUOUS, "y2");
        GRBVar y3 = model.addVar(0, y_bd, 0, GRB_CONTINUOUS, "y3");

        std::vector<GRBVar> z(6);
        for (int i = 0; i < 6; i++) {
            z[i] = model.addVar(0, z_bd, 0, GRB_CONTINUOUS, "z" + std::to_string(i + 1));
        }

        std::vector<GRBVar> s(10);
        for (int i = 0; i < 10; i++) {
            s[i] = model.addVar(0, s_bd, 0, GRB_CONTINUOUS, "s" + std::to_string(i + 1));
        }
        // // Create an auxiliary variable for constants * x_i
        GRBVar x11 = model.addVar(-2.0*x_bd, 2.0*x_bd, 0, GRB_CONTINUOUS, "x11");
        GRBVar x12 = model.addVar(-2.0*x_bd, 2.0*x_bd, 0, GRB_CONTINUOUS, "x12");
        GRBVar x13 = model.addVar(-2.0*x_bd, 2.0*x_bd, 0, GRB_CONTINUOUS, "x13");
        GRBVar x22 = model.addVar(-2.0*x_bd, 2.0*x_bd, 0, GRB_CONTINUOUS, "x22");
        GRBVar x23 = model.addVar(-2.0*x_bd, 2.0*x_bd, 0, GRB_CONTINUOUS, "x23");
        GRBVar x33 = model.addVar(-2.0*x_bd, 2.0*x_bd, 0, GRB_CONTINUOUS, "x33");
        GRBVar x111 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x111");
        GRBVar x112 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x112");
        GRBVar x113 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x113");
        GRBVar x122 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x122");
        GRBVar x123 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x123");
        GRBVar x222 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x222");
        GRBVar x133 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x133");
        GRBVar x233 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x233");
        GRBVar x223 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x223");
        GRBVar x333 = model.addVar(-3.0*x_bd, 3.0*x_bd, 0, GRB_CONTINUOUS, "x333");

        // Update the model
        model.update();

        // Objective function
        GRBQuadExpr objective = 0.0;
        // Calculate the exponential term
        double exp_term = std::exp(-1 * gamma);
        for (size_t j = 0; j < J; j++) {
            double m11 = M1_vector[j][0][0];
            double m12 = M1_vector[j][0][1];
            double m13 = M1_vector[j][0][2];
            double m21 = M1_vector[j][1][0];
            double m22 = M1_vector[j][1][1];
            double m23 = M1_vector[j][1][2];
            double m31 = M1_vector[j][2][0];
            double m32 = M1_vector[j][2][1];
            double m33 = M1_vector[j][2][2];

            double p1 = M2_vector[j][0];
            double p2 = M2_vector[j][1];
            double p3 = M2_vector[j][2];

            // Add the quadratic terms to the objective expression
            objective += p1 * p2 * p3 +
                        p1 * p2 * (m31 * y1) + p1 * p2 * (m32 * y2) + p1 * p2 * (m33 * y3) +
                        p1 * (m21 * y1) * p3 + p1 * m21 * (m31 * z[0]) + p1 * m21 * (m32 * z[1]) + p1 * m21 * (m33 * z[2]) +
                        p1 * (m22 * y2) * p3 + p1 * (m22 * z[1]) * m31 + p1 * m22 * (m32 * z[3]) + p1 * m22 * (m33 * z[4]) +
                        p1 * (m23 * y3) * p3 + p1 * (m23 * z[2]) * m31 + p1 * m23 * (m32 * z[4]) + p1 * m23 * (m33 * z[5]) +
                        (m11 * y1) * p2 * p3 + m11 * p2 * (m31 * z[0]) + m11 * p2 * (m32 * z[1]) + m11 * p2 * (m33 * z[2]) +
                        m11 * (m21 * z[0]) * p3 + m11 * m21 * (m31 * s[0]) + m11 * m21 * (m32 * s[1]) + m11 * m21 * (m33 * s[2]) +
                        m11 * (m22 * z[1]) * p3 + m11 * (m22 * s[1]) * m31 + m11 * m22 * (m32 * s[3]) + m11 * m22 * (m33 * s[4]) +
                        m11 * (m23 * z[2]) * p3 + m11 * (m23 * s[2]) * m31 + m11 * (m23 * s[4]) * m32 + (m11 * s[5]) * m23 * m33 +
                        (m12 * y2) * p2 * p3 + (m12 * z[1]) * p2 * m31 + m12 * p2 * (m32 * z[3]) + m12 * p2 * (m33 * z[4]) +
                        (m12 * z[1]) * m21 * p3 + (m12 * s[1]) * m21 * m31 + m12 * (m21 * s[3]) * m32 + m12 * (m21 * s[4]) * m33 +
                        m12 * (m22 * z[3]) * p3 + m12 * (m22 * s[3]) * m31 + m12 * (m22 * s[6]) * m32 + m12 * m22 * (m33 * s[7]) +
                        m12 * (m23 * z[4]) * p3 + (m12 * s[4]) * m23 * m31 + m12 * m23 * (m32 * s[7]) + m12 * (m23 * s[8]) * m33 +
                        (m13 * y3) * p2 * p3 + (m13 * z[2]) * p2 * m31 + (m13 * z[4]) * p2 * m32 + m13 * p2 * (m33 * z[5]) +
                        (m13 * z[2]) * m21 * p3 + (m13 * s[2]) * m21 * m31 + (m13 * s[4]) * m21 * m32 + m13 * (m21 * s[5]) * m33 +
                        m13 * (m22 * z[4]) * p3 + m13 * (m22 * s[4]) * m31 + (m13 * s[7]) * m22 * m32 + (m13 * s[8]) * m22 * m33 +
                        m13 * (m23 * z[5]) * p3 + m13 * (m23 * s[5]) * m31 + m13 * (m23 * s[8]) * m32 + m13 * m23 * (m33 * s[9]);

        }
        objective *= exp_term;
        model.setObjective(objective, GRB_MINIMIZE);

        // Constraint: x1 + x2 + x3 = gamma
        GRBLinExpr constraint123 = x1 + x2 + x3;
        model.addConstr(constraint123 == gamma);
        // constraints due to the auxilliary var
        model.addConstr(x11 == 2.0 * x1); // Add the constraint for x1_double = 2 * x1
        model.addConstr(x22 == 2.0 * x2);
        model.addConstr(x33 == 2.0 * x3);
        model.addConstr(x111 == 3.0 * x1);
        model.addConstr(x12 == x1 + x2);
        model.addConstr(x13 == x1 + x3);
        model.addConstr(x23 == x2 + x3);
        model.addConstr(x112 == x1 + x1 + x2);
        model.addConstr(x113 == x1 + x1 + x3);
        model.addConstr(x122 == x1 + x2 + x2);
        model.addConstr(x123 == x1 + x2 + x3);
        model.addConstr(x222 == x2 + x2 + x2);
        model.addConstr(x133 == x1 + x3 + x3);
        model.addConstr(x223 == x2 + x2 + x3);
        model.addConstr(x233 == x2 + x3 + x3);
        model.addConstr(x333 == x3 + x3 + x3);
        model.addGenConstrExp(x1, y1);
        model.addGenConstrExp(x2, y2);
        model.addGenConstrExp(x3, y3);
        model.addGenConstrExp(x11, z[0]);
        model.addGenConstrExp(x12, z[1]);
        model.addGenConstrExp(x13, z[2]);
        model.addGenConstrExp(x22, z[3]);
        model.addGenConstrExp(x23, z[4]);
        model.addGenConstrExp(x33, z[5]);
        model.addGenConstrExp(x111, s[0]);
        model.addGenConstrExp(x112, s[1]);
        model.addGenConstrExp(x113, s[2]);
        model.addGenConstrExp(x122, s[3]);
        model.addGenConstrExp(x123, s[4]);
        model.addGenConstrExp(x133, s[5]);
        model.addGenConstrExp(x222, s[6]);
        model.addGenConstrExp(x223, s[7]);
        model.addGenConstrExp(x233, s[8]);
        model.addGenConstrExp(x333, s[9]);
        
        tt.tic();
        // Optimize the model
        model.optimize();

        tt.toc();
        // Check the optimization status and retrieve the optimized values
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            double opt_x1 = x1.get(GRB_DoubleAttr_X);
            double opt_x2 = x2.get(GRB_DoubleAttr_X);
            double opt_x3 = x3.get(GRB_DoubleAttr_X);
            double eta_x_1 = y1.get(GRB_DoubleAttr_X);
            double eta_x_2 = y2.get(GRB_DoubleAttr_X);
            double eta_x_3 = y3.get(GRB_DoubleAttr_X);

            std::cout << "Optimal values: x1 = " << opt_x1 << ", x2 = " << opt_x2 << ", x3 = " << opt_x3 << std::endl;
            std::cout << "Optimal values: eta_x_1 = " << eta_x_1 << ", eta_x_2 = " << eta_x_2 << ", eta_x_3 = " << eta_x_3 << std::endl;
            std::cout << "for eta_x" << eta_x[0] << "," << eta_x[1] << "," << eta_x[2] << std::endl;
        } else {
            std::cout << "Optimization failed." << std::endl;
        }
    } catch (GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch (...) {
        std::cout << "Exception during optimization" << std::endl;
    }

    return 1;
}


scots::SymbolicSet vehicleCreateStateSpace(Cudd &mgr) {

    /* setup the workspace of the synthesis problem and the uniform grid */   
    // using pre-defined parameters
    scots::SymbolicSet ss(mgr,sDIM,lb,ub,eta_x);

    /* add the grid points to the SymbolicSet ss */
    ss.addGridPoints();

 return ss;
}

scots::SymbolicSet vehicleCreateInputSpace(Cudd &mgr) {

  /* lower bounds of the hyper rectangle */
  double lb[iDIM]={-1,-1};  
  /* upper bounds of the hyper rectangle */
  double ub[iDIM]={1,1}; 
  /* grid node distance diameter */
  double eta[iDIM]={.3,.3};   

  scots::SymbolicSet is(mgr,iDIM,lb,ub,eta);
  is.addGridPoints();

  return is;
}
