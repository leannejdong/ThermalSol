#include <Eigen/QR>
#include <Eigen/Dense>
#include <cmath>
#include "eigenData.h"
#include "Timer.h"
#include "util.h"
#include <typeinfo>
#include <stdexcept>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::cerr;
using Pair = std::vector<std::pair<int, int>>;
using std::pow;
const int Tg = 3;
const int T_initial = 80;

/// properties of water
const int cp = 4200;
const int rho = 980; const double PI{3.14159265};

namespace ThermalSol{
    void ThermalSolve(){
        std::cerr << LoadFile("data/verona.csv").value_or("File could not be opened!") << "\n";
        ifstream input("data/verona.csv");
        vector<vector<string>> data;
        string line, word;
        while (getline(input, line)) {
            stringstream ss(line);
            vector<string> row;
            while (getline(ss, word, ',')) {
                row.push_back(word);
            }
            data.push_back(row);
        }

//        for(const auto &row : data){
//            for(const auto &col: row){
//                std::cout << col << "\n";
//            }
//        }
//        std::cout << "\n";

        //! define the number of Nodes and Pipes
//        try{
//            const size_t n = stoul(data[0][0]);
//        } catch (std::invalid_argument const& ex){
//            std::cout << " a standard exception was caught, with message : " << ex.what() << "\n";
//        }
//
//        try{
//            const size_t m = stoul(data[0][1]);
//            std::cout << m << "\n";
//        } catch (std::invalid_argument const& ex){
//            std::cout << " a standard exception was caught, with message : " << ex.what() << "\n";
//        }
        const size_t n = stoul(data[0][0]);
        const size_t m = stoul(data[0][1]);


        //! Source nodes
        const size_t source = stoul(data[0][5]);
        vector<string> col4, col5;
        std::cerr << "the number of rows is " << m << "\n";
        for (size_t i = 0; i < m; ++i) {
            col4.push_back(data[i][3]);
            col5.push_back(data[i][4]);
        }

        vector<double> col4_int, col5_int;
        for (auto c2 : col4) {
            col4_int.push_back(stoi(c2));
        }
        //! c++ index from 0, so subtract 1 from each element of our vector
        for (auto &e2 : col4_int) {
            e2 -= 1.00;
        }

        for (auto c3 : col5) {
            col5_int.push_back(stoi(c3));
        }

        for (auto &e5 : col5_int) {
            e5 -= 1.00;
        }

        //!  create a matrix of in and out nodes for each time step
        VectorXd in = makeEigenVectorFromVectors(col4_int);
        VectorXd out = makeEigenVectorFromVectors(col5_int);
        MatrixXd in_node = MatrixXd::Zero(168, m);
        MatrixXd out_node = MatrixXd::Zero(168, m);
        in_node.rowwise() =   in.transpose();
        out_node.rowwise() = out.transpose();
        MatrixXd data_2 = openData("data/M.csv");
        int step =  data_2.rows();

        //! read pipes' lengths, diameter and U value from the Pipes_information file

        MatrixXd data_3 = openData("data/Pipes_information0.csv");
 //       VectorXd length = data_3.col(0);
 //       VectorXd diameters = data_3.col(2);
//        VectorXd U = data_3.col(4);
//        VectorXd A = pow(diameters.array(), 2)*PI/4;
//
//        //! time and spatial steps
//        const int dx = 20;
//        VectorXd delta_x = length/dx;
//        const int delta_t = 60;
//
//        //! Solution
//        MatrixXd T_deltat = MatrixXd::Zero(159, dx+1);
//        T_deltat.setConstant(80);
//        VectorXd T_t = VectorXd::Zero( dx+1);
//        VectorXd T_input = VectorXd::Zero(n);
//        T_input.setConstant(T_initial);
//        MatrixXd T_output = MatrixXd::Zero(168, m);
//        MatrixXd T_final_outlet = VectorXd::Zero( n);
//        MatrixXd results = MatrixXd::Zero(168, n);
//        MatrixXi in_node_Xi = in_node.cast<int>();
//        MatrixXi out_node_Xi = out_node.cast<int>();
//        for(int i{0}; i < step; ++i){
//            for(int j{0}; j < m; ++j){
//                if (data_2(i, j) < 0){
//                    int S = out_node_Xi(i, j);
//                    int N = in_node_Xi(i, j);
//                    out_node_Xi(i, j) = N;
//                    in_node_Xi(i, j) = S;
//                 }
//                double C1 = 2*delta_t*U(j)/A(j)/rho/cp;
//                double C2 = 2*std::fabs(data_2(i, j))*delta_t/rho/A(j)/delta_x(j);
//                double C = 1/(1+C1+C2);
//                T_t(0) = T_input(in_node_Xi(i, j));
//                for(int z{0}; z < 60; ++z){
//                    for(int k{0}; k < dx; ++k){
//                        T_t( k+1) = C*(T_deltat(j, k+1)+ C1*Tg + C2*T_t(k));
//                    }
//                    T_deltat.row(j) = T_t;
//                }
//                T_output(i, j) = T_t(dx);
//            }
//            for(int nodes{0}; nodes < n; ++nodes ){
//                int T_tot, G_tot, Count = 0;
//                for(int pipes{0}; pipes < m; ++ pipes){
//                    if(out_node(i, pipes) == nodes){
//                        ++Count;
//                        //T_final_outlet(nodes) = T_output(i, pipes);
//                    }
//                }
//                if(Count > 1){
//                    for(int pipes{0}; pipes < m; ++ pipes){
//                        if(out_node(i, pipes) == nodes){
//                            T_tot = + T_output(i, pipes)*std::fabs(data_2(i, pipes));
//                            G_tot = + std::fabs(data_2(i, pipes));
//                        }
//                    }
//                    T_final_outlet(nodes) = T_tot/G_tot;
//                }else{
//                    for(int pipes{0}; pipes < m; ++ pipes){
//                        if(out_node(i, pipes) == nodes){
//                            T_final_outlet(nodes) = T_output(i, pipes);
//                        }
//                    }
//                }
//            }
//            T_final_outlet(148) = 80;
//            T_final_outlet(149) = 80;
//            T_final_outlet(150) = 80;
//            T_input = T_final_outlet;
//            results.row(i-1) = T_final_outlet;
//        }








    }



}