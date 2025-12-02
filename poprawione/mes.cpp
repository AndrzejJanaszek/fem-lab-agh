#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <string>
#include <iomanip>
#include <cmath>
#include <map>
#include <algorithm>


#include "model.h"
#include "gauss.h"
#include "solver.h"

// założenia programu:
// element 4 węzłowy
// liczba możliwych punktów całkowania: 2,3,4

void calculate_jacobian_and_dN_dX_dY(Grid &grid, int npc = 2){
    for(auto &element : grid.elements){
        element.jacobian.J.clear();
        element.jacobian.J.reserve(4);
        element.jacobian.J_1.clear();
        element.jacobian.J_1.reserve(4);

        for(int i = 0; i < npc*npc; i++){
            element.jacobian.J.push_back(std::array<double, 4>());
            element.jacobian.J_1.push_back(std::array<double, 4>());
        }

        for(int i_ksi = 0; i_ksi < npc; i_ksi++){
            for(int i_eta = 0; i_eta < npc; i_eta++){
                // dla każdego punktu całkowania [npc*i_ksi + i_eta]
                
                double dx_dKsi = 0;
                double dx_dEta = 0;
                // dla każdego weżła w elemencie
                for(int i = 0; i < 4; i++){
                    // x * dN_i/dKsi
                    dx_dKsi += grid.nodes[element.node_ids[i]-1].x * GaussQuad::uniEl.dN_dKsi[npc*i_ksi + i_eta][i];
                    // x * dN_i/dEta
                    dx_dEta += grid.nodes[element.node_ids[i]-1].x * GaussQuad::uniEl.dN_dEta[npc*i_ksi + i_eta][i];
                }

                double dy_dKsi = 0;
                double dy_dEta = 0;
                // dla każdego weżła w elemencie
                for(int i = 0; i < 4; i++){
                    // y * dN_i/dKsi
                    dy_dKsi += grid.nodes[element.node_ids[i]-1].y * GaussQuad::uniEl.dN_dKsi[npc*i_ksi + i_eta][i];
                    // y * dN_i/dEta
                    dy_dEta += grid.nodes[element.node_ids[i]-1].y * GaussQuad::uniEl.dN_dEta[npc*i_ksi + i_eta][i];
                }

                // jakobian
                element.jacobian.J[npc*i_ksi + i_eta][0] = dx_dKsi;
                element.jacobian.J[npc*i_ksi + i_eta][1] = dy_dKsi;
                element.jacobian.J[npc*i_ksi + i_eta][2] = dx_dEta;
                element.jacobian.J[npc*i_ksi + i_eta][3] = dy_dEta;

                // det J
                element.jacobian.detJ = dx_dKsi * dy_dEta - (dy_dKsi * dx_dEta);
                double detJ = dx_dKsi * dy_dEta - (dy_dKsi * dx_dEta);
                // odwórcoy
                element.jacobian.J_1[npc*i_ksi + i_eta][0] = dy_dEta / detJ;
                element.jacobian.J_1[npc*i_ksi + i_eta][1] = -dy_dKsi / detJ;
                element.jacobian.J_1[npc*i_ksi + i_eta][2] = -dx_dEta / detJ;
                element.jacobian.J_1[npc*i_ksi + i_eta][3] = dx_dKsi / detJ;
            }
        }
    }

    for(auto &element: grid.elements){
        // dla każdego punktu całkowania
        for(int i_ksi = 0; i_ksi < npc; i_ksi++){
            for(int i_eta = 0; i_eta < npc; i_eta++){
                std::array<double, 4> dNi_dx;
                std::array<double, 4> dNi_dy;

                for(int i = 0; i < 4; i++){
                    // dN/dx
                    double dN_dx = 
                    element.jacobian.J_1[npc*i_ksi + i_eta][0] * GaussQuad::uniEl.dN_dKsi[npc*i_ksi + i_eta][i] +
                    element.jacobian.J_1[npc*i_ksi + i_eta][1] * GaussQuad::uniEl.dN_dEta[npc*i_ksi + i_eta][i];
                    dNi_dx[i] = dN_dx;

                    // dN/dy
                    double dN_dy = 
                    element.jacobian.J_1[npc*i_ksi + i_eta][2] * GaussQuad::uniEl.dN_dKsi[npc*i_ksi + i_eta][i] +
                    element.jacobian.J_1[npc*i_ksi + i_eta][3] * GaussQuad::uniEl.dN_dEta[npc*i_ksi + i_eta][i];
                    dNi_dy[i] = dN_dy;
                }

                element.dN_dX[npc*i_ksi + i_eta] = dNi_dx;
                element.dN_dY[npc*i_ksi + i_eta] = dNi_dy;
                // printf("dla %d punktu całkowania:\n", npc*i_ksi + i_eta);
                // printf("dNi/dx:\n");
                // printf("%lf, %lf, %lf, %lf\n", dNi_dx[0], dNi_dx[1], dNi_dx[2], dNi_dx[3]);
                // printf("dNi/dy:\n");
                // printf("%lf, %lf, %lf, %lf\n", dNi_dy[0], dNi_dy[1], dNi_dy[2], dNi_dy[3]);
            }
        }
    }
}

template<size_t N>
void calculate_H_matrix(Grid &grid, GlobalData &global_data, const std::array<double, N> &points, int npc = 2)
{
    for (auto &element : grid.elements)
    {
        double H[4][4] = {0};
        double C[4][4] = {0};

        // dla kazdego punktu calkowania
        for (int i_pc = 0; i_pc < npc; i_pc++){
            for (int j_pc = 0; j_pc < npc; j_pc++){
                // todo liczyc raz gdzieś indziej
                double N_arr[4] = {0};
                double x_ksi = points[i_pc*2];
                double y_eta = points[j_pc*2];
                N_arr[0] = 0.25*(1-x_ksi)*(1-y_eta);
                N_arr[1] = 0.25*(1+x_ksi)*(1-y_eta);
                N_arr[2] = 0.25*(1+x_ksi)*(1+y_eta);
                N_arr[3] = 0.25*(1-x_ksi)*(1+y_eta);

                for (int row = 0; row < 4; row++){
                    for (int col = 0; col < 4; col++){
                        // macierz H
                        H[row][col] +=
                        (element.dN_dX[npc*i_pc + j_pc][row] * element.dN_dX[npc*i_pc + j_pc][col] +
                        element.dN_dY[npc*i_pc + j_pc][row] * element.dN_dY[npc*i_pc + j_pc][col]) *
                        global_data.conductivity * element.jacobian.detJ *
                        GaussQuad::points_2[npc*i_pc+1] * GaussQuad::points_2[npc*j_pc+1];
                        
                        // macierz C
                        C[row][col] +=
                        (N_arr[row] * N_arr[col] +
                        N_arr[row] * N_arr[col]) *
                        global_data.specific_heat * global_data.density *
                        element.jacobian.detJ *
                        GaussQuad::points_2[npc*i_pc+1] * GaussQuad::points_2[npc*j_pc+1];
                    }
                }
            }
        }

        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++){
                element.H[i][j] = H[i][j];
                element.C[i][j] = C[i][j];
            }
        }
    }
}

void calculate_HG_matrix(Grid &grid){    
    // inicjalizacja GH z wartościami 0
    std::vector<std::vector<double>> GH(grid.node_number, std::vector<double>(grid.node_number, 0));

    for (auto &element : grid.elements)
    {
        // dla każdego node'a elementu
        for (int row = 0; row < 4; row++){
            for (int col = 0; col < 4; col++){
                int row_node_id = element.node_ids[row]-1;
                int col_node_id = element.node_ids[col]-1;

                GH[row_node_id][col_node_id] += element.H[row][col];
            }
        }
    }

    // printf("\n");
    // for(int i = 0; i < grid.node_number; i++){
    //     for(int j = 0; j < grid.node_number; j++){
    //         printf("%lf ", GH[i][j]);
    //     }
    //     printf("\n");
    // }
}
template<size_t N>
void calculate_HBC(Grid &grid, GlobalData& globalData, const std::array<double, N> &points, int npc = 2){
    
    for(auto& element : grid.elements){
        for(auto& edge : element.bc_edges){
            // dla każdego punjktu całkowania na krawędzi (1D wiec tylko npc a nie npc*npc)
            for(int i = 0; i < npc; i++){
                // waga punktu całkowania
                double w = points[2*i+1];                

                //wyzerowanie macierzyt HBC elemntu
                for(int row = 0; row < 4; row++){
                    for(int col = 0; col < 4; col++){
                        element.HBC[row][col] = 0;
                }}

                for(int row = 0; row < 4; row++){
                    for(int col = 0; col < 4; col++){
                        element.HBC[row][col] += 
                        // arrrr[row][col] += 
                        // printf("edge: %lf %lf %lf %lf \n")
                        UniversalElement4::edges_N_values[edge.edge_index][i][row] *
                        UniversalElement4::edges_N_values[edge.edge_index][i][col] *
                        w *
                        globalData.alpha *
                        element.jacobian.detJ;
                    }
                }
            }
        }
    }
}

template<size_t N>
void calculate_P(Grid &grid, GlobalData& globalData, const std::array<double, N> &points, int npc = 2){
    for(auto& element : grid.elements){
        for(auto& edge : element.bc_edges){
            // dla każdego punjktu całkowania na krawędzi (1D wiec tylko npc a nie npc*npc)
            for(int i = 0; i < npc; i++){
                // waga punktu całkowania
                double w = points[2*i+1];

                // wyzerowanie macierzy (wektora P)
                for(int row = 0; row < 4; row++)
                    element.P[row] = 0;

                for(int row = 0; row < 4; row++){
                    element.P[row] += 
                    UniversalElement4::edges_N_values[edge.edge_index][i][row] *
                    globalData.tot *
                    w *
                    globalData.alpha *
                    element.jacobian.detJ;
                }
            }
        }
    }
}

// agregere local H, C and P; add HBC to H
void agregate(Grid &grid, GlobalData& globalData, EquationData& eqData){
    // set every elemnt of H,C and P to zero
    for (auto& row : eqData.H)
        std::fill(row.begin(), row.end(), 0.0);
    for (auto& row : eqData.C)
        std::fill(row.begin(), row.end(), 0.0);
    std::fill(eqData.P.begin(), eqData.P.end(), 0.0);

    for(auto& element : grid.elements){
        for (int row = 0; row < 4; row++){
            for (int col = 0; col < 4; col++){
                eqData.H[element.node_ids[row]-1][element.node_ids[col]-1] += element.H[row][col];
                // dodanie HBC
                eqData.H[element.node_ids[row]-1][element.node_ids[col]-1] += element.HBC[row][col];

                eqData.C[element.node_ids[row]-1][element.node_ids[col]-1] += element.C[row][col];
            }
            eqData.P[element.node_ids[row]-1] += element.P[row];
        }
    }
}

void print_HG(EquationData ed){
    printf("HG\n");
    for(auto& row : ed.H){
        for(double& el : row){
            printf("%lf ", el);
        }
        printf("\n");
    }
}

void print_P(EquationData ed){
    printf("P:\n");
    for(double& el : ed.P){
        printf("%lf ", el);
    }
    printf("\n");

}


int main(int argc, char const *argv[])
{
    GaussQuad::init_universal_element(GaussQuad::points_2, 2);
    const std::string data_file_path = "./grid_data/Test1_4_4.txt";
    GaussQuad::uniEl.print();
    
    GlobalData global_data;
    Grid grid;
    EquationData equationData;
    load_data_from_file(data_file_path, global_data, grid);
    
    equationData.initMatrixes(global_data.node_number);
    
    
    init_univElem_bc_edges_N_values(GaussQuad::points_2, 2);
    
    global_data.print();
    grid.print();
    
    calculate_jacobian_and_dN_dX_dY(grid);
    
    calculate_H_matrix(grid, global_data, GaussQuad::points_2, 2);
    
    calculate_HG_matrix(grid);
    
    calculate_HBC(grid, global_data, GaussQuad::points_2, 2);

    calculate_P(grid, global_data, GaussQuad::points_2, 2);
    
    agregate(grid, global_data, equationData);
    
    // #######################################################
    
    printf("\n");
    print_HG(equationData);
    print_P(equationData);

    auto ttt = solveLinearSystem(equationData.H, equationData.P);

    for(double &t : ttt){
        printf("%lf ", t);
    }

    return 0;
}
