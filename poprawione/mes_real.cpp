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
#include <sstream>


#include "model.h"
#include "gauss.h"
#include "solver.h"
#include "vtu_parser.h"


// miedź
const double COPPER_CONDUCTIVITY = 400.0;
const double COPPER_DENSITY = 8940.0;
const double COPPER_SPECIFIC_HEAT = 385.0;
// powietrze
const double AIR_CONDUCTIVITY = 0.025;
const double AIR_DENSITY = 1.25;
const double AIR_SPECIFIC_HEAT = 1004.0;

const double HEAT_GENERATION = 65 * 1000000; // W/m3

// IHS I RESZTA RADIATORA
std::vector<int> COPPER_ELEMENTS = {95,96,97,101,102,103,107,108,109,113,114,115,125,126,127,131,132,133,137,138,139,143,144,145,155,156,157,161,162,163,167,168,169,173,174,175,185,186,187,191,192,193,197,198,199,203,204,205,215,216,217,221,222,223,227,228,229,233,234,235,245,246,247,251,252,253,257,258,259,263,264,265,275,276,277,281,282,283,287,288,289,293,294,295,305,306,307,311,312,313,317,318,319,323,324,325,335,336,337,341,342,343,347,348,349,353,354,355,365,366,367,371,372,373,377,378,379,383,384,385,395,396,397,401,402,403,407,408,409,413,414,415,425,426,427,431,432,433,437,438,439,443,444,445,455,456,457,461,462,463,467,468,469,473,474,475,485,486,487,491,492,493,497,498,499,503,504,505,515,516,517,521,522,523,527,528,529,533,534,535,545,546,547,551,552,553,557,558,559,563,564,565,575,576,577,581,582,583,587,588,589,593,594,595,605,606,607,611,612,613,617,618,619,623,624,625,635,636,637,641,642,643,647,648,649,653,654,655,665,666,667,671,672,673,677,678,679,683,684,685,695,696,697,701,702,703,707,708,709,713,714,715,725,726,727,731,732,733,737,738,739,743,744,745,755,756,757,761,762,763,767,768,769,773,774,775,785,786,787,791,792,793,797,798,799,803,804,805,811,812,813,814,815,816,817,818,819,820,821,822,823,824,825,826,827,828,829,830,831,832,833,834,835,836,837,838,839,840,841,842,843,844,845,846,847,848,849,850,851,852,853,854,855,856,857,858,859,860,861,862,863,864,865,866,867,868,869,870,871,872,873,874,875,876,877,878,879,880,881,882,883,884,885,886,887,888,889,890,891,892,893,894,895,896,897,898,899,900};
// std::vector<int> COPPER_ELEMENTS = {1
// };

// CPU DIE - CZESC PROCESORA, POLPRZEWODNIK
std::vector<int> INITIAL_HOT_ELEMENTS = {881,882,883,884,885,886,887,888,889,890};
// std::vector<int> INITIAL_HOT_ELEMENTS = {
// 404,405,406,434,435,436,464,465,466
// };

// założenia programu:
// element 4 węzłowy
// liczba możliwych punktów całkowania: 2,3,4

void calculate_jacobian_and_dN_dX_dY(Grid &grid, int npc){
    for(auto &element : grid.elements){
        element.jacobian.J.resize(npc*npc);
        element.jacobian.J_1.resize(npc*npc);
        element.jacobian.detJ.resize(npc*npc);

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
                element.jacobian.detJ[npc*i_ksi + i_eta] = dx_dKsi * dy_dEta - (dy_dKsi * dx_dEta);
                double detJ = element.jacobian.detJ[npc*i_ksi + i_eta];
                // odwórcoy
                element.jacobian.J_1[npc*i_ksi + i_eta][0] = dy_dEta / detJ;
                element.jacobian.J_1[npc*i_ksi + i_eta][1] = -dy_dKsi / detJ;
                element.jacobian.J_1[npc*i_ksi + i_eta][2] = -dx_dEta / detJ;
                element.jacobian.J_1[npc*i_ksi + i_eta][3] = dx_dKsi / detJ;
            }
        }
    }

    for(auto &element: grid.elements){
        element.dN_dX.resize(npc*npc);
        element.dN_dY.resize(npc*npc);

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

// CALCULATE H, C FOR EACH ELEMENT
template<size_t N>
void calculate_H_and_C_matrix(Grid &grid, GlobalData &global_data, const std::array<double, N> &points, int npc)
{
    for (auto &element : grid.elements)
    {
        double H[4][4] = {0};
        double C[4][4] = {0};

        double mat_conductivity = AIR_CONDUCTIVITY;
        double mat_specific_heat = AIR_SPECIFIC_HEAT;
        double mat_density = AIR_DENSITY;

        if(std::find(COPPER_ELEMENTS.begin(), COPPER_ELEMENTS.end(), element.id) != COPPER_ELEMENTS.end()){
            mat_conductivity = global_data.conductivity;
            mat_specific_heat = global_data.specific_heat;
            mat_density = global_data.density;
        }

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
                        mat_conductivity * element.jacobian.detJ[npc*i_pc + j_pc] *
                        points[2*i_pc+1] * points[2*j_pc+1];
                        
                        // macierz C
                        C[row][col] +=
                        (N_arr[row] * N_arr[col]) *
                        mat_specific_heat * mat_density *
                        element.jacobian.detJ[npc*i_pc + j_pc] *
                        points[2*i_pc+1] * points[2*j_pc+1];
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

// CALCULATE FOR EACH ELEMENT HBC
template<size_t N>
void calculate_HBC(Grid &grid, GlobalData& globalData, const std::array<double, N> &points, int npc){
    
    for(auto& element : grid.elements){
        //wyzerowanie macierzyt HBC elemntu
        for(int row = 0; row < 4; row++){
            for(int col = 0; col < 4; col++){
                element.HBC[row][col] = 0;
        }}

        for(auto& edge : element.bc_edges){
            // dla każdego punjktu całkowania na krawędzi (1D wiec tylko npc a nie npc*npc)
            for(int i = 0; i < npc; i++){
                // waga punktu całkowania
                double w = points[2*i+1];                

                
                double dx = grid.nodes[edge.node_id_2-1].x - grid.nodes[edge.node_id_1-1].x;
                double dy = grid.nodes[edge.node_id_2-1].y - grid.nodes[edge.node_id_1-1].y;
                double detJ = sqrt(dx*dx + dy*dy)/2.0;
                for(int row = 0; row < 4; row++){
                    for(int col = 0; col < 4; col++){
                        element.HBC[row][col] += 
                        // arrrr[row][col] += 
                        // printf("edge: %lf %lf %lf %lf \n")
                        UniversalElement4::edges_N_values[edge.edge_index][i][row] *
                        UniversalElement4::edges_N_values[edge.edge_index][i][col] *
                        w *
                        globalData.alpha *
                        detJ;
                    }
                }
            }
        }
    }
}

// calculate C for each element 
template<size_t N>
void calculate_P(Grid &grid, GlobalData& globalData, const std::array<double, N> &points, int npc){
    for(auto& element : grid.elements){
        // wyzerowanie macierzy (wektora P)
        for(int row = 0; row < 4; row++)
            element.P[row] = 0;

        for(auto& edge : element.bc_edges){
            // dla każdego punjktu całkowania na krawędzi (1D wiec tylko npc a nie npc*npc)
            for(int i = 0; i < npc; i++){
                // waga punktu całkowania
                double w = points[2*i+1];

                // dla każdego N
                double dx = grid.nodes[edge.node_id_2-1].x - grid.nodes[edge.node_id_1-1].x;
                double dy = grid.nodes[edge.node_id_2-1].y - grid.nodes[edge.node_id_1-1].y;
                double detJ = sqrt(dx*dx + dy*dy)/2.0;
                for(int row = 0; row < 4; row++){
                    element.P[row] += 
                    UniversalElement4::edges_N_values[edge.edge_index][i][row] *
                    globalData.tot *
                    w *
                    globalData.alpha *
                    detJ;
                    // element.jacobian.detJ;
                }
            }
        }

        // WEWNETRZNE XRODLO CIEPLA Q
        if(std::find(INITIAL_HOT_ELEMENTS.begin(), INITIAL_HOT_ELEMENTS.end(), element.id) != INITIAL_HOT_ELEMENTS.end()){
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
                        element.P[row] +=
                        N_arr[row] *
                        HEAT_GENERATION *
                        element.jacobian.detJ[npc*i_pc + j_pc] *
                        points[2*i_pc+1] * points[2*j_pc+1];
                    }
                }
            }
        }

    }
}

// agregate local H, C and P; add HBC to H
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

void agregate_time_part(Grid &grid, GlobalData& globalData, EquationData& eqData, std::vector<double> t0){
    for (int row = 0; row < grid.node_number; row++){
        for (int col = 0; col < grid.node_number; col++){
            // dodajemy C/dt
            
            // eqData.H[row][col] += eqData.C[row][col] * (1.0/ 50.0);
            eqData.H[row][col] += eqData.C[row][col] / globalData.simulation_step_time;
        }

        double tmp = 0; // C/dTeta * t0
        for(int i = 0; i < globalData.node_number; i++){
            tmp += eqData.C[row][i] * t0[i]; 
        }
        tmp /= globalData.simulation_step_time;

        eqData.P[row] += tmp;
    }
}

void print_H_from_elements(Grid grid){
    std::vector<std::vector<double>> H;
    H.assign(grid.node_number, std::vector<double>(grid.node_number, 0.0));

    for(auto& element : grid.elements){
        for (int row = 0; row < 4; row++){
            for (int col = 0; col < 4; col++){
                H[element.node_ids[row]-1][element.node_ids[col]-1] += element.H[row][col];
            }
        }
    }

    printf("H:\n");
    for(auto& row : H){
        for(double& el : row){
            printf("%lf ", el);
        }
        printf("\n");
    }
}

void print_C_from_elements(Grid grid){
    std::vector<std::vector<double>> H;
    H.assign(grid.node_number, std::vector<double>(grid.node_number, 0.0));

    for(auto& element : grid.elements){
        for (int row = 0; row < 4; row++){
            for (int col = 0; col < 4; col++){
                H[element.node_ids[row]-1][element.node_ids[col]-1] += element.C[row][col];
            }
        }
    }

    printf("C:\n");
    for(auto& row : H){
        for(double& el : row){
            printf("%lf ", el);
        }
        printf("\n");
    }
}

void print_P_from_elements(Grid grid){
    std::vector<double> P;
    P.assign(grid.node_number, 0.0);

    for(auto& element : grid.elements){
        for (int row = 0; row < 4; row++){
            P[element.node_ids[row]-1] += element.P[row];
        }
    }

    printf("P:\n");
    for(double& el : P){
        printf("%lf ", el);
    }
    printf("\n");
}



void print_H(EquationData ed){
    printf("H:\n");
    for(auto& row : ed.H){
        for(double& el : row){
            printf("%lf ", el);
        }
        printf("\n");
    }
}

void print_C(EquationData ed){
    printf("C:\n");
    for(auto& row : ed.C){
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
    const int GAUSS_I_POINTS = 2;
    auto GAUSS_POINTS_ARRAY = GaussQuad::points_2;

    // std::string RESULT_PATH = "results/1_4_4/";
    // std::string RESULT_PATH = "results/1_4_4_mix/";
    // std::string RESULT_PATH = "results/31_31/";
    std::string RESULT_PATH = "results/real_1/";
    // ###############################################################
    // #                            INIT
    // ###############################################################
    GaussQuad::init_universal_element(GAUSS_POINTS_ARRAY, GAUSS_I_POINTS);
    // const std::string data_file_path = "./grid_data/Test1_4_4.txt";
    // const std::string data_file_path = "./grid_data/Test2_4_4_MixGrid.txt";
    // const std::string data_file_path = "./grid_data/Test3_31_31_kwadrat.txt";
    const std::string data_file_path = "./real_problem/real_4_kolumny.txt";
    GaussQuad::uniEl.print(GAUSS_I_POINTS);
    
    GlobalData global_data;
    Grid grid;
    EquationData equationData;
    load_data_from_file(data_file_path, global_data, grid);
    
    equationData.initMatrixes(global_data.node_number);
    
    init_univElem_bc_edges_N_values(GAUSS_POINTS_ARRAY, GAUSS_I_POINTS);
    
    calculate_jacobian_and_dN_dX_dY(grid, GAUSS_I_POINTS);

    // ###############################################################
    // #                            SIM
    // ###############################################################

    // petla symulacji

    //* SETUP INITIAL TEMPERATURE
    // std::vector<double> temperature_v_initial = std::vector<double>(global_data.node_number,global_data.initial_temperature);
    std::vector<double> temperature_v_initial = std::vector<double>(global_data.node_number, global_data.tot);

    std::vector<double> temperature_v = std::vector<double>(global_data.node_number, 0);
    printf("step: initial\n");
        for(double &t : temperature_v_initial){
            printf("%lf ", t);
        }
    printf("\n");

    // global_data.simulation_step_time = 100;
    
    global_data.simulation_time = 40;
    global_data.simulation_step_time = 1;

    // global_data.print();
    int step = 0;
    for(int stime = 0; stime <= global_data.simulation_time; stime+=global_data.simulation_step_time){
        // // CAŁY CZAS USTALA TEMPERATURE OD PROCESORA
        // for(int& el_i : INITIAL_HOT_ELEMENTS){
        //     for(int& node_id : grid.elements[el_i-1].node_ids){
        //         temperature_v_initial[node_id-1] = 100;
        //     }
        // }

        calculate_H_and_C_matrix(grid, global_data, GAUSS_POINTS_ARRAY, GAUSS_I_POINTS);
        // print_H_from_elements(grid);
        // print_C_from_elements(grid);

        calculate_HBC(grid, global_data, GAUSS_POINTS_ARRAY, GAUSS_I_POINTS);

        calculate_P(grid, global_data, GAUSS_POINTS_ARRAY, GAUSS_I_POINTS);
        // print_P_from_elements(grid);

        agregate(grid, global_data, equationData);

        // print_H(equationData);

        
        // print_C(equationData);
        
        agregate_time_part(grid, global_data, equationData, temperature_v_initial);

        // #######################
        // printf("Po agregacji z czasem\n");
        // print_H(equationData);
        // print_P(equationData);
        // #######################

        
        temperature_v = solveLinearSystem(equationData.H, equationData.P);
        
        // printf("step: %d\n", stime);
        // for(double &t : temperature_v){
        //     printf("%lf ", t);
        // }
        // printf("\n");

        temperature_v_initial = temperature_v;

        // //* JEDNA ITERACJA BREAK
        // break;

        //* ZAPIS SYMULACJI DO PLIKU - ZAPIS KORKU SYMULACJI
        std::stringstream ss;
        ss << RESULT_PATH
        << "step_"
        << std::setw(3) << std::setfill('0') << step
        << ".vtu";

        writeVTU(ss.str(), grid, temperature_v);
        step++;

        printf("step: %d\n", step);
    }

    //* ZAPIS SYMULACJI DO PLIKU - OPIS KROKÓW SYMULACJI
    writePVD(
        RESULT_PATH + "simulation.pvd",
        step,
        global_data.simulation_step_time
    );

    return 0;
}
