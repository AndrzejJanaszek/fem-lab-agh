#pragma once

#include <vector>
#include <array>

class GaussQuad;

class UniversalElement4{
public:
    std::vector<std::array<double, 4>>dN_dKsi;   // X
    std::vector<std::array<double, 4>>dN_dEta;   // Y
    // double dN_dKsi[4][4] = {0};   // X
    // double dN_dEta[4][4] = {0};   // Y


    static std::array< std::vector< std::array<double,4> >,4>edges_N_values;

    // edge_index reprezentuje globalny numer krawędzi
    // 0: ma node pierwszy i drugi
    // 1: ma node drugi i trzeci itd.

    // class Surface{
    //     public:
    //     // todo dla zdefiniowanej liczby npc
    //     // dla 2 punktów całkowania
    //     std::array<std::vector<double>, 2> N_values_for_npc;
        
        
    // }

    // key: id edge'a, key: id punktu calkowania, wektor Ni dla danego punktu całkowania na krawędzi
    // static std::map<int, std::map<int, std::vector<double>>> N_values_for_edge_points;

    void print(){
        printf("dN_dKsi:\n");
        for(int i = 0; i < 4; i++)
            printf("%lf, %lf, %lf, %lf\n", dN_dKsi[i][0], dN_dKsi[i][1], dN_dKsi[i][2], dN_dKsi[i][3]);

        printf("dN_dEta:\n");
        for(int i = 0; i < 4; i++)
            printf("%lf, %lf, %lf, %lf\n", dN_dEta[i][0], dN_dEta[i][1], dN_dEta[i][2], dN_dEta[i][3]);
    }
};

std::array< std::vector<std::array<double,4>>, 4 > UniversalElement4::edges_N_values = {};

class GaussQuad
{
public:
    static const std::array<double, 4> points_2;
    static const std::array<double, 6> points_3;
    static const std::array<double, 8> points_4;

    static UniversalElement4 uniEl;

    template<size_t N>
    double integrate1D(double (*f)(double), const std::array<double, N> &points) {
        static_assert(N % 2 == 0, "Tablica z wartościami punktów całkowania musi być podzielna przez dwa (wrtość, waga)");

        double result = 0;
        constexpr int k = N / 2;

        for(int i = 0; i < k; i++)
            result += f(points[2*i]) * points[2*i+1];

        return result;
    }

    template<size_t N>
    double integrate2D(double (*f)(double, double), const std::array<double, N> &points) {
        static_assert(N % 2 == 0, "Tablica z wartościami punktów całkowania musi być podzielna przez dwa (wrtość, waga)");

        double result = 0;
        constexpr int k = N / 2;

        for(int i = 0; i < k; i++)
            for(int j = 0; j < k; j++)
                result += f(points[2*i], points[2*j]) *
                        points[2*i+1] * points[2*j+1];

        return result;
    }

    template<size_t N>
    static void init_universal_element(const std::array<double, N> &points, int npc = 2){
        // init vectors
        uniEl.dN_dKsi.clear();
        uniEl.dN_dEta.clear();
        uniEl.dN_dKsi.resize(npc * npc);
        uniEl.dN_dEta.resize(npc * npc);
        
        // set data
        for(int i_ksi = 0; i_ksi < 2; i_ksi++){
            for(int i_eta = 0; i_eta < 2; i_eta++){
                double ksi = points[2*i_ksi];
                double eta = points[2*i_eta];
                uniEl.dN_dKsi[npc*i_ksi+i_eta][0] = -0.25*(1-eta);
                uniEl.dN_dKsi[npc*i_ksi+i_eta][1] = 0.25*(1-eta);
                uniEl.dN_dKsi[npc*i_ksi+i_eta][2] = 0.25*(1+eta);
                uniEl.dN_dKsi[npc*i_ksi+i_eta][3] = -0.25*(1+eta);
    
                uniEl.dN_dEta[npc*i_ksi+i_eta][0] = -0.25*(1-ksi);
                uniEl.dN_dEta[npc*i_ksi+i_eta][1] = -0.25*(1+ksi);
                uniEl.dN_dEta[npc*i_ksi+i_eta][2] = 0.25*(1+ksi);
                uniEl.dN_dEta[npc*i_ksi+i_eta][3] = 0.25*(1-ksi);
            }
        }
    }
};

const std::array<double, 4> GaussQuad::points_2 = {
        -0.57735026918962584, 1, 
        0.57735026918962584, 1
};

const std::array<double, 6> GaussQuad::points_3 = {
    -0.77459666924148340, (5.0 / 9.0),
    0, 8.0 / 9.0,
    0.77459666924148340, (5.0 / 9.0)
};

const std::array<double, 8> GaussQuad::points_4 = {
    0.861136, 0.347855,
    0.339981, 0.652145,
    -0.339981, 0.652145,
    -0.861136, 0.347855
};

UniversalElement4 GaussQuad::uniEl;

template<size_t N>
void init_univElem_bc_edges_N_values(const std::array<double, N> &points, int npc = 2){
    // edge 0: node 1 and 2
    // edge 1: node 2 and 3 etc.
    for(int edge_index = 0; edge_index < 4; edge_index++){
        // dla kazdego punktu całkowania
        for(int i = 0; i < npc; i++){
            std::array<double,4> N_arr;
            // wartość dla nie zerowej współrzędnej
            double x_ksi = points[i*2];
            double y_eta = points[i*2];
            // dla edge_index 0 i 2 zerują się wartości pionowe (y)
            // dla edge_index 1 i 3 zerują sie wartości horyzontalne (x)
            if(edge_index % 2 == 0){
                y_eta = 0;
            }
            else{
                x_ksi = 0;
            }

            N_arr[0] = 0.25*(1-x_ksi)*(1-y_eta);
            N_arr[1] = 0.25*(1+x_ksi)*(1-y_eta);
            N_arr[2] = 0.25*(1+x_ksi)*(1+y_eta);
            N_arr[3] = 0.25*(1-x_ksi)*(1+y_eta);

            UniversalElement4::edges_N_values[edge_index].push_back(N_arr);
        }
    }
}