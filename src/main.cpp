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

std::string trim(const std::string & source) {
    std::string s(source);
    s.erase(0,s.find_first_not_of(" \n\r\t"));
    s.erase(s.find_last_not_of(" \n\r\t")+1);
    return s;
}

std::vector<std::string> split_str(std::string line, char delimiter = ' '){
    std::vector<std::string> words;
    std::string tmp = "";

    int word_index = 0;
    for(char &c : line){
        if(c == delimiter){
            // next word
            word_index++;
            words.push_back(tmp);
            tmp = "";
        }
        else{
            tmp += c;
        }
    }

    words.push_back(tmp);

    return words;
}

class Jacobian{
public:
    std::vector<std::array<double, 4>> J;
    std::vector<std::array<double, 4>> J_1;
    double detJ;

    void print(){
        printf("Jacobian:\n");
        printf("dx_dKsi, dy_dKsi, dx_dEta, dy_dEta\n");
        for(auto &el : J){
            // dla każdego punktu całkowania
            printf("%lf, %lf, %lf, %lf\n", el[0], el[1], el[2], el[3]);
        }

        printf("Jacobian Odwrocony:\n");
        printf("dy_dEta, -dy_dKsi, -dx_dEta, dx_dKsi\n");
        for(auto &el : J_1){
            // dla każdego punktu całkowania
            printf("%lf, %lf, %lf, %lf\n", el[0], el[1], el[2], el[3]);
        }

        printf("Wyznacznik Jacobianu:\n");
        printf("%lf:\n", detJ);
    }
};

namespace fem{
    // 2D point
    // id = index-in-array + 1 (counter starting from 1)
    class Node
    {
        public:
        int id;

        double x;
        double y;
        bool isBC = false;

        void print() {
            std::cout << "Node: { ";
            std::cout << "id: " << std::setw(5) << std::right << this->id << ", ";
            std::cout << "x: " << std::setw(10) << std::fixed << std::setprecision(4) << this->x << ", ";
            std::cout << "y: " << std::setw(10) << std::fixed << std::setprecision(4) << this->y;
            std::cout << "isBC: " << std::setw(10) << this->isBC;
            std::cout << " }\n";
        }
    };
    
    // Element made of 4 nodes (stores node ids); 
    // id = index-in-array + 1 (counter starting from 1)
    class Element4
    {
        public:
        
        class Edge{
            public:
            int node_id_1;
            int node_id_2;
            int edge_index;

            Edge(){}

            Edge(int node_id_1,int node_id_2,int edge_index){
                this->node_id_1 = node_id_1;
                this->node_id_2 = node_id_2;
                this->edge_index = edge_index;
            }
        };
        int id;

        int node_ids[4];

        Jacobian jacobian;

        double H[4][4];

        std::vector<Edge> bc_edges;

        std::array<std::array<double, 4>, 4> dN_dX;
        std::array<std::array<double, 4>, 4> dN_dY;

        // todo N1 .. N4 dla punktów na edgach
        // 4 * npc * 4
        // std::vector< std::vector< double > >

        void print(){
            std::cout << "Element4: {\nid: " << this->id << ",\n";
            std::cout << "node_ids: [";
            for(int id : this->node_ids){
                std::cout << id << ", ";
            }
            std::cout << "]\n}\n";
        }
    };
    
    // Grid made of Elements4
    // stores node and elements objects in vectors
    class Grid{
        public:
        int node_number;
        int element_number;
        
        std::vector<Node> nodes;
        std::vector<Element4> elements;

        std::vector<int> border_condition_node_ids;

        std::vector<std::vector<double>> HG;

        void print(){
            std::cout << " --- GRID DATA PRINT --- : " << "\n";
            std::cout << "node_number: " << this->node_number<< "\n";
            std::cout << "element_number: " << this->element_number<< "\n";
            std::cout << "\n";
            std::cout << "NODES LIST: \n";
            for(Node &node : this->nodes){
                node.print();
            }
            std::cout << "\n";
            std::cout << "ELEMENTS LIST: \n";
            for(Element4 &element : this->elements){
                element.print();
            }
            std::cout << "BC NODES: ";
            for(int &id : border_condition_node_ids){
                std::cout << id << ", ";
            }
        }
    };

    // stores simulation physical data & variables
    // and grid size (number of nodes and elements)
    class GlobalData{
        public:
        double simulation_time;
        double simulation_step_time;
        double conductivity;
        double alpha;
        double tot;
        double initial_temperature;
        double density;
        double specific_heat;
        
        int node_number;
        int element_number;

        // array of pointers to double class members
        static constexpr double GlobalData::* fields[] = {
            &GlobalData::simulation_time,
            &GlobalData::simulation_step_time,
            &GlobalData::conductivity,
            &GlobalData::alpha,
            &GlobalData::tot,
            &GlobalData::initial_temperature,
            &GlobalData::density,
            &GlobalData::specific_heat
        };

        void print(){
            std::cout << " --- GLOBAL DATA PRINT --- \n";
            std::cout << "simulation_time: " << this->simulation_time << "\n";
            std::cout << "simulation_step_time: " << this->simulation_step_time << "\n";
            std::cout << "conductivity: " << this->conductivity << "\n";
            std::cout << "alpha: " << this->alpha << "\n";
            std::cout << "tot: " << this->tot << "\n";
            std::cout << "initial_temperature: " << this->initial_temperature << "\n";
            std::cout << "density: " << this->density << "\n";
            std::cout << "specific_heat: " << this->specific_heat << "\n";

            std::cout << "node_number: " << this->node_number << "\n";
            std::cout << "element_number: " << this->element_number << "\n";
        }
    };

    constexpr double GlobalData::* GlobalData::fields[];

    // loads data from file (path) into already initialised objects (GlobalData & Grid)
    void load_data_from_file(std::string path, GlobalData &global_data, Grid &grid){
        // first 10 lines -> global variables [key value]
        std::ifstream file(path);
        
        if (!file.is_open()) {
            std::cout << "Error: failed to open file!" << std::endl;
            std::cerr << "Error: failed to open file!" << std::endl;
            exit(-1);
        }

        std::string line = "";

        // *read 8 doubles
        for(int i = 0; i < 8; i++){
            std::getline(file, line);

            // get value from "key value"
            global_data.*GlobalData::fields[i] = std::stod(split_str(line)[1]);
        }
        
        // *read 2 ints
        std::getline(file, line);
        global_data.node_number = std::stoi(split_str(line)[1]);
        grid.node_number = std::stoi(split_str(line)[1]);
        std::getline(file, line);
        global_data.element_number = std::stoi(split_str(line)[1]);
        grid.element_number = std::stoi(split_str(line)[1]);
        
        // *read Nodes
        std::getline(file, line);   // skip line ("*Node\n")
        grid.nodes.reserve(global_data.node_number);    // reserve space in vector
        for(int i = 0; i < global_data.node_number; i++){
            std::getline(file, line);

            std::vector<std::string> words = split_str(line, ',');

            Node node;
            node.id = stoi(trim(words[0]));
            node.x = stod(trim(words[1]));
            node.y = stod(trim(words[2]));

            grid.nodes.push_back(node);
        }


        // *read Elements
        std::getline(file, line);   // skip line ("*Elements\n")
        grid.elements.reserve(global_data.element_number);   // reserve space in vector
        for(int i = 0; i < global_data.element_number; i++){
            std::getline(file, line);

            std::vector<std::string> words = split_str(line, ',');

            Element4 element;
            element.id = std::stoi(trim(words[0]));
            element.node_ids[0] = std::stoi(trim(words[1]));
            element.node_ids[1] = std::stoi(trim(words[2]));
            element.node_ids[2] = std::stoi(trim(words[3]));
            element.node_ids[3] = std::stoi(trim(words[4]));

            grid.elements.push_back(element);
        }

        // *read Border Conditions
        std::getline(file, line);   // skip line ("*BC\n")

        std::getline(file, line);   // read BC nodes (ids)
        std::vector<std::string> words = split_str(line, ',');
        for(std::string &word : words){
            grid.border_condition_node_ids.push_back( 
                std::stoi( trim(word) )
            );
        }

        // *close file
        file.close();

        // set BC nodes
        for (auto& node : grid.nodes) {
            node.isBC = std::find(
                grid.border_condition_node_ids.begin(),
                grid.border_condition_node_ids.end(),
                node.id
            ) != grid.border_condition_node_ids.end();
        }


        for(auto& element : grid.elements){
            // dla każdej krawędzi
            for(int edge_index = 0; edge_index < 4; edge_index++){
                // czy node 1 i node 2 są BC
                int n1 = element.node_ids[edge_index];
                int n2 = element.node_ids[(edge_index+1)%4];
                bool is_node_1_bc = std::find(grid.border_condition_node_ids.begin(), grid.border_condition_node_ids.end(), n1) != grid.border_condition_node_ids.end();
                bool is_node_2_bc = std::find(grid.border_condition_node_ids.begin(), grid.border_condition_node_ids.end(), n2) != grid.border_condition_node_ids.end();

                // jeżeli oba są => dodaj krawędz do listy krawędzi BC w elemencie
                if(is_node_1_bc && is_node_2_bc){
                    element.bc_edges.push_back(Element4::Edge(n1, n2, edge_index));
                }
            }

        }
    }
    
}

class UniversalElement4{
public:
    double dN_dKsi[4][4] = {0};   // X
    double dN_dEta[4][4] = {0};   // Y


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

class GaussQuad
{
// private:
public:
    static const std::array<double, 4> points_2;
    static const std::array<double, 6> points_3;
    static const std::array<double, 8> points_4;

    static UniversalElement4 uniEl;

    static double dim1_pts2(double (*f)(double)){
        double result = 0;
        for(int i = 0; i < 2; i++){
            // result += f(x_i) * w_i
            result += f(points_2[2*i]) * points_2[2*i+1];
        }

        return result;
    };

    static double dim1_pts3(double (*f)(double)){
        double result = 0;
        for(int i = 0; i < 3; i++){
            // result += f(x_i) * w_i
            result += f(points_3[2*i]) * points_3[2*i+1];
        }
        return result;
    };

    static double dim1_pts4(double (*f)(double)){
        double result = 0;
        for(int i = 0; i < 4; i++){
            // result += f(x_i) * w_i
            result += f(points_4[2*i]) * points_4[2*i+1];
        }
        return result;
    };

    static double dim2_pts2(double (*f)(double, double)){
        double result = 0;
        for(int i = 0; i < 2; i++){
            for(int j = 0; j < 2; j++){
                //              value    *   area
                // result += f(x_i, x_j) * w_i * w_j
                result += f(points_2[2*i], points_2[2*j]) * points_2[2*i+1] * points_2[2*j+1];
            }
        }
        return result;
    };

    static double dim2_pts3(double (*f)(double, double)){
        double result = 0;
        for(int i = 0; i < 3; i++){
            for(int j = 0; j < 3; j++){
                //              value    *   area
                // result += f(x_i, x_j) * w_i * w_j
                result += f(points_3[2*i], points_3[2*j]) * points_3[2*i+1] * points_3[2*j+1];
            }
        }
        return result;
    };

    static double dim2_pts4(double (*f)(double, double)){
        double result = 0;
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++){
                //              value    *   area
                // result += f(x_i, x_j) * w_i * w_j
                result += f(points_4[2*i], points_4[2*j]) * points_4[2*i+1] * points_4[2*j+1];
            }
        }
        return result;
    };

    static void init_universal_element(){
        for(int i_ksi = 0; i_ksi < 2; i_ksi++){
            for(int i_eta = 0; i_eta < 2; i_eta++){
                double ksi = points_2[2*i_ksi];
                double eta = points_2[2*i_eta];
                uniEl.dN_dKsi[2*i_ksi+i_eta][0] = -0.25*(1-eta);
                uniEl.dN_dKsi[2*i_ksi+i_eta][1] = 0.25*(1-eta);
                uniEl.dN_dKsi[2*i_ksi+i_eta][2] = 0.25*(1+eta);
                uniEl.dN_dKsi[2*i_ksi+i_eta][3] = -0.25*(1+eta);
    
                uniEl.dN_dEta[2*i_ksi+i_eta][0] = -0.25*(1-ksi);
                uniEl.dN_dEta[2*i_ksi+i_eta][1] = -0.25*(1+ksi);
                uniEl.dN_dEta[2*i_ksi+i_eta][2] = 0.25*(1+ksi);
                uniEl.dN_dEta[2*i_ksi+i_eta][3] = 0.25*(1-ksi);
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

double f1(double x){
    return 5*pow(x,2) + 3*x + 6;
}

double f2(double x, double y){
    return 5*pow(x,2)*pow(y,2) + 3*x*y + 6;
}

double f_test(double x){
    return 1;
}

// todo refactor
void calculate_jacobian_and_dN_dX_dY(fem::Grid &grid){
    for(auto &element : grid.elements){
        element.jacobian.J.clear();
        element.jacobian.J.reserve(4);
        element.jacobian.J_1.clear();
        element.jacobian.J_1.reserve(4);

        for(int i = 0; i < 4; i++){
            element.jacobian.J.push_back(std::array<double, 4>());
            element.jacobian.J_1.push_back(std::array<double, 4>());
        }

        for(int i_ksi = 0; i_ksi < 2; i_ksi++){
            for(int i_eta = 0; i_eta < 2; i_eta++){
                // dla każdego punktu całkowania [2*i_ksi + i_eta]
                double dx_dKsi = 0;
                double dx_dEta = 0;
                for(int i = 0; i < 4; i++){
                    // x * dN_i/dKsi
                    dx_dKsi += grid.nodes[element.node_ids[i]-1].x * GaussQuad::uniEl.dN_dKsi[2*i_ksi + i_eta][i];
                    // x * dN_i/dEta
                    dx_dEta += grid.nodes[element.node_ids[i]-1].x * GaussQuad::uniEl.dN_dEta[2*i_ksi + i_eta][i];
                }

                double dy_dKsi = 0;
                double dy_dEta = 0;
                for(int i = 0; i < 4; i++){
                    // y * dN_i/dKsi
                    dy_dKsi += grid.nodes[element.node_ids[i]-1].y * GaussQuad::uniEl.dN_dKsi[2*i_ksi + i_eta][i];
                    // y * dN_i/dEta
                    dy_dEta += grid.nodes[element.node_ids[i]-1].y * GaussQuad::uniEl.dN_dEta[2*i_ksi + i_eta][i];
                }

                // jakobian
                element.jacobian.J[2*i_ksi + i_eta][0] = dx_dKsi;
                element.jacobian.J[2*i_ksi + i_eta][1] = dy_dKsi;
                element.jacobian.J[2*i_ksi + i_eta][2] = dx_dEta;
                element.jacobian.J[2*i_ksi + i_eta][3] = dy_dEta;

                // det J
                element.jacobian.detJ = dx_dKsi * dy_dEta - (dy_dKsi * dx_dEta);
                double detJ = dx_dKsi * dy_dEta - (dy_dKsi * dx_dEta);
                // odwórcoy
                element.jacobian.J_1[2*i_ksi + i_eta][0] = dy_dEta / detJ;
                element.jacobian.J_1[2*i_ksi + i_eta][1] = -dy_dKsi / detJ;
                element.jacobian.J_1[2*i_ksi + i_eta][2] = -dx_dEta / detJ;
                element.jacobian.J_1[2*i_ksi + i_eta][3] = dx_dKsi / detJ;

            }
        }
    }

    for(auto &element: grid.elements){
        // dla każdego punktu całkowania
        for(int i_ksi = 0; i_ksi < 2; i_ksi++){
            for(int i_eta = 0; i_eta < 2; i_eta++){
                std::array<double, 4> dNi_dx;
                std::array<double, 4> dNi_dy;

                for(int i = 0; i < 4; i++){
                    // dN/dx
                    double dN_dx = 
                    element.jacobian.J_1[2*i_ksi + i_eta][0] * GaussQuad::uniEl.dN_dKsi[2*i_ksi + i_eta][i] +
                    element.jacobian.J_1[2*i_ksi + i_eta][1] * GaussQuad::uniEl.dN_dEta[2*i_ksi + i_eta][i];
                    dNi_dx[i] = dN_dx;

                    // dN/dy
                    double dN_dy = 
                    element.jacobian.J_1[2*i_ksi + i_eta][2] * GaussQuad::uniEl.dN_dKsi[2*i_ksi + i_eta][i] +
                    element.jacobian.J_1[2*i_ksi + i_eta][3] * GaussQuad::uniEl.dN_dEta[2*i_ksi + i_eta][i];
                    dNi_dy[i] = dN_dy;
                }

                element.dN_dX[2*i_ksi + i_eta] = dNi_dx;
                element.dN_dY[2*i_ksi + i_eta] = dNi_dy;
                // printf("dla %d punktu całkowania:\n", 2*i_ksi + i_eta);
                // printf("dNi/dx:\n");
                // printf("%lf, %lf, %lf, %lf\n", dNi_dx[0], dNi_dx[1], dNi_dx[2], dNi_dx[3]);
                // printf("dNi/dy:\n");
                // printf("%lf, %lf, %lf, %lf\n", dNi_dy[0], dNi_dy[1], dNi_dy[2], dNi_dy[3]);
            }
        }
    }
}

void calculate_H_matrix(fem::Grid &grid, fem::GlobalData &global_data)
{
    for (auto &element : grid.elements)
    {
        double H[4][4] = {0};
        // dla kazdego punktu calkowania
        for (int i_pc = 0; i_pc < 2; i_pc++){
            for (int j_pc = 0; j_pc < 2; j_pc++){

                for (int row = 0; row < 4; row++){
                    for (int col = 0; col < 4; col++){
                        H[row][col] +=
                        (element.dN_dX[2*i_pc + j_pc][row] * element.dN_dX[2*i_pc + j_pc][col] +
                        element.dN_dY[2*i_pc + j_pc][row] * element.dN_dY[2*i_pc + j_pc][col]) *
                        global_data.conductivity * element.jacobian.detJ *
                        GaussQuad::points_2[2*i_pc+1] * GaussQuad::points_2[2*j_pc+1];
                    }
                }
            }
        }

        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++){
                element.H[i][j] = H[i][j];
            }
        }
    }
}

void calculate_HG_matrix(fem::Grid &grid){    
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

    printf("\n");
    for(int i = 0; i < grid.node_number; i++){
        for(int j = 0; j < grid.node_number; j++){
            printf("%lf ", GH[i][j]);
        }
        printf("\n");
    }
}

void calculate_HBC(fem::Grid &grid){
    for(auto& element : grid.elements){
        for(auto& edge : element.bc_edges){
            
        }
    }
}

void calculate_P(fem::Grid &grid, fem::GlobalData &globalData){
    for(auto& element : grid.elements){
        for(auto& edge : element.bc_edges){
            
        }
    }
}

void init_univElem_bc_edges_N_values(){
    // edge 0: node 1 and 2
    // edge 1: node 2 and 3 etc.
    for(int edge_index = 0; edge_index < 4; edge_index++){
        // dla kazdego punktu całkowania
        for(int npc = 0; npc < 2; npc++){
            std::array<double,4> N_arr;
            // wartość dla nie zerowej współrzędnej
            double x_ksi = GaussQuad::points_2[npc*2];
            double y_eta = GaussQuad::points_2[npc*2];
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

int main(int argc, char const *argv[])
{
    printf("test");
    GaussQuad::init_universal_element();
    GaussQuad::uniEl.print();

    fem::GlobalData global_data;
    fem::Grid grid;

    const std::string data_file_path = "./grid_data/Test1_4_4.txt";

    fem::load_data_from_file(data_file_path, global_data, grid);

    init_univElem_bc_edges_N_values();

    global_data.print();
    grid.print();

    calculate_jacobian_and_dN_dX_dY(grid);

    calculate_H_matrix(grid, global_data);

    calculate_HG_matrix(grid);

    
    // for(int i = 0; i < 4; i++){
    //     for(int j = 0; j < 4; j++){
    //         printf("%lf ", grid.elements[0].H[i][j]);
    //     }
    //     printf("\n");
    // }

    return 0;
}
