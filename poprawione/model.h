#pragma once

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


#include "utils.h"


class Jacobian{
public:
    std::vector<std::array<double, 4>> J;
    std::vector<std::array<double, 4>> J_1;
    std::vector<double> detJ;

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
    double HBC[4][4];
    double C[4][4];
    double P[4];

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

    class EquationData {
    public:
        std::vector<std::vector<double>> H;
        std::vector<std::vector<double>> C;
        std::vector<double> P;

        void initMatrixes(int n) {
            H.assign(n, std::vector<double>(n, 0.0));
            C.assign(n, std::vector<double>(n, 0.0));
            P.assign(n, 0.0);
        }
    };