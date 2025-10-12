#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <string>
#include <iomanip>
#include <cmath>
#include <map>

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

namespace fem{
    // 2D point
    // id = index-in-array + 1 (counter starting from 1)
    class Node
    {
        public:
        int id;

        double x;
        double y;

        void print() {
            std::cout << "Node: { ";
            std::cout << "id: " << std::setw(5) << std::right << this->id << ", ";
            std::cout << "x: " << std::setw(10) << std::fixed << std::setprecision(4) << this->x << ", ";
            std::cout << "y: " << std::setw(10) << std::fixed << std::setprecision(4) << this->y;
            std::cout << " }\n";
        }
    };
    
    // Element made of 4 nodes (stores node ids); 
    // id = index-in-array + 1 (counter starting from 1)
    class Element4
    {
        public:
        int id;

        int node_ids[4];

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
    }
    
}

// gauss quadrature
// 1 dimension
// 2 points
// domain [-1; 1]
inline double gauss_quad_1d_2p(double (*f)(double)){
    // 1/sqrt(3)
    constexpr double x0 = 0.57735026918962584;
    // weights = 1
    return (f(x0) + f(-x0));
}

// gauss quadrature
// 1 dimension
// 3 points
// domain [-1; 1]
inline double gauss_quad_1d_3p(double (*f)(double)){
    // sqrt(3/5)
    constexpr double x0 = 0.77459666924148340;
    constexpr double x1 = 0;
    // weights
    constexpr double w0 = 5.0 / 9.0;
    constexpr double w1 = 8.0 / 9.0;

    return w0*(f(x0) + f(-x0)) + w1*f(x1);
}



inline double gauss_quad_2d_2p(double (*f)(double, double)){
    // 1/sqrt(3)
    constexpr double x0 = 0.57735026918962584;
    // weights = 1
    return  f(x0, x0) + 
            f(-x0, -x0) + 
            f(-x0, x0) + 
            f(x0, -x0);
}

const double gauss_quad_points_2[] = {
    -0.57735026918962584, 1, 
    0.57735026918962584, 1
};

const double gauss_quad_points_3[] = {
    -0.77459666924148340, (5.0 / 9.0),
    0, 8.0 / 9.0,
    0.77459666924148340, (5.0 / 9.0)
};

const double* gauss_quad_point_arrays[]={
    gauss_quad_points_2,
    gauss_quad_points_3
};

class GaussQuad
{
private:
    static const std::array<double, 4> points_2;
    static const std::array<double, 6> points_3;
public:
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


double f1(double x){
    return 5*pow(x,2) + 3*x + 6;
}

double f2(double x, double y){
    return 5*pow(x,2)*pow(y,2) + 3*x*y + 6;
}

double f_test(double x){
    return 1;
}

int main(int argc, char const *argv[])
{

    // fem::GlobalData global_data;
    // fem::Grid grid;

    // const std::string data_file_path = "./grid_data/Test1_4_4.txt";

    // fem::load_data_from_file(data_file_path, global_data, grid);

    // global_data.print();
    // grid.print();

    // printf("%lf\n", GaussQuad::dim1_pts2(f_test));
    printf("%lf\n", GaussQuad::dim1_pts2(f1));
    printf("%lf\n", GaussQuad::dim1_pts3(f1));
    printf("%lf\n", GaussQuad::dim2_pts2(f2));
    printf("%lf\n", GaussQuad::dim2_pts3(f2));

    return 0;
}
