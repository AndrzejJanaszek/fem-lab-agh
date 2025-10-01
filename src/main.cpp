#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace fem{
    struct Node
    {
        int id;

        double x;
        double y;
    };
    
    struct Element4
    {
        int id;

        // todo refs or copies
        std::vector<Node> nodes;
    };
    
    struct Grid
    {
        int node_number;
        int element_number;
        
        // todo refs or copies
        std::vector<Node> nodes;
        std::vector<Element4> elements;
    };

    struct GlobalData
    {
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
    };
    
    void split_str_spaces(std::string line){
        std::vector<std::string> sd;
        int word_index = 0;
        for(char &c : line){
            if(c == ' '){
                word_index++;
            }
            else{

            }
        }
    }

    void load_data_from_file(std::string path){
        // first 10 lines -> global variables [name value]

        
        std::ifstream file(path);
        
        if (!file.is_open()) {
            std::cout << "Error: failed to open file!" << std::endl;
            std::cerr << "Error: failed to open file!" << std::endl;
            exit(-1);
        }

        std::string linia;
        int line_index = 0;
        while (std::getline(file, linia)) {
            

            std::cout << linia << std::endl;
        }

        file.close();
    }
    
    
}

int main(int argc, char const *argv[])
{
    /* code */
    std::cout << "Asdasds";
    return 0;
}
