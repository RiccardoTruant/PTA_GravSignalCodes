#include <iostream>
#include <fstream>
#include <string>
using namespace std;
typedef double my_type;
typedef vector<my_type> T1;//tensor 1 indices -> vector



int countLinesInFile(const string& path, const string& model_name) {
    int N_sources = 0; // Line counter
    string line;
    ifstream read_model(path + model_name);

    // Check if the file was opened successfully
    if (!read_model.is_open()) {
        cerr << "Error: file " << path + model_name << " not found!" << endl;
        exit(-1);
    }

    // Count the lines in the file
    while (getline(read_model, line)) {
        N_sources++;
    }

    read_model.close();

    // Output the result
    return N_sources;
}






