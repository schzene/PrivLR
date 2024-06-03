#include "data.h"

void load_dataset(vector<vector<double>> &datas, vector<int> &label, const string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }
    string line;
    while (std::getline(file, line)) {
        std::istringstream record(line);
        vector<double> data;
        data.push_back(0.5);
        double temp;
        while (record >> temp) {
            data.push_back(temp);
        }
        label.push_back(int(temp));
        data.pop_back();
        datas.push_back(data);
    }
}