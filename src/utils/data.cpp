#include "data.h"

void load_dataset(vector<vector<double>> &datas, vector<int> &label, const string &filename) {
    std::ifstream file(filename);
    string line;
    while (std::getline(file, line)) {
        std::istringstream record(line);
        vector<double> data;
        data.push_back(1.0);
        double temp;
        while (record >> temp) {
            data.push_back(temp);
        }
        label.push_back(int(temp));
        data.pop_back();
        datas.push_back(data);
    }
}