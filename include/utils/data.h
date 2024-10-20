#ifndef PRIV_LR_DATA_H__
#define PRIV_LR_DATA_H__

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using std::string;
using std::vector;

void load_dataset(vector<vector<double>>& datas, vector<int>& label, const string& filename);

#endif