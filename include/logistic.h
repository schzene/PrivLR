#ifndef PRIV_LR_LOGIST_REGRESSION_H_
#define PRIV_LR_LOGIST_REGRESSION_H_

#include <algorithm>

#include "protocols/linear.h"
#include "protocols/non-linear.h"

namespace PrivLR {
    class Logistic {
        Linear *linear = nullptr;
        NonLinear *non_linear = nullptr;
        vector<double> weight;

    public:
        Logistic(int party, IOPack *io_pack);
        ~Logistic();
        void gradAscent(vector<vector<double>> &datas, vector<int> &label,
                        int max_cycles = 500, const double alpha = 0.001);
        void stocGradAscent(vector<vector<double>> &datas, vector<int> &labelMat, 
                            int num_iter = 150, double alpha = 0.01);

        inline double classify(const vector<double> &data) const {
            return linear->dot_product(data, weight);
        }

        inline vector<double> classify(const vector<vector<double>> &datas) const {
            return linear->dot_product(datas, weight);
        }
    };
}

#endif