#ifndef PRIV_LR_LOGIST_REGRESSION_H_
#define PRIV_LR_LOGIST_REGRESSION_H_

#include "protocols/linear-bfv.h"
#include "protocols/non-linear-bfv.h"

namespace PrivLR_BFV {
class Logistic {
public:
    timestamp time_cost  = 0;
    timestamp start_time = 0;
    IOPack* io_pack;
    Linear* linear;
    NonLinear* non_linear;

    vector<double> weight;
    Logistic(BFVKey* party, IOPack* io_pack);
    ~Logistic();
    void gradAscent(vector<vector<double>>& datas, vector<int>& label, int max_cycles = 500,
                    const double alpha = 0.001);

    inline double classify(const vector<double>& data) const {
        return non_linear->sigmoid(linear->dot_product(data, weight));
    }

    inline vector<double> classify(const vector<vector<double>>& datas) const {
        return non_linear->sigmoid(linear->dot_product(datas, weight));
    }
};
}  // namespace PrivLR_BFV

#endif