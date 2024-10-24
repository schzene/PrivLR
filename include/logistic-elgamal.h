#ifndef PRIV_LR_LOGIST_REGRESSION_H_
#define PRIV_LR_LOGIST_REGRESSION_H_

#include "protocols/linear-elgamal.h"
#include "protocols/non-linear-elgamal.h"

namespace PrivLR_Elgamal {
class Logistic {
    IOPack* io_pack       = nullptr;
    Linear* linear        = nullptr;
    NonLinear* non_linear = nullptr;

public:
    vector<double> weight;
    Logistic(int party, IOPack* io_pack, ec_elgamal_secret_key sk, ec_elgamal_public_key pk);
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
}  // namespace PrivLR_Elgamal

#endif