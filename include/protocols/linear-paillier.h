#ifndef PRIV_LR_MUL_H__
#define PRIV_LR_MUL_H__

#include "protocol-paillier.h"

namespace PrivLR_Paillier {
class Linear : public Protocol {
public:
    Linear(int party, IOPack* io_pack) : Protocol(party, io_pack) {}
    double dot_product(const vector<double>& in_a, const vector<double>& in_b) const;
    vector<double> dot_product(const vector<vector<double>>& in_a, const vector<double>& in_b,
                               double transpose = false) const;
};
}  // namespace PrivLR_Paillier

#endif