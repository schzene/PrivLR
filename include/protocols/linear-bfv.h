#ifndef PRIV_LR_MUL_BFV_H__
#define PRIV_LR_MUL_BFV_H__

#include "protocol-bfv.h"

#ifdef USE_TIME_COUNT
timestamp linear_time = 0;
timestamp l_start_time = 0;
timestamp l_end_time = 0;
#endif

namespace PrivLR_BFV {
class Linear : public Protocol {
public:
    Linear(BFVKey* party, IOPack* io_pack) : Protocol(party, io_pack) {}
    double dot_product(const vector<double>& in_a, const vector<double>& in_b) const;
    vector<double> dot_product(const vector<vector<double>>& in_a, const vector<double>& in_b,
                               double transpose = false) const;
};
}  // namespace PrivLR_BFV

#endif