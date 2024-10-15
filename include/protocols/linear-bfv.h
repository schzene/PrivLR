#ifndef PRIV_LR_MUL_H__
#define PRIV_LR_MUL_H__

#include "protocol-bfv.h"

namespace PrivLR_BFV {
    class Linear : public Protocol {
    public:
        Linear(BFVKey* party, IOPack *io_pack) : Protocol(party, io_pack) {}
        double dot_product(const vector<double> &in_a, const vector<double> &in_b) const;
        vector<double> dot_product(const vector<vector<double>> &in_a, const vector<double> &in_b, double transpose = false) const;
    };
}

#endif