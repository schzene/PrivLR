#ifndef PRIV_LR_SIGMOID_H__
#define PRIV_LR_SIGMOID_H__

#include "protocol.h"

namespace PrivLR {
    class NonLinear : public Protocol {
    public:
        NonLinear(int party, IOPack *io_pack) : Protocol(party, io_pack) {}
        double mul2add(const double in) const;
        vector<double> mul2add(const vector<double> &in) const;
        double sigmoid(const double in) const;
        vector<double> sigmoid(const vector<double> &in) const;
    };
}

#endif