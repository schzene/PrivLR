#ifndef PRIV_LR_SIGMOID_BFV_H__
#define PRIV_LR_SIGMOID_BFV_H__

#include "protocol-bfv.h"

namespace PrivLR_BFV {
class NonLinear : public Protocol {
public:
    NonLinear(const BFVKey* party, const IOPack* io_pack) : Protocol(party, io_pack) {}
    double mul2add(const double in);
    vector<double> mul2add(const vector<double>& in);
    double sigmoid(const double in);
    vector<double> sigmoid(const vector<double>& in);
};
}  // namespace PrivLR_BFV

#endif