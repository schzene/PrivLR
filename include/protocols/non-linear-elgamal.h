#ifndef PRIV_LR_SIGMOID_H__
#define PRIV_LR_SIGMOID_H__

#include "protocol-elgamal.h"

namespace PrivLR_Elgamal {
class NonLinear : public Protocol {
public:
    NonLinear(int party, IOPack* io_pack, ec_elgamal_secret_key sk, ec_elgamal_public_key pk)
        : Protocol(party, io_pack, sk, pk) {}
    double mul2add(const double in) const;
    vector<double> mul2add(const vector<double>& in) const;
    double sigmoid(const double in) const;
    vector<double> sigmoid(const vector<double>& in) const;
};
}  // namespace PrivLR_Elgamal

#endif