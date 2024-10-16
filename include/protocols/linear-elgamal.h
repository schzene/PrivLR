#ifndef PRIV_LR_MUL_ELGAMAL_H__
#define PRIV_LR_MUL_ELGAMAL_H__

#include "protocol-elgamal.h"

namespace PrivLR_Elgamal {
class Linear : public Protocol {
public:
    Linear(int party, IOPack* io_pack, ec_elgamal_secret_key sk, ec_elgamal_public_key pk)
        : Protocol(party, io_pack, sk, pk) {}
    double dot_product(const vector<double>& in_a, const vector<double>& in_b) const;
    vector<double> dot_product(const vector<vector<double>>& in_a, const vector<double>& in_b,
                               double transpose = false) const;
};
}  // namespace PrivLR_Elgamal

#endif