#ifndef PRIV_LR_PROTOCOL_ELGAMAL_H__
#define PRIV_LR_PROTOCOL_ELGAMAL_H__

#define BIT_LENGTH 37
#define SCALE      12

#include <cassert>
#include <random>
#include <vector>

#include <utils.h>

using std::vector;

namespace PrivLR_Elgamal {
using namespace EC_Elgamal;

class Protocol {
protected:
    ec_elgamal_secret_key sk;
    ec_elgamal_public_key pk;

public:
    int party;
    IOPack* io_pack;

    Protocol(int party, IOPack* io_pack, ec_elgamal_secret_key sk, ec_elgamal_public_key pk) {
        assert(io_pack != nullptr);
        this->party   = party;
        this->io_pack = io_pack;
        this->sk      = sk;
        this->pk      = pk;
    }
};
}  // namespace PrivLR_Elgamal

#endif