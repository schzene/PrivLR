#ifndef PRIV_LR_PROTOCOL_H__
#define PRIV_LR_PROTOCOL_H__

#define BIT_LENGTH 37
#define SCALE      12

#include <cassert>
#include <random>
#include <vector>

#include <utils.h>

using std::vector;

namespace PrivLR_Paillier {
using namespace paillier;

class Protocol {
protected:
    paillier::PublicKey pk;
    paillier::PrivateKey sk;

public:
    int party;
    IOPack* io_pack;

    Protocol(int party, IOPack* io_pack) {
        assert(party == ALICE || party == BOB);
        assert(io_pack != nullptr);
        this->party   = party;
        this->io_pack = io_pack;
        keygen(pk, sk);
    }
};
}  // namespace PrivLR_Paillier

#endif