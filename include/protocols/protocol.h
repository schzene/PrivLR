#ifndef PRIV_LR_PROTOCOL_H__
#define PRIV_LR_PROTOCOL_H__

#define BIT_LENGTH 31

#include <cassert>
#include <random>
#include <vector>

#include <utils.h>

using std::vector;

namespace PrivLR {
    using namespace paillier;

    class Protocol {
    protected:
        int party;
        IOPack *io_pack;
        PublicKey pk;
        PrivateKey sk;

    public:
        Protocol(int party, IOPack *io_pack) {
            assert(party == ALICE || party == BOB);
            assert(io_pack != nullptr);
            this->party = party;
            this->io_pack = io_pack;
            keygen(pk, sk);
        }
    };
}

#endif