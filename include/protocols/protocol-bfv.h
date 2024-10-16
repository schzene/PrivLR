#ifndef PRIV_LR_PROTOCOL_BFV_H__
#define PRIV_LR_PROTOCOL_BFV_H__

#define BIT_LENGTH 37
#define SCALE      12

#include <cassert>
#include <random>
#include <vector>

#include <utils.h>

using std::vector;

namespace PrivLR_BFV {
class Protocol {
protected:
    BFVKey* party;

public:
    IOPack* io_pack;

    Protocol(BFVKey* party, IOPack* io_pack) {
        assert(io_pack != nullptr);
        this->party   = party;
        this->io_pack = io_pack;
    }
};
}  // namespace PrivLR_BFV

#endif