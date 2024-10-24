#ifndef PRIV_LR_PROTOCOL_BFV_H__
#define PRIV_LR_PROTOCOL_BFV_H__

#define SCALE 12

#include <cassert>
#include <random>
#include <vector>

#include <utils.h>

using std::vector;

namespace PrivLR_BFV {
class Protocol {
protected:
    const BFVKey* party;

public:
    const IOPack* io_pack;
    timestamp time_cost  = 0;
    timestamp start_time = 0;

    Protocol(const BFVKey* _party, const IOPack* _io_pack) : party(_party), io_pack(_io_pack) {
        assert(party != nullptr);
        assert(io_pack != nullptr);
    }
};
}  // namespace PrivLR_BFV

#endif