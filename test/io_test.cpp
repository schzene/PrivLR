#include <utils.h>

int main(int argc, const char **argv) {
    int party_ = argv[1][0] - '0';
    if (party_ == ALICE) {
        std::cout << "Party: ALICE"
                  << "\n";
    } else {
        party_ = BOB;
        std::cout << "Party: BOB"
                  << "\n";
    }
    IOPack *io_pack = new IOPack(party_);
    if (party_ == ALICE) {
        ZZ n(100);
        // string n_str = NTL::
        io_pack->send_data(&n, sizeof(ZZ));
    } else {
        ZZ n;
        io_pack->recv_data(&n, sizeof(ZZ));
        std::cout << n << "\n";
    }
}