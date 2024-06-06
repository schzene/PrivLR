#include <protocols/non-linear.h>
using namespace PrivLR;

int main(int argc, const char **argv) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-1, 1);
    std::uniform_real_distribution<> positive_dist(0, 1);
    size_t size = 10;

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

    NonLinear *non_linear = new NonLinear(party_, io_pack);
    vector<double> in(size), in_remote(size);
    for (size_t i = 0; i < size; i++) {
        in[i] = dist(gen);
    }

    vector<double> res = non_linear->sigmoid(in), res_remote(size);
    io_pack->send_data(res.data(), sizeof(double) * size);
    io_pack->recv_data(res_remote.data(), sizeof(double) * size);
    io_pack->send_data(in.data(), sizeof(double) * size);
    io_pack->recv_data(in_remote.data(), sizeof(double) * size);
    std::cout << "error:\n";
    for (size_t i = 0; i < size; i++) {
        std::cout << res[i] + res_remote[i] - 1 / (1 + exp(-in[i] - in_remote[i])) << "\n";
    }

    delete io_pack;
    delete non_linear;
}