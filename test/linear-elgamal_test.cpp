#include <protocols/linear-elgamal.h>
#include "utils/ec-elgamal.h"
#include "utils/ec-elgamal.hpp"
#include "utils/time-count.h"

using namespace PrivLR_Elgamal;

int main(int argc, const char** argv) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-1, 1);
    std::uniform_real_distribution<> positive_dist(0, 1);
    int party_ = argv[1][0] - '0';
    if (party_ == ALICE) {
        std::cout << "Party: ALICE"
                  << "\n";
    }
    else {
        party_ = BOB;
        std::cout << "Party: BOB"
                  << "\n";
    }
    IOPack* io_pack = new IOPack(party_);

    ec_elgamal_init(DEFAULT_CURVE);
    bsgs_table_t bsgs_table;
    bsgs_table_init(bsgs_table, 1ULL << (MAX_BITS / 2));
    ec_elgamal_secret_key sk = ec_elgamal_keygen();
    ec_elgamal_public_key pk = ec_elgamal_from_secret_key(sk);

    size_t size      = 10;
    size_t data_size = 20;
    vector<vector<double>> in_a(size, vector<double>(data_size));
    vector<double> in_b(size);
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < data_size; j++) {
            in_a[i][j] = dist(gen);
        }
        in_b[i] = dist(gen);
    }
    Linear* linear = new Linear(party_, io_pack, sk, pk);
    auto res       = linear->dot_product(in_a, in_b, true);

    if (party_ == ALICE) {
        for (size_t i = 0; i < in_a.size(); i++) {
            io_pack->send_data(in_a[i].data(), sizeof(double) * data_size);
        }
        io_pack->send_data(in_b.data(), sizeof(double) * size);
        io_pack->send_data(res.data(), sizeof(double) * data_size);
    }
    else {
        vector<double> res_a(data_size), true_res(data_size);
        vector<vector<double>> in_a_a(size, vector<double>(data_size));
        vector<double> in_b_a(size);
        for (size_t i = 0; i < in_a.size(); i++) {
            io_pack->recv_data(in_a_a[i].data(), sizeof(double) * data_size);
        }
        io_pack->recv_data(in_b_a.data(), sizeof(double) * size);
        io_pack->recv_data(res_a.data(), sizeof(double) * data_size);
        for (size_t j = 0; j < data_size; j++) {
            for (size_t i = 0; i < size; i++) {
                true_res[j] += (in_a[i][j] + in_a_a[i][j]) * (in_b[i] + in_b_a[i]);
            }
        }
        std::cout << "true_res: "
                  << "\n";
        for (size_t j = 0; j < data_size; j++) {
            cout << true_res[j] << "\n";
        }

        std::cout << "res: "
                  << "\n";
        for (size_t j = 0; j < data_size; j++) {
            cout << res[j] + res_a[j] << "\n";
        }
    }

    // size_t size = 2;
    // vector<uint64_t> a(size);
    // for (size_t i = 0; i < size; i++) {
    //     a[i] = uint64_t(positive_dist(gen) * (1ULL << 24));
    // }
    // CipherVector cta(a, pk);

    // INIT_TIMER
    // START_TIMER
    // vector<uint64_t> res = cta.decrypt(sk, bsgs_table);
    // STOP_TIMER("elgamal dcrypt")

    // for(size_t i = 0; i < size; i++) {
    //     std::cout << res[i] << "\n";
    // }

    ec_elgamal_free_public_key(pk);
    ec_elgamal_free_secret_key(sk);
    bsgs_table_free(bsgs_table);
    ec_elgamal_deinit();
    delete io_pack;
}