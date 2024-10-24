#include "utils/bfv-tools.h"
#include "utils/time-count.h"
#include <protocols/non-linear-bfv.h>
#include <cstdint>

using namespace PrivLR_BFV;

#ifdef USE_TIME_COUNT
timestamp bfv_time     = 0;
timestamp b_start_time = 0;
timestamp b_end_time   = 0;
#endif

class NonLinear_BFV : public Protocol {
public:
    NonLinear_BFV(BFVKey* party, IOPack* io_pack) : Protocol(party, io_pack) {}
    BFVLongCiphertext f(BFVLongCiphertext& __x, BFVLongCiphertext& __x3, BFVLongCiphertext& __x5) const;
    vector<uint64_t> sigmoid(const vector<double>& in) const;
};

uint64_t mod_inverse(uint64_t x, uint64_t mod) {
    int128_t originalMod = mod, temp, quotient;
    int128_t inverse = 0, result = 1;

    if (mod == 1) {
        return 0;
    }

    while (x > 1) {
        quotient = x / mod;
        temp     = mod;
        mod      = x % mod;
        x        = temp;
        temp     = inverse;
        inverse  = result - quotient * inverse;
        result   = temp;
    }

    if (result < 0)
        result += originalMod;

    return static_cast<uint64_t>(result);
}

inline int64_t neg_mod(int64_t val, int64_t mod) {
    return ((val % mod) + mod) % mod;
}

BFVLongCiphertext NonLinear_BFV::f(BFVLongCiphertext& __x1, BFVLongCiphertext& __x3, BFVLongCiphertext& __x5) const {
    uint64_t prime_parm0 = neg_mod(static_cast<int64_t>(0.5 * (1ULL << (2 * SCALE))), party->parm->plain_mod),
             prime_parm1 = neg_mod(static_cast<int64_t>(0.231141280873647 * (1ULL << SCALE)), party->parm->plain_mod),
             prime_parm3 = neg_mod(static_cast<int64_t>(-0.011201270430657 * (1ULL << SCALE)), party->parm->plain_mod),
             prime_parm5 = neg_mod(static_cast<int64_t>(0.000318118716777 * (1ULL << SCALE)), party->parm->plain_mod);
    BFVLongPlaintext parm0(party->parm, prime_parm0), parm1(party->parm, prime_parm1), parm3(party->parm, prime_parm3),
        parm5(party->parm, prime_parm5);
    BFVLongCiphertext __f1 = __x1.multiply_plain(parm1, party->parm->evaluator);
    BFVLongCiphertext __f3 = __x3.multiply_plain(parm3, party->parm->evaluator);
    BFVLongCiphertext __f5 = __x5.multiply_plain(parm5, party->parm->evaluator);
    __f1.add_inplace(__f3, party->parm->evaluator);
    __f1.add_inplace(__f5, party->parm->evaluator);
    __f1.add_plain_inplace(parm0, party->parm->evaluator);
    return __f1;
}

vector<uint64_t> NonLinear_BFV::sigmoid(const vector<double>& in) const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-1, 1);
    std::uniform_real_distribution<> dist1(0, 1);
    size_t size = in.size();
    vector<uint64_t> res_prime(size);
    if (*party == ALICE) {
        vector<double> in_blind(size);
        double r = 1, r3 = r * r * r, r5 = r3 * r * r;
#ifdef USE_TIME_COUNT
        b_start_time = TIME_STAMP;
#endif
        for (size_t i = 0; i < size; i++) {
            in_blind[i] = in[i] * r;
        }
        uint64_t r_prime  = uint64_t(1 / r * (1ULL << SCALE)) % party->parm->plain_mod;
        uint64_t r3_prime = uint64_t(1 / r3 * (1ULL << SCALE)) % party->parm->plain_mod;
        uint64_t r5_prime = uint64_t(1 / r5 * (1ULL << SCALE)) % party->parm->plain_mod;
        BFVLongCiphertext r1_a(party->parm, r_prime, party), r3_a(party->parm, r3_prime, party),
            r5_a(party->parm, r5_prime, party);

        io_pack->send_data(in_blind.data(), sizeof(double) * size);
        BFVLongCiphertext::send(io_pack->io, &r1_a);
        BFVLongCiphertext::send(io_pack->io, &r3_a);
        BFVLongCiphertext::send(io_pack->io, &r5_a);

#ifdef USE_TIME_COUNT
        bfv_time += TIME_STAMP - b_start_time;
#endif

        BFVLongCiphertext ret;
        BFVLongCiphertext::recv(io_pack->io_rev, &ret, party->parm->context);
#ifdef USE_TIME_COUNT
        b_start_time = TIME_STAMP;
#endif
        auto ret_plain = ret.decrypt(party);
        auto ret_prime = ret_plain.decode_uint(party->parm);
#ifdef USE_TIME_COUNT
        bfv_time += TIME_STAMP - b_start_time;
#endif
        return ret_prime;
    }
    else {
        vector<uint64_t> in_prime(size);
        vector<uint64_t> in3_prime(size), in5_prime(size);
        vector<double> in_blind_a(size);
        io_pack->recv_data(in_blind_a.data(), sizeof(double) * size);
        BFVLongCiphertext in1_a, in3_a, in5_a;
        BFVLongCiphertext::recv(io_pack->io_rev, &in1_a, party->parm->context);
        BFVLongCiphertext::recv(io_pack->io_rev, &in3_a, party->parm->context);
        BFVLongCiphertext::recv(io_pack->io_rev, &in5_a, party->parm->context);
#ifdef USE_TIME_COUNT
        b_start_time = TIME_STAMP;
#endif
        for (size_t i = 0; i < size; i++) {
            res_prime[i] = static_cast<uint64_t>(dist1(gen) * (1ULL << (2 * SCALE))) % party->parm->plain_mod;
            in_prime[i]  = static_cast<uint64_t>(in[i] * in_blind_a[i] * (1ULL << SCALE)) % party->parm->plain_mod;
            in3_prime[i] = static_cast<uint64_t>(in[i] * in[i] * in[i] * in_blind_a[i] * in_blind_a[i] * in_blind_a[i] *
                                                 (1ULL << SCALE)) %
                           party->parm->plain_mod;
            in5_prime[i] = static_cast<uint64_t>(in[i] * in[i] * in[i] * in[i] * in[i] * in_blind_a[i] * in_blind_a[i] *
                                                 in_blind_a[i] * in_blind_a[i] * in_blind_a[i] * (1ULL << SCALE)) %
                           party->parm->plain_mod;
        }
        BFVLongPlaintext scale_inv(party->parm, mod_inverse(1ULL << SCALE, party->parm->plain_mod)),
            res_p(party->parm, res_prime), in1_plain(party->parm, in_prime), in3_plain(party->parm, in3_prime),
            in5_plain(party->parm, in5_prime);
        in1_a.multiply_plain_inplace(in1_plain, party->parm->evaluator);
        in3_a.multiply_plain_inplace(in3_plain, party->parm->evaluator);
        in5_a.multiply_plain_inplace(in5_plain, party->parm->evaluator);
        in1_a.multiply_plain_inplace(scale_inv, party->parm->evaluator);
        in3_a.multiply_plain_inplace(scale_inv, party->parm->evaluator);
        in5_a.multiply_plain_inplace(scale_inv, party->parm->evaluator);
        BFVLongCiphertext res_sec_a = f(in1_a, in3_a, in5_a);
        res_sec_a.sub_plain_inplace(res_p, party->parm->evaluator);

        BFVLongCiphertext::send(io_pack->io, &res_sec_a);
#ifdef USE_TIME_COUNT
        bfv_time += TIME_STAMP - b_start_time;
#endif
    }
    return res_prime;
}

int party_;

void test_mpc(NonLinear* non_linear, const vector<double>& data) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-3, 3);
    size_t size = data.size();
    vector<double> in(size);
    if (party_ == ALICE) {
        vector<double> in_b(size);
        for (size_t i = 0; i < size; i++) {
            in[i]   = dist(gen);
            in_b[i] = data[i] - in[i];
        }
        non_linear->io_pack->send_data(in_b.data(), sizeof(double) * size);
    }
    else {
        non_linear->io_pack->recv_data(in.data(), sizeof(double) * size);
    }

    vector<double> res = non_linear->sigmoid(in), res_remote(size);

    // non_linear->io_pack->send_data(res.data(), sizeof(double) * size);
    // non_linear->io_pack->recv_data(res_remote.data(), sizeof(double) * size);
    // std::cout << "error:\n";
    // for (size_t i = 0; i < size; i++) {
    //     std::cout << res[i] + res_remote[i] - 1 / (1 + exp(-data[i])) << "\n";
    // }
}

void test_bfv(NonLinear_BFV* non_linear, const vector<double>& data, uint64_t plain_mod) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-3, 3);
    size_t size = data.size();
    vector<double> in(size), in_b(size);
    if (party_ == ALICE) {
        for (size_t i = 0; i < size; i++) {
            in[i]   = dist(gen);
            in_b[i] = data[i] / in[i];
        }
        non_linear->io_pack->send_data(in_b.data(), sizeof(double) * size);
    }
    else {
        non_linear->io_pack->recv_data(in.data(), sizeof(double) * size);
    }

    vector<uint64_t> res_prime = non_linear->sigmoid(in);

    // vector<double> res(size), res_remote(size);
    // for (size_t i = 0; i < size; i++) {
    //     res[i] = res_prime[i];
    //     if (res[i] > plain_mod / 2) {
    //         res[i] -= plain_mod;
    //     }
    //     res[i] *= (1. / (1ULL << (2 * SCALE)));
    // }
    // non_linear->io_pack->send_data(res.data(), sizeof(double) * size);
    // non_linear->io_pack->recv_data(res_remote.data(), sizeof(double) * size);
    // std::cout << "error:\n";
    // for (size_t i = 0; i < size; i++) {
    //     std::cout << res[i] + res_remote[i] - 1 / (1 + exp(-data[i])) << "\n";
    // }
}

int main(int argc, const char** argv) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-3, 3);

    size_t size = 15;
    size_t round = 25;

    party_ = argv[1][0] - '0';
    if (party_ == ALICE) {
        std::cout << "Party: ALICE"
                  << "\n";
    }
    else {
        party_ = BOB;
        std::cout << "Party: BOB"
                  << "\n";
    }
    BFVParm* parm   = new BFVParm(8192, default_prime_mod.at(29));
    BFVKey* party   = new BFVKey(party_, parm);
    IOPack* io_pack = new IOPack(party_);

    NonLinear* non_linear_mpc     = new NonLinear(party, io_pack);
    NonLinear_BFV* non_linear_bfv = new NonLinear_BFV(party, io_pack);
    vector<double> in(size);
    if (party_ == ALICE) {
        for (size_t i = 0; i < size; i++) {
            in[i] = dist(gen);
        }
        io_pack->send_data(in.data(), sizeof(double) * size);
    }
    else {
        io_pack->recv_data(in.data(), sizeof(double) * size);
    }

    for (size_t i = 0; i < round; i++) {
        test_mpc(non_linear_mpc, in);
        test_bfv(non_linear_bfv, in, party->parm->plain_mod);
    }
    std::cout << "mpc time count: " << non_linear_time << "\n";
    std::cout << "bfv time count: " << bfv_time << "\n";

    delete io_pack;
    delete non_linear_mpc;
}