#define USE_TIME_COUNT
#include <protocols/non-linear-bfv.h>
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
    if (party->party == ALICE) {
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
        io_pack->send_data(in_blind.data(), sizeof(double) * size);
#ifdef USE_TIME_COUNT
        bfv_time += TIME_STAMP - b_start_time;
#endif

        BFVLongCiphertext r1_a(r_prime, party), r3_a(r3_prime, party), r5_a(r5_prime, party);
        BFVLongCiphertext::send(io_pack->io, &r1_a);
        BFVLongCiphertext::send(io_pack->io, &r3_a);
        BFVLongCiphertext::send(io_pack->io, &r5_a);

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

void Logistic::gradAscent(vector<vector<double>>& datas, vector<int>& label, int max_cycles, const double alpha) {
    const size_t data_size = datas.size(), size = datas[0].size();
    weight = vector<double>(size, .5);
    while (max_cycles > 0) {
        vector<double> classified = linear->dot_product(datas, weight);
        vector<double> h          = non_linear->sigmoid(classified);
        vector<double> error(data_size);
        start_time = TIME_STAMP;
        for (int i = 0; i < data_size; i++) {
            double dist = label[i] - h[i];
            if (abs(dist) < 1e-10) {
                dist = 0;
            }
            error[i] = dist;
        }
        time_cost += TIME_STAMP - start_time;
        vector<double> delta_weight = linear->dot_product(datas, error, true);
        start_time                  = TIME_STAMP;
        for (int i = 0; i < size; i++) {
            weight[i] += alpha * delta_weight[i];
        }
        max_cycles--;
        time_cost += TIME_STAMP - start_time;
        // Not a protocol content, only for statistical purposes
        // {
        //     printf("Cycle remain: %3d", max_cycles);
        //     double sum_error = 0.;
        //     vector<double> h_remote(data_size);
        //     vector<int> label_remote(data_size);
        //     io_pack->send_data(h.data(), sizeof(double) * data_size, false);
        //     io_pack->send_data(label.data(), sizeof(int) * data_size, false);
        //     io_pack->recv_data(h_remote.data(), sizeof(double) * data_size, false);
        //     io_pack->recv_data(label_remote.data(), sizeof(int) * data_size, false);
        //     for (size_t i = 0; i < data_size; i++) {
        //         sum_error += -1 * (label[i] + label_remote[i]) * log(h[i] + h_remote[i]) -
        //                      (1 - label[i] - label_remote[i]) * log(1 - h[i] - h_remote[i]);
        //     }
        //     printf(", loss: %.10lf\n", sum_error / data_size);
        // }
    }
    time_cost += time_cost + linear->time_cost + non_linear->time_cost;
}

int party_;
string ip = "127.0.0.1";