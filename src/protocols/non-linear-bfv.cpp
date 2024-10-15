#include "non-linear-bfv.h"
#include "utils/bfv-tools.h"

namespace PrivLR_BFV {
double NonLinear::mul2add(const double in) const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-1, 1);
    std::uniform_real_distribution<> positive_dist1(0, 0.5);
    std::uniform_real_distribution<> positive_dist2(0.5, 1);
    std::uniform_real_distribution<> negative_dist(-0.5, 0);
    double r1, r2;
    if (*party == ALICE) {
        r1           = positive_dist1(gen);
        r2           = positive_dist2(gen);
        double r1_in = r1 * in, r2_in = r2 * in;
        double r1b_inb;
        BFVLongCiphertext r1b_secret_b;
        BFVLongCiphertext::recv(io_pack->io_rev, &r1b_secret_b, party->parm->context);
        io_pack->recv_data(&r1b_inb, sizeof(double));

        uint64_t div_r2_prime = static_cast<uint64_t>(1 / r2 * (1ULL << (2 * SCALE))) % party->parm->plain_mod;
        uint64_t r1_div_r2_prime = static_cast<uint64_t>(r1 / r2 * (1ULL << SCALE)) % party->parm->plain_mod;
        BFVLongPlaintext div_r2_plain(party->parm, div_r2_prime), r1_div_r2_plain(party->parm, r1_div_r2_prime);
        r1b_secret_b.multiply_plain_inplace(r1_div_r2_plain, party->parm->evaluator);
        r1b_secret_b.add_plain_inplace(div_r2_plain, party->parm->evaluator);
        BFVLongCiphertext::send(io_pack->io, &r1b_secret_b);
        io_pack->send_data(&r2_in, sizeof(double));

        return r1b_inb * r1_in;
    } else {
        r1 = negative_dist(gen);
        double r1_in = r1 * in, r2a_ina;
        uint64_t r1_prime = static_cast<uint64_t>(-r1 * (1ULL << SCALE)) % party->parm->plain_mod;
        BFVLongPlaintext r1_plain(party->parm, r1_prime);
        BFVLongCiphertext r1_secret_b(r1_plain, party);
        BFVLongCiphertext::send(io_pack->io, &r1_secret_b);
        io_pack->send_data(&r1_in, sizeof(double));

        BFVLongCiphertext r2_secret_b;
        BFVLongCiphertext::recv(io_pack->io_rev, &r2_secret_b, party->parm->context);
        io_pack->recv_data(&r2a_ina, sizeof(double));
        BFVLongPlaintext r2_plain = r2_secret_b.decrypt(party);
        vector<uint64_t> r2_primes = r2_plain.decode_uint(party->parm);
        uint64_t r2_prime = r2_primes[0];
        r2 = static_cast<double>(r2_prime) / (1ULL << (2 * SCALE));

        return r2 * in * r2a_ina;
    }
    return 0;
}

vector<double> NonLinear::mul2add(const vector<double>& in) const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-1, 1);
    std::uniform_real_distribution<> positive_dist1(0, 0.5);
    std::uniform_real_distribution<> positive_dist2(0.5, 1);
    std::uniform_real_distribution<> negative_dist(-0.5, 0);
    size_t size = in.size();
    vector<double> r1(size), r2(size), res(size);
    if (*party == ALICE) {
        vector<double> r1_in(size), r2_in(size);
        vector<uint64_t> div_r2_prime(size), r1_div_r2_prime(size);
        for (size_t i = 0; i < size; i++) {
            r1[i] = positive_dist1(gen);
            r2[i] = positive_dist2(gen);
            r1_in[i] = r1[i] * in[i];
            r2_in[i] = r2[i] * in[i];
            div_r2_prime[i] = static_cast<uint64_t>(1 / r2[i] * (1ULL << (2 * SCALE))) % party->parm->plain_mod;
            r1_div_r2_prime[i] = static_cast<uint64_t>(r1[i] / r2[i] * (1ULL << SCALE)) % party->parm->plain_mod;
        }
        vector<double> r1b_inb(size);
        BFVLongCiphertext r1b_secret_b;
        BFVLongCiphertext::recv(io_pack->io_rev, &r1b_secret_b, party->parm->context);
        io_pack->recv_data(r1b_inb.data(), sizeof(double) * size);
        
        BFVLongPlaintext div_r2_plain(party->parm, div_r2_prime), r1_div_r2_plain(party->parm, r1_div_r2_prime);
        r1b_secret_b.multiply_plain_inplace(r1_div_r2_plain, party->parm->evaluator);
        r1b_secret_b.add_plain_inplace(div_r2_plain, party->parm->evaluator);
        BFVLongCiphertext::send(io_pack->io, &r1b_secret_b);
        io_pack->send_data(r2_in.data(), sizeof(double) * size);

        for (size_t i = 0; i < size; i++) {
            res[i] = r1b_inb[i] * r1_in[i];
        }
    } else {
        vector<double> r1_in(size), r2a_ina(size);
        vector<uint64_t> r1_prime(size);
        for (size_t i = 0; i < size; i++) {
            r1[i] = negative_dist(gen);
            r1_in[i] = r1[i] * in[i];
            r1_prime[i] = static_cast<uint64_t>(-r1[i] * (1ULL << SCALE)) % party->parm->plain_mod;
        }
        BFVLongPlaintext r1_plain(party->parm, r1_prime);
        BFVLongCiphertext r1_secret_b(r1_plain, party);
        BFVLongCiphertext::send(io_pack->io, &r1_secret_b);
        io_pack->send_data(r1_in.data(), sizeof(double) * size);

        BFVLongCiphertext r2_secret_b;
        BFVLongCiphertext::recv(io_pack->io_rev, &r2_secret_b, party->parm->context);
        io_pack->recv_data(r2a_ina.data(), sizeof(double) * size);
        BFVLongPlaintext r2_plain = r2_secret_b.decrypt(party);
        vector<uint64_t> r2_prime = r2_plain.decode_uint(party->parm);
        for (size_t i = 0; i < size; i++) {
            r2[i] = static_cast<double>(r2_prime[i]) / (1ULL << (2 * SCALE));
            res[i] = r2[i] * in[i] * r2a_ina[i];
        }
    }

    return res;
}

double NonLinear::sigmoid(const double in) const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-1, 1);
    double r = dist(gen), exp_neg_in = exp(-in) * r;
    double r_add = mul2add(r), exp_neg_in_add = mul2add(exp_neg_in);
    double theta = r_add + exp_neg_in_add, theta_remote;
    io_pack->send_data(&theta, sizeof(double));
    io_pack->recv_data(&theta_remote, sizeof(double));
    return r_add / (theta + theta_remote);
}

vector<double> NonLinear::sigmoid(const vector<double>& in) const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-1, 1);
    double size = in.size();
    vector<double> res(size), r(size), exp_neg_in(size), theta(size), theta_remote(size);
    for (size_t i = 0; i < size; i++) {
        r[i]          = dist(gen);
        exp_neg_in[i] = exp(-in[i]) * r[i];
    }
    vector<double> r_add = mul2add(r), exp_neg_in_add = mul2add(exp_neg_in);
    for (size_t i = 0; i < size; i++) {
        theta[i] = r_add[i] + exp_neg_in_add[i];
    }
    io_pack->send_data(theta.data(), sizeof(double) * size);
    io_pack->recv_data(theta_remote.data(), sizeof(double) * size);
    for (size_t i = 0; i < size; i++) {
        r_add[i] = r_add[i] / (theta[i] + theta_remote[i]);
    }
    return r_add;
}
}  // namespace PrivLR_BFV