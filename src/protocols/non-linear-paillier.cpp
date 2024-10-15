#include "non-linear-paillier.h"

namespace PrivLR_Paillier {
double NonLinear::mul2add(const double in) const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-1, 1);
    std::uniform_real_distribution<> positive_dist1(0, 0.5);
    std::uniform_real_distribution<> positive_dist2(0.5, 1);
    std::uniform_real_distribution<> negative_dist(-0.5, 0);
    double r1, r2;
    if (party == ALICE) {
        r1           = positive_dist1(gen);
        r2           = positive_dist2(gen);
        double r1_in = r1 * in, r2_in = r2 * in;
        double r1b_inb;
        size_t pkb_str_size;
        io_pack->recv_data(&pkb_str_size, sizeof(size_t));
        string pkb_str;
        pkb_str.resize(pkb_str_size);
        io_pack->recv_data(pkb_str.data(), sizeof(char) * pkb_str_size);
        size_t r1b_str_size;
        io_pack->recv_data(&r1b_str_size, sizeof(size_t));
        string r1b_str;
        r1b_str.resize(r1b_str_size);
        io_pack->recv_data(r1b_str.data(), sizeof(char) * r1b_str_size);
        io_pack->recv_data(&r1b_inb, sizeof(double));
        paillier::PublicKey pkb(string2ZZ(pkb_str));
        paillier::Ciphertext r1b_secret_b(string2ZZ(r1b_str));

        ZZ div_r2_fixed    = ZZ(uint64_t((1 / r2) * (1ULL << (2 * SCALE))) % (1ULL << BIT_LENGTH));
        ZZ r1_div_r2_fixed = ZZ(uint64_t((r1 / r2) * (1ULL << SCALE)) % (1ULL << BIT_LENGTH));
        mul_plain_inplace(r1b_secret_b, r1_div_r2_fixed, pkb);
        add_plain_inplace(r1b_secret_b, div_r2_fixed, pkb);
        string ct_str      = ZZ2string(r1b_secret_b.data);
        size_t ct_str_size = ct_str.size();
        io_pack->send_data(&ct_str_size, sizeof(size_t));
        io_pack->send_data(ct_str.data(), sizeof(char) * ct_str_size);
        io_pack->send_data(&r2_in, sizeof(double));

        return r1b_inb * r1_in;
    }
    else {
        r1                               = negative_dist(gen);
        double r1_in                     = r1 * in, r2a_ina;
        ZZ r1_fixed                      = ZZ(uint64_t((-r1) * (1ULL << SCALE)));
        paillier::Ciphertext r1_secret_b = encrypt(r1_fixed, pk);
        string pk_str                    = ZZ2string(pk.n);
        size_t pk_str_size               = pk_str.size();
        string ct_str                    = ZZ2string(r1_secret_b.data);
        size_t ct_str_size               = ct_str.size();
        io_pack->send_data(&pk_str_size, sizeof(size_t));
        io_pack->send_data(pk_str.data(), sizeof(char) * pk_str_size);
        io_pack->send_data(&ct_str_size, sizeof(size_t));
        io_pack->send_data(ct_str.data(), sizeof(char) * ct_str_size);
        io_pack->send_data(&r1_in, sizeof(double));

        size_t r2_str_size;
        io_pack->recv_data(&r2_str_size, sizeof(size_t));
        string r2_str;
        r2_str.resize(r2_str_size);
        io_pack->recv_data(r2_str.data(), sizeof(char) * r2_str_size);
        io_pack->recv_data(&r2a_ina, sizeof(double));
        paillier::Ciphertext r2_secret_b(string2ZZ(r2_str));
        ZZ r2_fixed = decrypt(r2_secret_b, sk);
        r2          = NTL::to_double(r2_fixed) / (1ULL << (2 * SCALE));

        return r2 * in * r2a_ina;
    }
}

vector<double> NonLinear::mul2add(const vector<double>& in) const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-1, 1);
    std::uniform_real_distribution<> positive_dist1(0, 0.5);
    std::uniform_real_distribution<> positive_dist2(0.5, 1);
    std::uniform_real_distribution<> negative_dist(-0.5, 0);
    size_t size = in.size();
    if (party == ALICE) {
        vector<double> r1(size), r2(size), r1_in(size), r2_in(size), r1b_inb(size);
        vector<ZZ> div_r2_fixed(size), r1_div_r2_fixed(size);
        CipherVector r1b_secret_b;
        for (size_t i = 0; i < size; i++) {
            r1[i]              = positive_dist1(gen);
            r2[i]              = positive_dist2(gen);
            r1_in[i]           = r1[i] * in[i];
            r2_in[i]           = r2[i] * in[i];
            div_r2_fixed[i]    = ZZ(uint64_t((1 / r2[i]) * (1ULL << (2 * SCALE))) % (1ULL << BIT_LENGTH));
            r1_div_r2_fixed[i] = ZZ(uint64_t((r1[i] / r2[i]) * (1ULL << (SCALE))) % (1ULL << BIT_LENGTH));
        }
        CipherVector::recv(io_pack, &r1b_secret_b);
        io_pack->recv_data(r1b_inb.data(), sizeof(double) * size);

        r1b_secret_b.mul_plain_inplace(r1_div_r2_fixed);
        r1b_secret_b.add_plain_inplace(div_r2_fixed);
        CipherVector::send(io_pack, &r1b_secret_b);
        io_pack->send_data(r2_in.data(), sizeof(double) * size);

        for (size_t i = 0; i < size; i++) {
            r1b_inb[i] *= r1_in[i];
        }
        return r1b_inb;
    }
    else {
        vector<double> r1(size), r1_in(size), r2a_ina(size);
        vector<ZZ> r1_fixed(size);
        CipherVector r2_secret_b;
        for (size_t i = 0; i < size; i++) {
            r1[i]       = negative_dist(gen);
            r1_in[i]    = r1[i] * in[i];
            r1_fixed[i] = ZZ(uint64_t((-r1[i]) * (1ull << SCALE)) % (1ULL << BIT_LENGTH));
        }
        CipherVector r1_secret_b(r1_fixed, pk);
        CipherVector::send(io_pack, &r1_secret_b);
        io_pack->send_data(r1_in.data(), sizeof(double) * size);

        CipherVector::recv(io_pack, &r2_secret_b);
        io_pack->recv_data(r2a_ina.data(), sizeof(double) * size);
        vector<ZZ> r2_fixed = r2_secret_b.decrypt(sk);

        for (size_t i = 0; i < size; i++) {
            r2a_ina[i] = r2a_ina[i] * in[i] * NTL::to_double(r2_fixed[i] % (1ULL << BIT_LENGTH)) / (1ULL << (2 * SCALE));
        }

        return r2a_ina;
    }
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
}  // namespace PrivLR_Paillier