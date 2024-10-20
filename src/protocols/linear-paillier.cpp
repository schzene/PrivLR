#include "linear-paillier.h"
#include "utils/paillier.h"
#include <iostream>

namespace PrivLR_Paillier {
double Linear::dot_product(const vector<double>& in_a, const vector<double>& in_b) const {
    assert(in_a.size() == in_b.size());
    size_t size = in_a.size();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-1, 1);

    double res = 0;
    vector<ZZ> in_a_ring(size), in_b_ring(size);
    for (size_t i = 0; i < size; i++) {
        res += in_a[i] * in_b[i];
        in_a_ring[i] = ZZ(static_cast<uint64_t>(in_a[i] * (1ULL << SCALE)) % (1ULL << BIT_LENGTH));
        in_b_ring[i] = ZZ(static_cast<uint64_t>(in_b[i] * (1ULL << SCALE)) % (1ULL << BIT_LENGTH));
    }
    if (party == ALICE) {
        CipherVector in_a_secret_a(in_a_ring, pk), in_b_secret_a(in_b_ring, pk);
        CipherVector::send(io_pack, &in_a_secret_a);
        CipherVector::send(io_pack, &in_b_secret_a);

        size_t p1_size, p2_size;
        string res_part1_sec_str, res_part2_sec_str;
        io_pack->recv_data(&p1_size, sizeof(size_t));
        res_part1_sec_str.resize(p1_size);
        io_pack->recv_data(&res_part1_sec_str[0], sizeof(char) * p1_size);

        io_pack->recv_data(&p2_size, sizeof(size_t));
        res_part2_sec_str.resize(p2_size);
        io_pack->recv_data(&res_part2_sec_str[0], sizeof(char) * p2_size);

        ZZ res_part1_sec_zz = string2ZZ(res_part1_sec_str), res_part2_sec_zz = string2ZZ(res_part2_sec_str);
        paillier::Ciphertext res_part1_sec(res_part1_sec_zz), res_part2_sec(res_part2_sec_zz);
        uint64_t res_part1_ring = NTL::conv<uint64_t>(decrypt(res_part1_sec, sk) % (1ULL << BIT_LENGTH)),
                 res_part2_ring = NTL::conv<uint64_t>(decrypt(res_part2_sec, sk) % (1ULL << BIT_LENGTH));
        int64_t res_part1 = (res_part1_ring >= (1ULL << BIT_LENGTH) / 2) ? (res_part1_ring - (1ULL << BIT_LENGTH)) :
                                                                           res_part1_ring,
                res_part2 = (res_part2_ring >= (1ULL << BIT_LENGTH) / 2) ? (res_part2_ring - (1ULL << BIT_LENGTH)) :
                                                                           res_part2_ring;
        res += (res_part1 + res_part2) * 1. / (1ULL << (2 * SCALE));
    }
    else {
        CipherVector in_a_secret_a, in_b_secret_a;
        CipherVector::recv(io_pack, &in_a_secret_a);
        CipherVector::recv(io_pack, &in_b_secret_a);

        in_a_secret_a.mul_plain_inplace(in_b_ring);
        in_b_secret_a.mul_plain_inplace(in_a_ring);

        ZZ res_part1_sec_a = in_a_secret_a.data[0], res_part2_sec_a = in_b_secret_a.data[0],
           n_square = in_a_secret_a.public_key.n * in_a_secret_a.public_key.n;
        for (size_t i = 1; i < size; i++) {
            res_part1_sec_a = (res_part1_sec_a * in_a_secret_a.data[i]) % n_square;
            res_part2_sec_a = (res_part2_sec_a * in_b_secret_a.data[i]) % n_square;
        }
        double res_part1 = dist(gen), res_part2 = dist(gen);
        ZZ res_part1_ring = ZZ(static_cast<uint64_t>(res_part1 * (1ULL << (2 * SCALE))) % (1ULL << BIT_LENGTH)),
           res_part2_ring = ZZ(static_cast<uint64_t>(res_part2 * (1ULL << (2 * SCALE))) % (1ULL << BIT_LENGTH));

        paillier::Ciphertext res_part1_sec = encrypt(res_part1_ring, in_a_secret_a.public_key),
                             res_part2_sec = encrypt(res_part2_ring, in_a_secret_a.public_key);
        res_part1_sec_a                    = (res_part1_sec_a * res_part1_sec.data) % n_square;
        res_part2_sec_a                    = (res_part2_sec_a * res_part2_sec.data) % n_square;
        res -= (res_part1 + res_part2);

        string res_part1_sec_a_str = ZZ2string(res_part1_sec_a), res_part2_sec_a_str = ZZ2string(res_part2_sec_a);
        size_t p1_size = res_part1_sec_a_str.size(), p2_size = res_part2_sec_a_str.size();
        io_pack->send_data(&p1_size, sizeof(size_t));
        io_pack->send_data(res_part1_sec_a_str.c_str(), sizeof(char) * p1_size);
        io_pack->send_data(&p2_size, sizeof(size_t));
        io_pack->send_data(res_part2_sec_a_str.c_str(), sizeof(char) * p2_size);
    }
    return res;
}

vector<double> Linear::dot_product(const vector<vector<double>>& in_a, const vector<double>& in_b,
                                   double transpose) const {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-1, 1);
    size_t data_size = transpose ? in_a[0].size() : in_a.size(), size = in_b.size();
    vector<double> in_a_flatten(data_size * size);
    if (transpose) {
        for (size_t i = 0; i < data_size; i++) {
            for (size_t j = 0; j < size; j++) {
                in_a_flatten[i * size + j] = in_a[j][i];
            }
        }
    }
    else {
        for (size_t i = 0; i < data_size; i++) {
            for (size_t j = 0; j < size; j++) {
                in_a_flatten[i * size + j] = in_a[i][j];
            }
        }
    }
    vector<double> res(data_size);
    vector<ZZ> in_a_flatten_ring(data_size * size), in_b_ring(data_size * size);
    for (size_t j = 0; j < size; j++) {
        for (size_t i = 0; i < data_size; i++) {
            res[i] += in_a_flatten[i * size + j] * in_b[j];
            in_a_flatten_ring[i * size + j] =
                ZZ(static_cast<uint64_t>(in_a_flatten[i * size + j] * (1ULL << SCALE)) % (1ULL << BIT_LENGTH));
        }
        in_b_ring[j] = ZZ(static_cast<uint64_t>(in_b[j] * (1ULL << SCALE)) % (1ULL << BIT_LENGTH));
    }

    if (party == ALICE) {
        CipherVector in_a_secret_a(in_a_flatten_ring, pk), in_b_secret_a(in_b_ring, pk);
        CipherVector::send(io_pack, &in_a_secret_a);
        CipherVector::send(io_pack, &in_b_secret_a);

        CipherVector res_part1_secret_a, res_part2_secret_a;
        CipherVector::recv(io_pack, &res_part1_secret_a);
        CipherVector::recv(io_pack, &res_part2_secret_a);
        vector<ZZ> res_part1_ring = res_part1_secret_a.decrypt(sk), res_part2_ring = res_part2_secret_a.decrypt(sk);
        vector<int64_t> res_part1(data_size), res_part2(data_size);
        for (size_t i = 0; i < data_size; i++) {
            res_part1_ring[i] = res_part1_ring[i] % (1ULL << BIT_LENGTH);
            res_part2_ring[i] = res_part2_ring[i] % (1ULL << BIT_LENGTH);
            res_part1[i]      = NTL::conv<int64_t>((res_part1_ring[i] >= ((1ULL << BIT_LENGTH) / 2)) ?
                                                       (res_part1_ring[i] - ((1ULL << BIT_LENGTH))) :
                                                       res_part1_ring[i]);
            res_part2[i]      = NTL::conv<int64_t>((res_part2_ring[i] >= ((1ULL << BIT_LENGTH) / 2)) ?
                                                       (res_part2_ring[i] - ((1ULL << BIT_LENGTH))) :
                                                       res_part2_ring[i]);
            res[i] += (res_part1[i] + res_part2[i]) * 1. / (1ULL << (2 * SCALE));
        }
    }
    else {
        CipherVector in_a_secret_a, in_b_sec_a_compressed, in_b_secret_a;
        CipherVector::recv(io_pack, &in_a_secret_a);
        CipherVector::recv(io_pack, &in_b_sec_a_compressed);

        in_b_secret_a.public_key = in_b_sec_a_compressed.public_key;
        in_b_secret_a.data.resize(data_size * size);
        vector<ZZ> in_b_flatten_ring(data_size * size);
        for (size_t i = 0; i < data_size; i++) {
            for (size_t j = 0; j < size; j++) {
                in_b_flatten_ring[i * size + j]  = in_b_ring[j];
                in_b_secret_a.data[i * size + j] = in_b_sec_a_compressed.data[j];
            }
        }
        in_a_secret_a.mul_plain_inplace(in_b_flatten_ring);
        in_b_secret_a.mul_plain_inplace(in_a_flatten_ring);

        vector<double> res_part1(data_size * size), res_part2(data_size * size);
        vector<ZZ> res_part1_ring(data_size * size), res_part2_ring(data_size * size);
        for (size_t i = 0; i < data_size * size; i++) {
            res_part1[i]      = dist(gen);
            res_part2[i]      = dist(gen);
            res_part1_ring[i] = ZZ(static_cast<uint64_t>(res_part1[i] * (1ULL << (2 * SCALE))) % (1ULL << BIT_LENGTH));
            res_part2_ring[i] = ZZ(static_cast<uint64_t>(res_part2[i] * (1ULL << (2 * SCALE))) % (1ULL << BIT_LENGTH));
        }
        in_a_secret_a.add_plain_inplace(res_part1_ring);
        in_b_secret_a.add_plain_inplace(res_part2_ring);
        for (size_t i = 0; i < data_size; i++) {
            for (size_t j = 0; j < size; j++) {
                res[i] -= (res_part1[i * size + j] + res_part2[i * size + j]);
            }
        }
        ZZ n_square = in_a_secret_a.public_key.n * in_a_secret_a.public_key.n;
        CipherVector res_part1_secret_a, res_part2_secret_a;
        res_part1_secret_a.public_key = in_a_secret_a.public_key;
        res_part2_secret_a.public_key = in_b_secret_a.public_key;
        res_part1_secret_a.data       = vector<ZZ>(data_size, ZZ(1));
        res_part2_secret_a.data       = vector<ZZ>(data_size, ZZ(1));
        for (size_t i = 0; i < data_size; i++) {
            for (size_t j = 0; j < size; j++) {
                res_part1_secret_a.data[i] = (res_part1_secret_a.data[i] * in_a_secret_a.data[i * size + j]) % n_square;
                res_part2_secret_a.data[i] = (res_part2_secret_a.data[i] * in_b_secret_a.data[i * size + j]) % n_square;
            }
        }
        CipherVector::send(io_pack, &res_part1_secret_a);
        CipherVector::send(io_pack, &res_part2_secret_a);
    }
    return res;
}

// vector<double> Linear::dot_product(const vector<vector<double>>& in_a, const vector<double>& in_b,
//                                    double transpose) const {
//     size_t data_size = transpose ? in_a[0].size() : in_a.size(), size = in_b.size();
//     vector<vector<double>> in_a_t;
//     if (transpose) {
//         in_a_t = vector<vector<double>>(data_size, vector<double>(size));
//         for (size_t i = 0; i < data_size; i++) {
//             for (size_t j = 0; j < size; j++) {
//                 in_a_t[i][j] = in_a[j][i];
//             }
//         }
//     }

//     vector<double> res(data_size);
//     if (transpose) {
//         for (size_t i = 0; i < data_size; i++) {
//             res[i] = dot_product(in_a_t[i], in_b);
//         }
//     } else {
//         for (size_t i = 0; i < data_size; i++) {
//             res[i] = dot_product(in_a[i], in_b);
//         }
//     }
//     return res;
// }
}  // namespace PrivLR_Paillier