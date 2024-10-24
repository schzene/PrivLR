#include "protocols/linear-elgamal.h"

namespace PrivLR_Elgamal {
double Linear::dot_product(const vector<double>& in_a, const vector<double>& in_b) const {
    assert(in_a.size() == in_b.size());
    size_t size = in_a.size();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-1, 1);

    double res = 0;
    vector<uint64_t> in_a_ring(size), in_b_ring(size);
    for (size_t i = 0; i < size; i++) {
        res += in_a[i] * in_b[i];
        in_a_ring[i] = static_cast<uint64_t>(in_a[i] * (1ULL << SCALE)) % (1ULL << BIT_LENGTH);
        in_b_ring[i] = static_cast<uint64_t>(in_b[i] * (1ULL << SCALE)) % (1ULL << BIT_LENGTH);
    }
    if (party == ALICE) {
        CipherVector in_a_secret_a(in_a_ring, pk), in_b_secret_a(in_b_ring, pk);
        CipherVector::send(io_pack, &in_a_secret_a);
        CipherVector::send(io_pack, &in_b_secret_a);

        CipherVector res_part1_secret_a, res_part2_secret_a;
        CipherVector::recv(io_pack, &res_part1_secret_a);
        CipherVector::recv(io_pack, &res_part2_secret_a);
        vector<uint64_t> res_part1_ring = res_part1_secret_a.decrypt(sk),
                         res_part2_ring = res_part2_secret_a.decrypt(sk);
        vector<int64_t> res_part1(size), res_part2(size);
        for (size_t i = 0; i < size; i++) {
            res_part1_ring[i] = res_part1_ring[i] % (1ULL << BIT_LENGTH);
            res_part2_ring[i] = res_part2_ring[i] % (1ULL << BIT_LENGTH);
            res_part1[i] =
                ((res_part1_ring[i] >= (1ULL << BIT_LENGTH) / 2) ? (res_part1_ring[i] - (1ULL << BIT_LENGTH)) :
                                                                   res_part1_ring[i]);
            res_part2[i] =
                ((res_part2_ring[i] >= (1ULL << BIT_LENGTH) / 2) ? (res_part2_ring[i] - (1ULL << BIT_LENGTH)) :
                                                                   res_part2_ring[i]);
            res += (res_part1[i] + res_part2[i]) * 1. / (1ULL << (2 * SCALE));
        }
    }
    else {
        CipherVector in_a_secret_a, in_b_secret_a;
        CipherVector::recv(io_pack, &in_a_secret_a);
        CipherVector::recv(io_pack, &in_b_secret_a);

        in_a_secret_a.mul_plain_inplace(in_b_ring);
        in_b_secret_a.mul_plain_inplace(in_a_ring);

        vector<double> res_part1(size), res_part2(size);
        vector<uint64_t> res_part1_ring(size), res_part2_ring(size);
        for (size_t i = 0; i < size; i++) {
            res_part1[i]      = dist(gen);
            res_part2[i]      = dist(gen);
            res_part1_ring[i] = static_cast<uint64_t>(res_part1[i] * (1ULL << (2 * SCALE))) % (1ULL << BIT_LENGTH);
            res_part2_ring[i] = static_cast<uint64_t>(res_part2[i] * (1ULL << (2 * SCALE))) % (1ULL << BIT_LENGTH);
        }

        in_a_secret_a.add_plain_inplace(res_part1_ring);
        in_b_secret_a.add_plain_inplace(res_part2_ring);
        for (size_t i = 0; i < size; i++) {
            res -= (res_part1[i] + res_part2[i]);
        }
        CipherVector::send(io_pack, &in_a_secret_a);
        CipherVector::send(io_pack, &in_b_secret_a);
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
    vector<double> res(data_size), in_b_flatten(data_size * size);
    vector<uint64_t> in_a_flatten_ring(data_size * size), in_b_flatten_ring(data_size * size);
    for (size_t i = 0; i < data_size; i++) {
        for (size_t j = 0; j < size; j++) {
            res[i] += in_a_flatten[i * size + j] * in_b[j];
            in_b_flatten[i * size + j] = in_b[j];
            in_a_flatten_ring[i * size + j] =
                static_cast<uint64_t>(in_a_flatten[i * size + j] * (1ULL << SCALE)) % (1ULL << BIT_LENGTH);
            in_b_flatten_ring[i * size + j] =
                static_cast<uint64_t>(in_b_flatten[i * size + j] * (1ULL << SCALE)) % (1ULL << BIT_LENGTH);
        }
    }

    if (party == ALICE) {
        CipherVector in_a_secret_a(in_a_flatten_ring, pk), in_b_secret_a(in_b_flatten_ring, pk);
        CipherVector::send(io_pack, &in_a_secret_a);
        CipherVector::send(io_pack, &in_b_secret_a);

        CipherVector res_part1_secret_a, res_part2_secret_a;
        CipherVector::recv(io_pack, &res_part1_secret_a);
        CipherVector::recv(io_pack, &res_part2_secret_a);
        vector<uint64_t> res_part1_ring = res_part1_secret_a.decrypt(sk),
                         res_part2_ring = res_part2_secret_a.decrypt(sk);
        vector<int64_t> res_part1(data_size * size), res_part2(data_size * size);
        for (size_t i = 0; i < data_size; i++) {
            for (size_t j = 0; j < size; j++) {
                res_part1_ring[i * size + j] = res_part1_ring[i * size + j] % (1ULL << BIT_LENGTH);
                res_part2_ring[i * size + j] = res_part2_ring[i * size + j] % (1ULL << BIT_LENGTH);
                res_part1[i * size + j]      = (res_part1_ring[i * size + j] >= ((1ULL << BIT_LENGTH) / 2)) ?
                                                   (res_part1_ring[i * size + j] - ((1ULL << BIT_LENGTH))) :
                                                   res_part1_ring[i * size + j];
                res_part2[i * size + j]      = (res_part2_ring[i * size + j] >= ((1ULL << BIT_LENGTH) / 2)) ?
                                                   (res_part2_ring[i * size + j] - ((1ULL << BIT_LENGTH))) :
                                                   res_part2_ring[i * size + j];
                res[i] += (res_part1[i * size + j] + res_part2[i * size + j]) * 1. / (1ULL << (2 * SCALE));
            }
        }
    }
    else {
        CipherVector in_a_secret_a, in_b_secret_a;
        CipherVector::recv(io_pack, &in_a_secret_a);
        CipherVector::recv(io_pack, &in_b_secret_a);

        in_a_secret_a.mul_plain_inplace(in_b_flatten_ring);
        in_b_secret_a.mul_plain_inplace(in_a_flatten_ring);

        vector<double> res_part1(data_size * size), res_part2(data_size * size);
        vector<uint64_t> res_part1_ring(data_size * size), res_part2_ring(data_size * size);
        for (size_t i = 0; i < data_size * size; i++) {
            res_part1[i]      = dist(gen);
            res_part2[i]      = dist(gen);
            res_part1_ring[i] = static_cast<uint64_t>(res_part1[i] * (1ULL << (2 * SCALE))) % (1ULL << BIT_LENGTH);
            res_part2_ring[i] = static_cast<uint64_t>(res_part2[i] * (1ULL << (2 * SCALE))) % (1ULL << BIT_LENGTH);
        }
        in_a_secret_a.add_plain_inplace(res_part1_ring);
        in_b_secret_a.add_plain_inplace(res_part2_ring);
        for (size_t i = 0; i < data_size; i++) {
            for (size_t j = 0; j < size; j++) {
                res[i] -= (res_part1[i * size + j] + res_part2[i * size + j]);
            }
        }
        CipherVector::send(io_pack, &in_a_secret_a);
        CipherVector::send(io_pack, &in_b_secret_a);
    }
    return res;
}
}  // namespace PrivLR_Elgamal