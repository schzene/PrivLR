#include "linear-bfv.h"
#include "utils/bfv-tools.h"

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

namespace PrivLR_BFV {
double Linear::dot_product(const vector<double>& in_a, const vector<double>& in_b) const {
    assert(in_a.size() == in_b.size());
    size_t size = in_a.size();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-1, 1);

    double res = 0;
    vector<uint64_t> in_a_ring(size), in_b_ring(size), in_a_prime(size), in_b_prime(size);
    for (size_t i = 0; i < size; i++) {
        res += in_a[i] * in_b[i];
        in_a_ring[i]  = static_cast<uint64_t>(in_a[i] * (1ULL << SCALE)) % (1ULL << BIT_LENGTH);
        in_b_ring[i]  = static_cast<uint64_t>(in_b[i] * (1ULL << SCALE)) % (1ULL << BIT_LENGTH);
        in_a_prime[i] = static_cast<uint64_t>(in_a[i] * (1ULL << SCALE)) % party->parm->plain_mod;
        in_b_prime[i] = static_cast<uint64_t>(in_b[i] * (1ULL << SCALE)) % party->parm->plain_mod;
    }
    BFVLongPlaintext in_a_plain(party->parm, in_a_prime), in_b_plain(party->parm, in_b_prime);
    if (*party == ALICE) {
        BFVLongCiphertext lct_a(in_a_plain, party), lct_b(in_b_plain, party);
        BFVLongCiphertext::send(io_pack->io, &lct_a);
        BFVLongCiphertext::send(io_pack->io, &lct_b);

        BFVLongCiphertext res_part1_secret_a, res_part2_secret_a;
        BFVLongCiphertext::recv(io_pack->io_rev, &res_part1_secret_a, party->parm->context);
        BFVLongCiphertext::recv(io_pack->io_rev, &res_part2_secret_a, party->parm->context);
        BFVLongPlaintext res_part1_plain = res_part1_secret_a.decrypt(party),
                         res_part2_plain = res_part2_secret_a.decrypt(party);
        vector<uint64_t> res_part1_prime = res_part1_plain.decode_uint(party->parm),
                         res_part2_prime = res_part2_plain.decode_uint(party->parm);
        vector<int64_t> res_part1(size), res_part2(size);
        for (size_t i = 0; i < size; i++) {
            res_part1[i] =
                ((res_part1_prime[i] >= party->parm->plain_mod / 2) ? (res_part1_prime[i] - party->parm->plain_mod) :
                                                                      res_part1_prime[i]);
            res_part2[i] =
                ((res_part2_prime[i] >= party->parm->plain_mod / 2) ? (res_part2_prime[i] - party->parm->plain_mod) :
                                                                      res_part2_prime[i]);
            res += (res_part1[i] + res_part2[i]) * 1. / (1ULL << (2 * SCALE));
        }
    }
    else {
        BFVLongCiphertext lct_a, lct_b;
        BFVLongCiphertext::recv(io_pack->io_rev, &lct_a, party->parm->context);
        BFVLongCiphertext::recv(io_pack->io_rev, &lct_b, party->parm->context);

        lct_a.multiply_plain_inplace(in_b_plain, party->parm->evaluator);
        lct_b.multiply_plain_inplace(in_a_plain, party->parm->evaluator);

        vector<double> res_part1(size), res_part2(size);
        vector<uint64_t> res_part1_prime(size), res_part2_prime(size);
        for (size_t i = 0; i < size; i++) {
            res_part1[i]       = dist(gen);
            res_part2[i]       = dist(gen);
            res_part1_prime[i] = static_cast<uint64_t>(res_part1[i] * (1ULL << (2 * SCALE))) % party->parm->plain_mod;
            res_part2_prime[i] = static_cast<uint64_t>(res_part2[i] * (1ULL << (2 * SCALE))) % party->parm->plain_mod;
        }
        BFVLongPlaintext res_part1_plain(party->parm, res_part1_prime), res_part2_plain(party->parm, res_part2_prime);
        lct_a.add_plain_inplace(res_part1_plain, party->parm->evaluator);
        lct_b.add_plain_inplace(res_part2_plain, party->parm->evaluator);
        for (size_t i = 0; i < size; i++) {
            res -= (res_part1[i] + res_part2[i]);
        }
        BFVLongCiphertext::send(io_pack->io, &lct_a);
        BFVLongCiphertext::send(io_pack->io, &lct_b);
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
    vector<uint64_t> in_a_flatten_ring(data_size * size), in_b_flatten_ring(data_size * size),
        in_a_flatten_prime(data_size * size), in_b_flatten_prime(data_size * size);
    for (size_t i = 0; i < data_size; i++) {
        for (size_t j = 0; j < size; j++) {
            res[i] += in_a_flatten[i * size + j] * in_b[j];
            in_b_flatten[i * size + j] = in_b[j];
            in_a_flatten_ring[i * size + j] =
                static_cast<uint64_t>(in_a_flatten[i * size + j] * (1ULL << SCALE)) % (1ULL << BIT_LENGTH);
            in_b_flatten_ring[i * size + j] =
                static_cast<uint64_t>(in_b_flatten[i * size + j] * (1ULL << SCALE)) % (1ULL << BIT_LENGTH);
            in_a_flatten_prime[i * size + j] =
                static_cast<uint64_t>(in_a_flatten[i * size + j] * (1ULL << SCALE)) % party->parm->plain_mod;
            in_b_flatten_prime[i * size + j] =
                static_cast<uint64_t>(in_b_flatten[i * size + j] * (1ULL << SCALE)) % party->parm->plain_mod;
        }
    }

    BFVLongPlaintext in_a_flatten_plain(party->parm, in_a_flatten_prime),
        in_b_flatten_plain(party->parm, in_b_flatten_prime);
    if (*party == ALICE) {
        BFVLongCiphertext lct_a(in_a_flatten_plain, party), lct_b(in_b_flatten_plain, party);
        BFVLongCiphertext::send(io_pack->io, &lct_a);
        BFVLongCiphertext::send(io_pack->io, &lct_b);

        BFVLongCiphertext res_part1_secret_a, res_part2_secret_a;
        BFVLongCiphertext::recv(io_pack->io_rev, &res_part1_secret_a, party->parm->context);
        BFVLongCiphertext::recv(io_pack->io_rev, &res_part2_secret_a, party->parm->context);
        BFVLongPlaintext res_part1_plain = res_part1_secret_a.decrypt(party),
                         res_part2_plain = res_part2_secret_a.decrypt(party);
        vector<uint64_t> res_part1_prime = res_part1_plain.decode_uint(party->parm),
                         res_part2_prime = res_part2_plain.decode_uint(party->parm);
        vector<int64_t> res_part1(data_size * size), res_part2(data_size * size);
        for (size_t i = 0; i < data_size; i++) {
            for (size_t j = 0; j < size; j++) {
                res_part1[i * size + j] = ((res_part1_prime[i * size + j] >= party->parm->plain_mod / 2) ?
                                               (res_part1_prime[i * size + j] - party->parm->plain_mod) :
                                               res_part1_prime[i * size + j]);
                res_part2[i * size + j] = ((res_part2_prime[i * size + j] >= party->parm->plain_mod / 2) ?
                                               (res_part2_prime[i * size + j] - party->parm->plain_mod) :
                                               res_part2_prime[i * size + j]);
                res[i] += (res_part1[i * size + j] + res_part2[i * size + j]) * 1. / (1ULL << (2 * SCALE));
            }
        }
    }
    else {
        BFVLongCiphertext lct_a, lct_b;
        BFVLongCiphertext::recv(io_pack->io_rev, &lct_a, party->parm->context);
        BFVLongCiphertext::recv(io_pack->io_rev, &lct_b, party->parm->context);

        lct_a.multiply_plain_inplace(in_b_flatten_plain, party->parm->evaluator);
        lct_b.multiply_plain_inplace(in_a_flatten_plain, party->parm->evaluator);

        vector<double> res_part1(data_size * size), res_part2(data_size * size);
        vector<uint64_t> res_part1_prime(data_size * size), res_part2_prime(data_size * size);
        for (size_t i = 0; i < data_size * size; i++) {
            res_part1[i]       = dist(gen);
            res_part2[i]       = dist(gen);
            res_part1_prime[i] = static_cast<uint64_t>(res_part1[i] * (1ULL << (2 * SCALE))) % party->parm->plain_mod;
            res_part2_prime[i] = static_cast<uint64_t>(res_part2[i] * (1ULL << (2 * SCALE))) % party->parm->plain_mod;
        }
        BFVLongPlaintext res_part1_plain(party->parm, res_part1_prime), res_part2_plain(party->parm, res_part2_prime);
        lct_a.add_plain_inplace(res_part1_plain, party->parm->evaluator);
        lct_b.add_plain_inplace(res_part2_plain, party->parm->evaluator);
        for (size_t i = 0; i < data_size; i++) {
            for (size_t j = 0; j < size; j++) {
                res[i] -= (res_part1[i * size + j] + res_part2[i * size + j]);
            }
        }
        BFVLongCiphertext::send(io_pack->io, &lct_a);
        BFVLongCiphertext::send(io_pack->io, &lct_b);
    }
    return res;
}
}  // namespace PrivLR_BFV