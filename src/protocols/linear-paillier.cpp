#include "linear-paillier.h"

namespace PrivLR_Paillier {
    double Linear::dot_product(const vector<double> &in_a, const vector<double> &in_b) const {
        assert(in_a.size() == in_b.size());
        size_t size = in_a.size();
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(-1, 1);
        std::uniform_real_distribution<> positive_dist(0, 1);

        double res = 0;
        vector<char> in_a_signed(size), in_b_signed(size);
        vector<ZZ> in_a_fixed(size), in_b_fixed(size);
        for (size_t i = 0; i < size; i++) {
            res += in_a[i] * in_b[i];
            in_a_signed[i] = in_a[i] < 0;
            in_b_signed[i] = in_b[i] < 0;
            in_a_fixed[i] = ZZ(uint64_t((in_a_signed[i] ? -in_a[i] : in_a[i]) * (1ull << BIT_LENGTH)));
            in_b_fixed[i] = ZZ(uint64_t((in_b_signed[i] ? -in_b[i] : in_b[i]) * (1ull << BIT_LENGTH)));
        }

        if (party == ALICE) {
            vector<char> in_a_b_signed(size), in_b_b_signed(size);
            CipherVector in_a_secret_a(in_a_fixed, pk), in_b_secret_a(in_b_fixed, pk), res_part1_secret_a, res_part2_secret_a;
            // send in_a_secret_a, in_b_secret_b to bob;
            CipherVector::send(io_pack, &in_a_secret_a);
            CipherVector::send(io_pack, &in_b_secret_a);
            io_pack->send_data(in_a_signed.data(), sizeof(char) * size);
            io_pack->send_data(in_b_signed.data(), sizeof(char) * size);

            // recv res_second from bob;
            CipherVector::recv(io_pack, &res_part1_secret_a);
            CipherVector::recv(io_pack, &res_part2_secret_a);
            io_pack->recv_data(in_a_b_signed.data(), sizeof(char) * size);
            io_pack->recv_data(in_b_b_signed.data(), sizeof(char) * size);
            auto res_part1_fixed = res_part1_secret_a.decrypt(sk);
            auto res_part2_fixed = res_part2_secret_a.decrypt(sk);
            for (size_t i = 0; i < size; i++) {
                const double flag1 = in_a_signed[i] == in_b_b_signed[i] ? 1. : -1.;
                const double flag2 = in_b_signed[i] == in_a_b_signed[i] ? 1. : -1.;
                res += (NTL::to_double(res_part1_fixed[i]) / (1ull << BIT_LENGTH * 2)) * flag1;
                res += (NTL::to_double(res_part2_fixed[i]) / (1ull << BIT_LENGTH * 2)) * flag2;
            }

        } else {
            vector<char> in_a_a_signed(size), in_b_a_signed(size);
            CipherVector in_a_secret_a, in_b_secret_a;
            // recv in_a_secret_a, in_b_secret_b from alice;
            CipherVector::recv(io_pack, &in_a_secret_a);
            CipherVector::recv(io_pack, &in_b_secret_a);
            io_pack->recv_data(in_a_a_signed.data(), sizeof(char) * size);
            io_pack->recv_data(in_b_a_signed.data(), sizeof(char) * size);

            in_a_secret_a.mul_plain_inplace(in_b_fixed);
            in_b_secret_a.mul_plain_inplace(in_a_fixed);
            vector<ZZ> res_second(size);
            for (size_t i = 0; i < size; i++) {
                res_second[i] = ZZ(uint64_t(positive_dist(gen) * (1ull << BIT_LENGTH * 2)));
            }
            in_a_secret_a.add_plain_inplace(res_second);
            in_b_secret_a.add_plain_inplace(res_second);
            for (size_t i = 0; i < size; i++) {
                const double flag1 = in_a_a_signed[i] == in_b_signed[i] ? 1. : -1.;
                const double flag2 = in_b_a_signed[i] == in_a_signed[i] ? 1. : -1.;
                res -= (NTL::to_double(res_second[i]) / (1ull << BIT_LENGTH * 2)) * flag1;
                res -= (NTL::to_double(res_second[i]) / (1ull << BIT_LENGTH * 2)) * flag2;
            }
            // send in_a_secret_a to alice
            CipherVector::send(io_pack, &in_a_secret_a);
            CipherVector::send(io_pack, &in_b_secret_a);
            io_pack->send_data(in_a_signed.data(), sizeof(char) * size);
            io_pack->send_data(in_b_signed.data(), sizeof(char) * size);
        }
        return res;
    }

    vector<double> Linear::dot_product(const vector<vector<double>> &in_a, const vector<double> &in_b, double transpose) const {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(-1, 1);
        std::uniform_real_distribution<> positive_dist(0, 1);
        size_t data_size = transpose ? in_a[0].size() : in_a.size(), size = in_b.size();
        vector<double> in_a_flatten(data_size * size);
        if (transpose) {
            for (size_t i = 0; i < data_size; i++) {
                for (size_t j = 0; j < size; j++) {
                    in_a_flatten[i * size + j] = in_a[j][i];
                }
            }
        } else {
            for (size_t i = 0; i < data_size; i++) {
                for (size_t j = 0; j < size; j++) {
                    in_a_flatten[i * size + j] = in_a[i][j];
                }
            }
        }
        vector<double> res(data_size), in_b_flatten(data_size * size);
        vector<char> in_a_flatten_signed(data_size * size), in_b_signed(size);
        vector<ZZ> in_a_flatten_fixed(data_size * size), in_b_flatten_fixed(data_size * size);
        for (size_t j = 0; j < size; j++) {
            in_b_signed[j] = in_b[j] < 0;
        }
        for (size_t i = 0; i < data_size; i++) {
            for (size_t j = 0; j < size; j++) {
                res[i] += in_a_flatten[i * size + j] * in_b[j];
                in_b_flatten[i * size + j] = in_b[j];
                in_a_flatten_signed[i * size + j] = in_a_flatten[i * size + j] < 0;
                in_a_flatten_fixed[i * size + j] = ZZ(uint64_t((in_a_flatten_signed[i * size + j] ? -in_a_flatten[i * size + j] : in_a_flatten[i * size + j]) * (1ull << BIT_LENGTH)));
                in_b_flatten_fixed[i * size + j] = ZZ(uint64_t((in_b_signed[j] ? -in_b_flatten[i * size + j] : in_b_flatten[i * size + j]) * (1ull << BIT_LENGTH)));
            }
        }

        if (party == ALICE) {
            vector<char> in_a_b_flatten_signed(data_size * size), in_b_b_signed(size);
            CipherVector in_a_flatten_secret_a(in_a_flatten_fixed, pk), in_b_flatten_secret_a(in_b_flatten_fixed, pk), res_part1_secret_a, res_part2_secret_a;
            // send in_a_secret_a, in_b_secret_b to bob;
            CipherVector::send(io_pack, &in_a_flatten_secret_a);
            CipherVector::send(io_pack, &in_b_flatten_secret_a);
            io_pack->send_data(in_a_flatten_signed.data(), sizeof(char) * data_size * size);
            io_pack->send_data(in_b_signed.data(), sizeof(char) * size);

            // recv res_second from bob;
            CipherVector::recv(io_pack, &res_part1_secret_a);
            CipherVector::recv(io_pack, &res_part2_secret_a);
            io_pack->recv_data(in_a_b_flatten_signed.data(), sizeof(char) * data_size * size);
            io_pack->recv_data(in_b_b_signed.data(), sizeof(char) * size);
            auto res_part1_fixed = res_part1_secret_a.decrypt(sk);
            auto res_part2_fixed = res_part2_secret_a.decrypt(sk);
            for (size_t i = 0; i < data_size; i++) {
                for (size_t j = 0; j < size; j++) {
                    const double flag1 = in_a_flatten_signed[i * size + j] == in_b_b_signed[j] ? 1. : -1.;
                    const double flag2 = in_b_signed[j] == in_a_b_flatten_signed[i * size + j] ? 1. : -1.;
                    res[i] += (NTL::to_double(res_part1_fixed[i * size + j]) / (1ull << BIT_LENGTH * 2)) * flag1;
                    res[i] += (NTL::to_double(res_part2_fixed[i * size + j]) / (1ull << BIT_LENGTH * 2)) * flag2;
                }
            }
        } else {
            vector<char> in_a_a_flatten_signed(data_size * size), in_b_a_signed(size);
            CipherVector in_a_flatten_secret_a, in_b_flatten_secret_a;
            // recv in_a_secret_a, in_b_secret_b from alice;
            CipherVector::recv(io_pack, &in_a_flatten_secret_a);
            CipherVector::recv(io_pack, &in_b_flatten_secret_a);
            io_pack->recv_data(in_a_a_flatten_signed.data(), sizeof(char) * data_size * size);
            io_pack->recv_data(in_b_a_signed.data(), sizeof(char) * size);

            in_a_flatten_secret_a.mul_plain_inplace(in_b_flatten_fixed);
            in_b_flatten_secret_a.mul_plain_inplace(in_a_flatten_fixed);
            vector<ZZ> res_second(data_size * size);
            for (size_t i = 0; i < data_size * size; i++) {
                res_second[i] = ZZ(uint64_t(positive_dist(gen) * (1ull << BIT_LENGTH * 2)));
            }
            in_a_flatten_secret_a.add_plain_inplace(res_second);
            in_b_flatten_secret_a.add_plain_inplace(res_second);
            for (size_t i = 0; i < data_size; i++) {
                for (size_t j = 0; j < size; j++) {
                    const double flag1 = in_a_a_flatten_signed[i * size + j] == in_b_signed[j] ? 1. : -1.;
                    const double flag2 = in_b_a_signed[j] == in_a_flatten_signed[i * size + j] ? 1. : -1.;
                    res[i] -= (NTL::to_double(res_second[i * size + j]) / (1ull << BIT_LENGTH * 2)) * flag1;
                    res[i] -= (NTL::to_double(res_second[i * size + j]) / (1ull << BIT_LENGTH * 2)) * flag2;
                }
            }
            // send in_a_secret_a to alice
            CipherVector::send(io_pack, &in_a_flatten_secret_a);
            CipherVector::send(io_pack, &in_b_flatten_secret_a);
            io_pack->send_data(in_a_flatten_signed.data(), sizeof(char) * data_size * size);
            io_pack->send_data(in_b_signed.data(), sizeof(char) * size);
        }
        return res;
    }
}