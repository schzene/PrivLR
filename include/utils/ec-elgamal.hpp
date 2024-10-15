#ifndef PRIV_LR_EC_ELGAMAL_HPP__
#define PRIV_LR_EC_ELGAMAL_HPP__
#include <cstddef>
#include <openssl/ec.h>
#include <vector>

#include "io.h"
extern "C" {
#include "ec-elgamal.h"
}

using std::vector;

namespace EC_Elgamal {
class CipherVector {
public:
    class lenth_error : public std::exception {
        const char *message;

    public:
        lenth_error(const char *msg) : message(msg) {}
        const char *what() const throw() override {
            return message;
        }
    };

    vector<ec_elgamal_ciphertext*> data;

    CipherVector() {}

    CipherVector(const vector<uint64_t>& messages, const ec_elgamal_public_key& public_key) {
        const EC_GROUP* init_group = ec_elgamal_get_group();
        BIGNUM *bn_plain = BN_new(), *ord = BN_new(), *rand = BN_new();
        BN_CTX* ctx = BN_CTX_new();

        EC_GROUP_get_order(init_group, ord, ctx);
        BN_rand_range(rand, ord);

        size_t size = messages.size();
        data.resize(size);
        for (size_t i = 0; i < size; i++) {
            data[i] = ec_elgamal_new_ciphertext();
            BN_set_word(bn_plain, messages[i]);

            EC_POINT_mul(init_group, data[i]->C1, NULL, EC_GROUP_get0_generator(init_group), rand, ctx);
            EC_POINT_mul(init_group, data[i]->C2, bn_plain, public_key, rand, ctx);
        }

        BN_clear_free(rand);
        BN_free(ord);
        BN_free(bn_plain);
        BN_CTX_free(ctx);
    }

    ~CipherVector() {
        size_t size = data.size();
        for (size_t i = 0; i < size; i++) {
            ec_elgamal_free_ciphertext(data[i]);
        }
    }

    vector<uint64_t> decrypt(const ec_elgamal_secret_key& secret_key, const bsgs_table_t table = NULL) const {
        const EC_GROUP* init_group = ec_elgamal_get_group();
        BN_CTX* ctx                = BN_CTX_new();
        EC_POINT* M                = EC_POINT_new(init_group);
        size_t size                = data.size();
        vector<uint64_t> res(size);

        for (size_t i = 0; i < size; i++) {
            EC_POINT_mul(init_group, M, NULL, data[i]->C1, secret_key, ctx);
            EC_POINT_invert(init_group, M, ctx);
            EC_POINT_add(init_group, M, data[i]->C2, M, ctx);

            if (table != NULL) {
                ecdlp_bsgs(table, M, &res[i], 1L << MAX_BITS);
            }
            else {
                ecdlp_brute(init_group, M, &res[i], 1L << MAX_BITS);
            }
        }

        BN_CTX_free(ctx);
        EC_POINT_clear_free(M);

        return res;
    }

    CipherVector add(const CipherVector& other) const {
        const EC_GROUP* init_group = ec_elgamal_get_group();
        BN_CTX *ctx = BN_CTX_new();
        CipherVector ret;
        const size_t this_size = data.size();
        const size_t other_size = other.data.size();
        if (this_size == 1) {
            ret.data.resize(other_size); 
            for (size_t i = 0; i < other_size; i++) {
                ret.data[i] = ec_elgamal_new_ciphertext();
                EC_POINT_add(init_group, ret.data[i]->C1, data[0]->C1, other.data[i]->C1, ctx);
                EC_POINT_add(init_group, ret.data[i]->C2, data[0]->C2, other.data[i]->C2, ctx);
            }
        } else if (other_size == 1) {
            ret.data.resize(this_size); 
            for (size_t i = 0; i < this_size; i++) {
                ret.data[i] = ec_elgamal_new_ciphertext();
                EC_POINT_add(init_group, ret.data[i]->C1, data[i]->C1, other.data[0]->C1, ctx);
                EC_POINT_add(init_group, ret.data[i]->C2, data[i]->C2, other.data[0]->C2, ctx);
            }
        } else if (this_size == other_size) {
            ret.data.resize(this_size); 
            for (size_t i = 0; i < this_size; i++) {
                ret.data[i] = ec_elgamal_new_ciphertext();
                EC_POINT_add(init_group, ret.data[i]->C1, data[i]->C1, other.data[i]->C1, ctx);
                EC_POINT_add(init_group, ret.data[i]->C2, data[i]->C2, other.data[i]->C2, ctx);
            }
        } else {
            char buf[100];
            sprintf(buf, "Length of CipherVector(%ld) and CipherVector(%ld) mismatch", this_size, other_size);
            throw lenth_error(buf);
        }
        BN_CTX_free(ctx);

        return ret;
    }

    void add_inplace(const CipherVector& other) {
        const EC_GROUP* init_group = ec_elgamal_get_group();
        BN_CTX *ctx = BN_CTX_new();
        const size_t this_size = data.size();
        const size_t other_size = other.data.size();
        if (this_size == 1) {
            EC_POINT* C1 = EC_POINT_new(init_group);
            EC_POINT_copy(C1, data[0]->C1);
            EC_POINT* C2 = EC_POINT_new(init_group);
            EC_POINT_copy(C2, data[0]->C2);
            data.resize(other_size);
            EC_POINT_add(init_group, data[0]->C1, C1, other.data[0]->C1, ctx);
            EC_POINT_add(init_group, data[0]->C2, C2, other.data[0]->C2, ctx);
            for (size_t i = 1; i < other_size; i++) {
                data[i] = ec_elgamal_new_ciphertext();
                EC_POINT_add(init_group, data[i]->C1, C1, other.data[i]->C1, ctx);
                EC_POINT_add(init_group, data[i]->C2, C2, other.data[i]->C2, ctx);
            }
            EC_POINT_free(C1);
            EC_POINT_free(C2);
        } else if (other_size == 1) {
            for (size_t i = 0; i < this_size; i++) {
                EC_POINT_add(init_group, data[i]->C1, data[i]->C1, other.data[0]->C1, ctx);
                EC_POINT_add(init_group, data[i]->C2, data[i]->C2, other.data[0]->C2, ctx);
            }
        } else if (this_size == other_size) {
            for (size_t i = 0; i < this_size; i++) {
                EC_POINT_add(init_group, data[i]->C1, data[i]->C1, other.data[i]->C1, ctx);
                EC_POINT_add(init_group, data[i]->C2, data[i]->C2, other.data[i]->C2, ctx);
            }
        } else {
            char buf[100];
            sprintf(buf, "Length of CipherVector(%ld) and CipherVector(%ld) mismatch", this_size, other_size);
            throw lenth_error(buf);
        }
        BN_CTX_free(ctx);
    }

    CipherVector add_plain(const vector<uint64_t>& other) const {
        const EC_GROUP* init_group = ec_elgamal_get_group();
        BN_CTX *ctx = BN_CTX_new();
        CipherVector ret;
        const size_t this_size = data.size();
        const size_t other_size = other.size();
        BIGNUM *bn_plain = BN_new();
        EC_POINT *mG = EC_POINT_new(init_group);
        if (this_size == 1) {
            ret.data.resize(other_size); 
            for (size_t i = 0; i < other_size; i++) {
                ret.data[i] = ec_elgamal_new_ciphertext();
                BN_set_word(bn_plain, other[i]);
                EC_POINT_mul(init_group, mG, NULL, EC_GROUP_get0_generator(init_group), bn_plain, ctx);
                EC_POINT_copy(ret.data[i]->C1, data[0]->C1);
                EC_POINT_add(init_group, ret.data[i]->C2, data[0]->C2, mG, ctx);
            }
        } else if (other_size == 1) {
            ret.data.resize(this_size); 
            for (size_t i = 0; i < this_size; i++) {
                ret.data[i] = ec_elgamal_new_ciphertext();
                BN_set_word(bn_plain, other[0]);
                EC_POINT_mul(init_group, mG, NULL, EC_GROUP_get0_generator(init_group), bn_plain, ctx);
                EC_POINT_copy(ret.data[i]->C1, data[i]->C1);
                EC_POINT_add(init_group, ret.data[i]->C2, data[i]->C2, mG, ctx);
            }
        } else if (this_size == other_size) {
            ret.data.resize(this_size); 
            for (size_t i = 0; i < this_size; i++) {
                ret.data[i] = ec_elgamal_new_ciphertext();
                BN_set_word(bn_plain, other[i]);
                EC_POINT_mul(init_group, mG, NULL, EC_GROUP_get0_generator(init_group), bn_plain, ctx);
                EC_POINT_copy(ret.data[i]->C1, data[i]->C1);
                EC_POINT_add(init_group, ret.data[i]->C2, data[i]->C2, mG, ctx);
            }
        } else {
            char buf[100];
            sprintf(buf, "Length of CipherVector(%ld) and Plaintext(%ld) mismatch", this_size, other_size);
            throw lenth_error(buf);
        }
        EC_POINT_free(mG);
        BN_free(bn_plain);
        BN_CTX_free(ctx);

        return ret;
    }

    void add_plain_inplace(const vector<uint64_t>& other) {
        const EC_GROUP* init_group = ec_elgamal_get_group();
        BN_CTX *ctx = BN_CTX_new();
        const size_t this_size = data.size();
        const size_t other_size = other.size();
        BIGNUM *bn_plain = BN_new();
        EC_POINT *mG = EC_POINT_new(init_group);
        if (this_size == 1) {
            EC_POINT* C1 = EC_POINT_new(init_group);
            EC_POINT_copy(C1, data[0]->C1);
            EC_POINT* C2 = EC_POINT_new(init_group);
            EC_POINT_copy(C2, data[0]->C2);
            data.resize(other_size);

            BN_set_word(bn_plain, other[0]);
            EC_POINT_mul(init_group, mG, NULL, EC_GROUP_get0_generator(init_group), bn_plain, ctx);
            EC_POINT_copy(data[0]->C1, C1);
            EC_POINT_add(init_group, data[0]->C2, C2, mG, ctx);
            for (size_t i = 1; i < other_size; i++) {
                data[i] = ec_elgamal_new_ciphertext();
                EC_POINT_mul(init_group, mG, NULL, EC_GROUP_get0_generator(init_group), bn_plain, ctx);
                EC_POINT_copy(data[i]->C1, C1);
                EC_POINT_add(init_group, data[i]->C2, C2, mG, ctx);
            }
            EC_POINT_free(C1);
            EC_POINT_free(C2);
        } else if (other_size == 1) {
            for (size_t i = 0; i < this_size; i++) {
                BN_set_word(bn_plain, other[0]);
                EC_POINT_mul(init_group, mG, NULL, EC_GROUP_get0_generator(init_group), bn_plain, ctx);
                EC_POINT_add(init_group, data[i]->C2, data[i]->C2, mG, ctx);
            }
        } else if (this_size == other_size) {
            for (size_t i = 0; i < this_size; i++) {
                BN_set_word(bn_plain, other[i]);
                EC_POINT_mul(init_group, mG, NULL, EC_GROUP_get0_generator(init_group), bn_plain, ctx);
                EC_POINT_add(init_group, data[i]->C2, data[i]->C2, mG, ctx);
            }
        } else {
            char buf[100];
            sprintf(buf, "Length of CipherVector(%ld) and Plaintext(%ld) mismatch", this_size, other_size);
            throw lenth_error(buf);
        }
        EC_POINT_free(mG);
        BN_free(bn_plain);
        BN_CTX_free(ctx);
    }

    CipherVector mul_plain(const vector<uint64_t>& other) const {
        const EC_GROUP* init_group = ec_elgamal_get_group();
        BN_CTX *ctx = BN_CTX_new();
        CipherVector ret;
        const size_t this_size = data.size();
        const size_t other_size = other.size();
        BIGNUM *bn_plain = BN_new();

        if (this_size == 1) {
            ret.data.resize(other_size); 
            for (size_t i = 0; i < other_size; i++) {
                ret.data[i] = ec_elgamal_new_ciphertext();
                BN_set_word(bn_plain, other[i]);
                EC_POINT_mul(init_group, ret.data[i]->C1, NULL, data[0]->C1, bn_plain, ctx);
                EC_POINT_mul(init_group, ret.data[i]->C2, NULL, data[0]->C2, bn_plain, ctx);
            }
        } else if (other_size == 1) {
            ret.data.resize(this_size); 
            for (size_t i = 0; i < this_size; i++) {
                ret.data[i] = ec_elgamal_new_ciphertext();
                BN_set_word(bn_plain, other[0]);
                EC_POINT_mul(init_group, ret.data[i]->C1, NULL, data[i]->C1, bn_plain, ctx);
                EC_POINT_mul(init_group, ret.data[i]->C2, NULL, data[i]->C2, bn_plain, ctx);
            }
        } else if (this_size == other_size) {
            ret.data.resize(this_size); 
            for (size_t i = 0; i < this_size; i++) {
                ret.data[i] = ec_elgamal_new_ciphertext();
                BN_set_word(bn_plain, other[i]);
                EC_POINT_mul(init_group, ret.data[i]->C1, NULL, data[i]->C1, bn_plain, ctx);
                EC_POINT_mul(init_group, ret.data[i]->C2, NULL, data[i]->C2, bn_plain, ctx);
            }
        } else {
            char buf[100];
            sprintf(buf, "Length of CipherVector(%ld) and Plaintext(%ld) mismatch", this_size, other_size);
            throw lenth_error(buf);
        }
        BN_free(bn_plain);
        BN_CTX_free(ctx);

        return ret;
    }

    void mul_plain_inplace(const vector<uint64_t>& other) {
        const EC_GROUP* init_group = ec_elgamal_get_group();
        BN_CTX *ctx = BN_CTX_new();
        const size_t this_size = data.size();
        const size_t other_size = other.size();
        BIGNUM *bn_plain = BN_new();

        if (this_size == 1) {
            EC_POINT* C1 = EC_POINT_new(init_group);
            EC_POINT_copy(C1, data[0]->C1);
            EC_POINT* C2 = EC_POINT_new(init_group);
            EC_POINT_copy(C2, data[0]->C2);
            data.resize(other_size);

            BN_set_word(bn_plain, other[0]);
            EC_POINT_mul(init_group, data[0]->C1, NULL, C1, bn_plain, ctx);
            EC_POINT_mul(init_group, data[0]->C2, NULL, C2, bn_plain, ctx);
            for (size_t i = 1; i < other_size; i++) {
                data[i] = ec_elgamal_new_ciphertext();
                BN_set_word(bn_plain, other[i]);
                EC_POINT_mul(init_group, data[i]->C1, NULL, C1, bn_plain, ctx);
                EC_POINT_mul(init_group, data[i]->C2, NULL, C2, bn_plain, ctx);
            }
            EC_POINT_free(C1);
            EC_POINT_free(C2);
        } else if (other_size == 1) {
            for (size_t i = 0; i < this_size; i++) {
                BN_set_word(bn_plain, other[0]);
                EC_POINT_mul(init_group, data[i]->C1, NULL, data[i]->C1, bn_plain, ctx);
                EC_POINT_mul(init_group, data[i]->C2, NULL, data[i]->C2, bn_plain, ctx);
            }
        } else if (this_size == other_size) {
            for (size_t i = 0; i < this_size; i++) {
                BN_set_word(bn_plain, other[i]);
                EC_POINT_mul(init_group, data[i]->C1, NULL, data[i]->C1, bn_plain, ctx);
                EC_POINT_mul(init_group, data[i]->C2, NULL, data[i]->C2, bn_plain, ctx);
            }
        } else {
            char buf[100];
            sprintf(buf, "Length of CipherVector(%ld) and Plaintext(%ld) mismatch", this_size, other_size);
            throw lenth_error(buf);
        }
        BN_free(bn_plain);
        BN_CTX_free(ctx);
    }

    static void send(IOPack* io_pack, CipherVector* cv);
    static void recv(IOPack* io_pack, CipherVector* cv);
};
}  // namespace EC_Elgamal
#endif