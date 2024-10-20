/*
 * Copyright (c) 2018, Institute for Pervasive Computing, ETH Zurich.
 * All rights reserved.
 *
 * Author:
 *       Lukas Burkhalter <lubu@inf.ethz.ch>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE
 * COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "ec-elgamal.h"
#include <openssl/bn.h>
#include <openssl/ec.h>

char* point_to_string(EC_GROUP* curve_group, const EC_POINT* point) {
    BN_CTX* ctx = BN_CTX_new();
    char* s;
    point_conversion_form_t form = POINT_CONVERSION_COMPRESSED;
    s                            = EC_POINT_point2hex(curve_group, point, form, ctx);
    BN_CTX_free(ctx);
    return s;
}

int add_value_to_table(bsgs_table_t table, EC_POINT* point, uint32_t value) {
    unsigned char* point_key;
    BN_CTX* ctx       = BN_CTX_new();
    size_t point_size = EC_POINT_point2oct(table->group, point, POINT_CONVERSION_COMPRESSED, NULL, 0, ctx);
    bsgs_hash_table_entry_t* new_entry = (bsgs_hash_table_entry_t*)malloc(sizeof(bsgs_hash_table_entry_t));
    point_key                          = (unsigned char*)malloc(point_size);
    EC_POINT_point2oct(table->group, point, POINT_CONVERSION_COMPRESSED, point_key, point_size, ctx);

    new_entry->key   = point_key;
    new_entry->value = value;
    HASH_ADD_KEYPTR(hh, table->table, point_key, point_size, new_entry);
    BN_CTX_free(ctx);
    return 0;
}

int get_power_from_table(uint64_t* power, const bsgs_table_t bsgs_table, const EC_POINT* lookup_point) {
    unsigned char* point_key;
    BN_CTX* ctx       = BN_CTX_new();
    size_t point_size = EC_POINT_point2oct(bsgs_table->group, lookup_point, POINT_CONVERSION_COMPRESSED, NULL, 0, ctx);
    bsgs_hash_table_entry_t* entry;
    point_key = (unsigned char*)malloc(point_size);
    EC_POINT_point2oct(bsgs_table->group, lookup_point, POINT_CONVERSION_COMPRESSED, point_key, point_size, ctx);

    HASH_FIND(hh, bsgs_table->table, point_key, point_size, entry);
    BN_CTX_free(ctx);
    free(point_key);

    if (entry == NULL)
        return -1;
    *power = (uint64_t)entry->value;
    return 0;
}

int ecdlp_bsgs(const bsgs_table_t bsgs_table, const EC_POINT* M, uint64_t* x, int64_t max_it) {
    uint64_t j = 0, i = 0;
    int ok;
    EC_GROUP* curve_group = bsgs_table->group;
    EC_POINT* curPoint    = EC_POINT_dup(M, curve_group);
    EC_POINT* curPointNeg = EC_POINT_dup(M, curve_group);
    BN_CTX* ctx           = BN_CTX_new();

    while (i <= max_it) {
        ok = get_power_from_table(&j, bsgs_table, curPoint);
        if (ok == 0) {
            *x = i * bsgs_table->tablesize + j;
            break;
        }
        EC_POINT_add(curve_group, curPoint, curPoint, bsgs_table->mG_inv, ctx);
        i = i + 1;
    }

    if (i > max_it) {
        return -1;
    }

    EC_POINT_free(curPoint);
    EC_POINT_free(curPointNeg);
    BN_CTX_free(ctx);
    return 0;
}

// Finds the value x with brute force s.t. M=xG
int ecdlp_brute(const EC_GROUP* curve_group, const EC_POINT* M, uint64_t* x, uint64_t max_it) {
    EC_POINT* cur;
    const EC_POINT* G;
    uint64_t max, x_local = 1;
    BN_CTX* ctx = BN_CTX_new();
    max         = (int64_t)max_it;

    cur = EC_POINT_new(curve_group);
    G   = EC_GROUP_get0_generator(curve_group);
    EC_POINT_set_to_infinity(curve_group, cur);

    if (EC_POINT_is_at_infinity(curve_group, M)) {
        *x = 0;
        return 0;
    }
    else {
        for (; x_local < max; (*x) = x_local++) {
            EC_POINT_add(curve_group, cur, cur, G, ctx);
            if (EC_POINT_cmp(curve_group, cur, M, ctx) == 0) {
                break;
            }
        }
        *x = x_local;
    }
    EC_POINT_free(cur);
    BN_CTX_free(ctx);
    return 0;
}

// API IMPLEMENTATION

//the ec group used
EC_GROUP* init_group = NULL;

int ec_elgamal_init(int curve_id) {
    init_group = EC_GROUP_new_by_curve_name(curve_id);
    return 0;
}

const EC_GROUP* ec_elgamal_get_group() {
    return init_group;
}

int ec_elgamal_deinit() {
    if (init_group != NULL) {
        EC_GROUP_free(init_group);
        init_group = NULL;
    }
    return 0;
}

int bsgs_table_init(bsgs_table_t table, uint64_t t_size) {
    uint64_t count      = 0;
    BIGNUM* bn_size     = BN_new();
    EC_POINT* cur_point = EC_POINT_new(init_group);
    const EC_POINT* gen = EC_GROUP_get0_generator(init_group);
    BN_CTX* ctx         = BN_CTX_new();

    table->table  = NULL;
    table->group  = init_group;
    table->mG     = EC_POINT_new(init_group);
    table->mG_inv = EC_POINT_new(init_group);

    //set Table metadata
    BN_set_word(bn_size, (BN_ULONG)t_size);
    table->tablesize = t_size;
    EC_POINT_mul(init_group, table->mG, NULL, gen, bn_size, ctx);
    BN_set_negative(bn_size, 1);
    EC_POINT_mul(init_group, table->mG_inv, NULL, gen, bn_size, ctx);
    BN_free(bn_size);

    EC_POINT_set_to_infinity(init_group, cur_point);
    for (; count <= t_size; count++) {
        add_value_to_table(table, cur_point, count);
        EC_POINT_add(init_group, cur_point, cur_point, gen, ctx);
    }

    EC_POINT_free(cur_point);
    BN_CTX_free(ctx);
    return 0;
}

int bsgs_table_free(bsgs_table_t bsgs_table) {
    bsgs_hash_table_entry_t *tmp, *current;
    HASH_ITER(hh, bsgs_table->table, current, tmp) {
        HASH_DEL(bsgs_table->table, current);
        free(current->key);
        free(current);
    }
    EC_POINT_free(bsgs_table->mG);
    EC_POINT_free(bsgs_table->mG_inv);
    return 0;
}

ec_elgamal_secret_key ec_elgamal_keygen() {
    ec_elgamal_secret_key sk = BN_new();
    BN_CTX* ctx              = BN_CTX_new();
    BIGNUM* ord              = BN_new();
    EC_GROUP_get_order(init_group, ord, ctx);
    BN_rand_range(sk, ord);
    BN_free(ord);
    BN_CTX_free(ctx);
    return sk;
}

ec_elgamal_public_key ec_elgamal_from_secret_key(const ec_elgamal_secret_key sk) {
    ec_elgamal_public_key pk = EC_POINT_new(init_group);
    BN_CTX* ctx              = BN_CTX_new();
    EC_POINT_mul(init_group, pk, NULL, EC_GROUP_get0_generator(init_group), sk, ctx);
    BN_CTX_free(ctx);
    return pk;
}

int ec_elgamal_free_public_key(ec_elgamal_public_key key) {
    if (key != NULL) {
        EC_POINT_clear_free(key);
    }
    return 0;
}

int ec_elgamal_free_secret_key(ec_elgamal_secret_key key) {
    if (key != NULL) {
        BN_clear_free(key);
    }
    return 0;
}

ec_elgamal_ciphertext* ec_elgamal_new_ciphertext() {
    ec_elgamal_ciphertext* ciphertext = (ec_elgamal_ciphertext*)malloc(sizeof(ec_elgamal_ciphertext));
    ciphertext->C1                    = EC_POINT_new(init_group);
    ciphertext->C2                    = EC_POINT_new(init_group);
    return ciphertext;
}

int ec_elgamal_free_ciphertext(ec_elgamal_ciphertext* ciphertext) {
    if (ciphertext != NULL) {
        EC_POINT_clear_free(ciphertext->C1);
        EC_POINT_clear_free(ciphertext->C2);
        free(ciphertext);
        ciphertext = NULL;
    }
    return 0;
}

int ec_elgamal_encrypt(ec_elgamal_ciphertext* ciphertext, const ec_elgamal_public_key pk, const uint64_t plaintext) {
    BN_CTX* ctx      = BN_CTX_new();
    BIGNUM *bn_plain = BN_new(), *ord = BN_new(), *rand = BN_new();

    EC_GROUP_get_order(init_group, ord, ctx);
    BN_rand_range(rand, ord);

    BN_set_word(bn_plain, plaintext);

    EC_POINT_mul(init_group, ciphertext->C1, NULL, EC_GROUP_get0_generator(init_group), rand, ctx);
    EC_POINT_mul(init_group, ciphertext->C2, bn_plain, pk, rand, ctx);

    BN_clear_free(rand);
    BN_free(ord);
    BN_free(bn_plain);
    BN_CTX_free(ctx);
    return 0;
}

// if table == NULL use bruteforce
int ec_elgamal_decrypt(uint64_t* res, const ec_elgamal_secret_key sk, const ec_elgamal_ciphertext* ciphertext,
                       bsgs_table_t table) {
    EC_POINT* M = EC_POINT_new(init_group);
    uint64_t plaintext;
    BN_CTX* ctx = BN_CTX_new();

    EC_POINT_mul(init_group, M, NULL, ciphertext->C1, sk, ctx);
    EC_POINT_invert(init_group, M, ctx);
    EC_POINT_add(init_group, M, ciphertext->C2, M, ctx);

    if (table != NULL) {
        ecdlp_bsgs(table, M, &plaintext, 1L << MAX_BITS);
    }
    else {
        ecdlp_brute(init_group, M, &plaintext, 1L << MAX_BITS);
    }
    *res = (uint64_t)plaintext;

    BN_CTX_free(ctx);
    EC_POINT_clear_free(M);
    return 0;
}

int ec_elgamal_add(ec_elgamal_ciphertext* res, const ec_elgamal_ciphertext* ciphertext1,
                   const ec_elgamal_ciphertext* ciphertext2) {
    BN_CTX* ctx = BN_CTX_new();
    EC_POINT_add(init_group, res->C1, ciphertext1->C1, ciphertext2->C1, ctx);
    EC_POINT_add(init_group, res->C2, ciphertext1->C2, ciphertext2->C2, ctx);
    BN_CTX_free(ctx);
    return 0;
}

int ec_elgamal_add_inplace(ec_elgamal_ciphertext* ciphertext1, const ec_elgamal_ciphertext* ciphertext2) {
    BN_CTX* ctx = BN_CTX_new();
    EC_POINT_add(init_group, ciphertext1->C1, ciphertext1->C1, ciphertext2->C1, ctx);
    EC_POINT_add(init_group, ciphertext1->C2, ciphertext1->C2, ciphertext2->C2, ctx);
    BN_CTX_free(ctx);
    return 0;
}

int ec_elgamal_add_plain(ec_elgamal_ciphertext* res, const ec_elgamal_ciphertext* ciphertext, const uint64_t m) {
    BN_CTX* ctx      = BN_CTX_new();
    BIGNUM* bn_plain = BN_new();
    BN_set_word(bn_plain, m);

    EC_POINT* mG = EC_POINT_new(init_group);
    EC_POINT_mul(init_group, mG, NULL, EC_GROUP_get0_generator(init_group), bn_plain, ctx);
    EC_POINT_copy(res->C1, ciphertext->C1);
    EC_POINT_add(init_group, res->C2, ciphertext->C2, mG, ctx);

    EC_POINT_free(mG);
    BN_free(bn_plain);
    BN_CTX_free(ctx);
    return 0;
}

int ec_elgamal_add_plain_inplace(ec_elgamal_ciphertext* ciphertext, const uint64_t m) {
    BN_CTX* ctx      = BN_CTX_new();
    BIGNUM* bn_plain = BN_new();
    BN_set_word(bn_plain, m);

    EC_POINT* mG = EC_POINT_new(init_group);
    EC_POINT_mul(init_group, mG, NULL, EC_GROUP_get0_generator(init_group), bn_plain, ctx);
    EC_POINT_add(init_group, ciphertext->C2, ciphertext->C2, mG, ctx);

    EC_POINT_free(mG);
    BN_free(bn_plain);
    BN_CTX_free(ctx);
    return 0;
}

int ec_elgamal_mul_plain(ec_elgamal_ciphertext* res, const ec_elgamal_ciphertext* ciphertext, const uint64_t m) {
    BN_CTX* ctx      = BN_CTX_new();
    BIGNUM* bn_plain = BN_new();
    BN_set_word(bn_plain, m);

    EC_POINT_mul(init_group, res->C1, NULL, ciphertext->C1, bn_plain, ctx);
    EC_POINT_mul(init_group, res->C2, NULL, ciphertext->C2, bn_plain, ctx);

    BN_free(bn_plain);
    BN_CTX_free(ctx);

    return 0;
}

int ec_elgamal_mul_plain_inplace(ec_elgamal_ciphertext* ciphertext, const uint64_t m) {
    BN_CTX* ctx      = BN_CTX_new();
    BIGNUM* bn_plain = BN_new();
    BN_set_word(bn_plain, m);

    EC_POINT_mul(init_group, ciphertext->C1, NULL, ciphertext->C1, bn_plain, ctx);
    EC_POINT_mul(init_group, ciphertext->C2, NULL, ciphertext->C2, bn_plain, ctx);

    BN_free(bn_plain);
    BN_CTX_free(ctx);

    return 0;
}

EC_GROUP* ec_elgamal_get_current_group() {
    return init_group;
}

int ec_elgamal_get_point_compressed_size() {
    BN_CTX* ctx = BN_CTX_new();
    int res     = (int)EC_POINT_point2oct(init_group, EC_GROUP_get0_generator(init_group), POINT_CONVERSION_COMPRESSED,
                                          NULL, 0, ctx);
    BN_CTX_free(ctx);
    return res;
}

// ENCODING + DECODING

void write_size(unsigned char* buffer, size_t size) {
    buffer[0] = (unsigned char)((size >> 8) & 0xFF);
    buffer[1] = (unsigned char)(size & 0xFF);
}

size_t read_size(unsigned char* buffer) {
    return ((uint8_t)buffer[0] << 8) | ((uint8_t)buffer[1]);
}

size_t get_encoded_ciphertext_size(ec_elgamal_ciphertext* ciphertext) {
    return (size_t)ec_elgamal_get_point_compressed_size() * 2;
}

int encode_ciphertext(unsigned char* buff, ec_elgamal_ciphertext* ciphertext) {
    unsigned char* cur_ptr = buff;
    size_t len_point, tmp;
    BN_CTX* ctx = BN_CTX_new();
    len_point   = (size_t)ec_elgamal_get_point_compressed_size();
    tmp         = EC_POINT_point2oct(init_group, ciphertext->C1, POINT_CONVERSION_COMPRESSED, cur_ptr, len_point, ctx);
    cur_ptr += len_point;
    if (tmp != len_point)
        return -1;
    tmp = EC_POINT_point2oct(init_group, ciphertext->C2, POINT_CONVERSION_COMPRESSED, cur_ptr, len_point, ctx);
    if (tmp != len_point)
        return -1;
    BN_CTX_free(ctx);
    return 0;
}

int decode_ciphertext(ec_elgamal_ciphertext* ciphertext, unsigned char* buff) {
    size_t len_point;
    BN_CTX* ctx            = BN_CTX_new();
    unsigned char* cur_ptr = buff;
    len_point              = (size_t)ec_elgamal_get_point_compressed_size();

    EC_POINT_oct2point(init_group, ciphertext->C1, cur_ptr, len_point, ctx);
    cur_ptr += len_point;

    EC_POINT_oct2point(init_group, ciphertext->C2, cur_ptr, len_point, ctx);

    BN_CTX_free(ctx);
    return 0;
}