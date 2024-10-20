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

#ifndef PRIV_LR_ECELGAMAL_H
#define PRIV_LR_ECELGAMAL_H

#include <openssl/ec.h>
#include <openssl/bn.h>
#include <openssl/objects.h>
#include <inttypes.h>
#include "uthash.h"

#define DEFAULT_CURVE NID_X9_62_prime192v1
#define CURVE_256_SEC NID_X9_62_prime256v1

#define MAX_BITS 32

typedef BIGNUM* ec_elgamal_secret_key;
typedef EC_POINT* ec_elgamal_public_key;

int ec_elgamal_free_public_key(ec_elgamal_public_key key);
int ec_elgamal_free_secret_key(ec_elgamal_secret_key key);

typedef struct {
    EC_POINT* C1;
    EC_POINT* C2;
} ec_elgamal_ciphertext;

ec_elgamal_ciphertext* ec_elgamal_new_ciphertext();
int ec_elgamal_free_ciphertext(ec_elgamal_ciphertext* ciphertext);
size_t get_encoded_ciphertext_size(ec_elgamal_ciphertext* ciphertext);
int encode_ciphertext(unsigned char* buff, ec_elgamal_ciphertext* ciphertext);
int decode_ciphertext(ec_elgamal_ciphertext* ciphertext, unsigned char* buff);

typedef struct bsgs_hash_table_entry {
    unsigned char* key;
    uint32_t value;
    UT_hash_handle hh;
} bsgs_hash_table_entry_t;

struct bsgs_table_s {
    bsgs_hash_table_entry_t* table;
    EC_POINT* mG;
    EC_POINT* mG_inv;
    EC_GROUP* group;
    uint64_t tablesize;
};
typedef struct bsgs_table_s* bsgs_table_ptr;
typedef struct bsgs_table_s bsgs_table_t[1];

int ecdlp_bsgs(const bsgs_table_t bsgs_table, const EC_POINT* M, uint64_t* x, int64_t max_it);
int ecdlp_brute(const EC_GROUP* curve_group, const EC_POINT* M, uint64_t* x, uint64_t max_it);

/**
 * Inits the library with the given curve
 * @param curve_id
 * @return
 */
int ec_elgamal_init(int curve_id);

const EC_GROUP* ec_elgamal_get_group();

/**
 * Deinits the library
 * @return
 */
int ec_elgamal_deinit();

/**
 * Inititlaizes the baby-step-giant-step table.
 * @param table
 * @param size number of elemnts to store in the table
 * @return
 */
int bsgs_table_init(bsgs_table_t table, size_t size);

/**
 * Frees the memory of the table
 * @param table
 * @return
 */
int bsgs_table_free(bsgs_table_t table);

/**
 * Generates an EC-Elgamal secret key
 * @return the secret key
 */
ec_elgamal_secret_key ec_elgamal_keygen();

/**
 * Generates an EC-Elgamal public key
 * @param sk the secret key
 * @return the public key
 */
ec_elgamal_public_key ec_elgamal_from_secret_key(const ec_elgamal_secret_key sk);

/**
 * Returns the EC_Group (elliptic curve group) struct if initialized
 */
EC_GROUP* ec_elgamal_get_current_group();

/**
 * Returns the encded size of an EC-Point in this group.
 */
int ec_elgamal_get_point_compressed_size();

/**
 * Encrypts an Integer with additadive homomorphic EC-ElGamal
 * @param ciphertext
 * @param key
 * @param plaintext
 * @return
 */
int ec_elgamal_encrypt(ec_elgamal_ciphertext* ciphertext, const ec_elgamal_public_key pk, const uint64_t plaintext);

/**
 * Decrypts an EC-Elgamal ciphertext
 * @param res the resulting plaintext integer
 * @param key
 * @param ciphertext
 * @param table if NULL bruteforce is used
 * @return
 */
int ec_elgamal_decrypt(uint64_t* res, const ec_elgamal_secret_key sk, const ec_elgamal_ciphertext* ciphertext,
                       bsgs_table_t table);

/**
 * Adds two EC-Elgamal ciphertext and stores it in res.
 * @param res the resulting ciphertext
 * @param ciphertext1
 * @param ciphertext2
 * @return
 */
int ec_elgamal_add(ec_elgamal_ciphertext* res, const ec_elgamal_ciphertext* ciphertext1,
                   const ec_elgamal_ciphertext* ciphertext2);

/**
 * Adds two EC-Elgamal ciphertext and stores it in ciphertext1.
 * @param ciphertext1
 * @param ciphertext2
 * @return
 */
int ec_elgamal_add_inplace(ec_elgamal_ciphertext* ciphertext1, const ec_elgamal_ciphertext* ciphertext2);

/**
 * Adds EC-Elgamal ciphertext and plaintext, and stores it in res.
 * @param res the resulting ciphertext
 * @param ciphertext
 * @param m
 * @return
 */
int ec_elgamal_add_plain(ec_elgamal_ciphertext* res, const ec_elgamal_ciphertext* ciphertext, const uint64_t m);

/**
 * Adds EC-Elgamal ciphertext and plaintext, and stores it in ciphertext.
 * @param ciphertext
 * @param m
 * @return
 */
int ec_elgamal_add_plain_inplace(ec_elgamal_ciphertext* ciphertext, const uint64_t m);

/**
 * Muls EC-Elgamal ciphertext and plaintext, and stores it in res.
 * @param res the resulting ciphertext
 * @param ciphertext
 * @param m
 * @return
 */
int ec_elgamal_mul_plain(ec_elgamal_ciphertext* res, const ec_elgamal_ciphertext* ciphertext, const uint64_t m);

/**
 * Muls EC-Elgamal ciphertext and plaintext, and stores it in ciphertext.
 * @param ciphertext
 * @param m
 * @return
 */
int ec_elgamal_mul_plain_inplace(ec_elgamal_ciphertext* ciphertext, const uint64_t m);

#endif  // PRIV_LR_ECELGAMAL_H
