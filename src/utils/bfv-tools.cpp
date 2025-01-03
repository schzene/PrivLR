#include "utils/bfv-tools.h"

#define MAX_SZ 10

void print_parameters(std::shared_ptr<seal::SEALContext> context) {
    // Verify parameters
    if (!context) {
        throw std::invalid_argument("context is not set");
    }
    auto& context_data = *context->key_context_data();

    /*
    Which scheme are we using?
    */
    std::string scheme_name;
    switch (context_data.parms().scheme()) {
        case seal::scheme_type::bfv:
            scheme_name = "BFV";
            break;
        case seal::scheme_type::ckks:
            scheme_name = "CKKS";
            break;
        default:
            throw std::invalid_argument("unsupported scheme");
    }
    std::cout << "/" << std::endl;
    std::cout << "| Encryption parameters :" << std::endl;
    std::cout << "|   scheme: " << scheme_name << std::endl;
    std::cout << "|   poly_modulus_degree: " << context_data.parms().poly_modulus_degree() << std::endl;

    /*
    Print the size of the true (product) coefficient modulus.
    */
    std::cout << "|   coeff_modulus size: ";
    std::cout << context_data.total_coeff_modulus_bit_count() << " (";
    auto coeff_modulus          = context_data.parms().coeff_modulus();
    std::size_t coeff_mod_count = coeff_modulus.size();
    for (std::size_t i = 0; i < coeff_mod_count - 1; i++) {
        std::cout << coeff_modulus[i].bit_count() << " + ";
    }
    std::cout << coeff_modulus.back().bit_count();
    std::cout << ") bits" << std::endl;

    /*
    For the BFV scheme print the plain_modulus parameter.
    */
    if (context_data.parms().scheme() == seal::scheme_type::bfv) {
        std::cout << "|   plain_modulus: " << context_data.parms().plain_modulus().value() << std::endl;
    }

    std::cout << "\\" << std::endl;
}

BFVParm::BFVParm(size_t poly_modulus_degree_, uint64_t plain_mod_)
    : poly_modulus_degree(poly_modulus_degree_), slot_count(poly_modulus_degree_), plain_mod(plain_mod_) {
    // Generate keys
    EncryptionParameters parms(scheme_type::bfv);

    parms.set_poly_modulus_degree(poly_modulus_degree);
    // parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_bit_sizes_));
    parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    parms.set_plain_modulus(plain_mod);

    context   = new SEALContext(parms, true, seal::sec_level_type::tc128);
    encoder   = new BatchEncoder(*context);
    evaluator = new Evaluator(*context);
}

BFVParm::~BFVParm() {
    delete context;
    delete encoder;
    delete evaluator;
}

BFVKey::BFVKey(int party_, BFVParm* parm_) : party(party_), parm(parm_) {
    assert(party == ALICE || party == BOB);

    KeyGenerator* keygen = new KeyGenerator(*(parm->context));
    keygen->create_public_key(public_key);
    keygen->create_relin_keys(relin_keys);
    keygen->create_galois_keys(galois_keys);
    encryptor = new Encryptor(*(parm->context), public_key);
    decryptor = new Decryptor(*(parm->context), keygen->secret_key());
    delete keygen;
}

BFVKey::~BFVKey() {
    delete encryptor;
    delete decryptor;
}

BFVLongPlaintext::BFVLongPlaintext(const Plaintext& pt) {
    len = 1;
    plain_data.push_back(pt);
}

/*
     BFVLongPlaintext:  for field [0,2p]
*/
BFVLongPlaintext::BFVLongPlaintext(BFVParm* parm, uint64_t data) {
    // TODO value len =1, use the BFV batchencoder to encode the palaintext
    len = 1;
    Plaintext pt;
    vector<uint64_t> temp(parm->slot_count, data);
    parm->encoder->encode(temp, pt);
    plain_data.push_back(pt);
}

BFVLongPlaintext::BFVLongPlaintext(BFVParm* parm, vector<uint64_t> data) {
    len               = data.size();
    size_t slot_count = parm->slot_count;  // TODO:: this slot_count use SEALcontext? BFVLongPlaintext contain it.
    size_t count      = len / slot_count;

    if (len % slot_count) {
        count++;
    }
    size_t i, j;
    if (slot_count >= len) {
        Plaintext pt;
        parm->encoder->encode(data, pt);
        plain_data.push_back(pt);
    }
    else {
        vector<uint64_t>::iterator curPtr = data.begin(), endPtr = data.end(), end;
        while (curPtr < endPtr) {
            end        = endPtr - curPtr > slot_count ? slot_count + curPtr : endPtr;
            slot_count = endPtr - curPtr > slot_count ? slot_count : endPtr - curPtr;
            vector<uint64_t> temp(curPtr, end);
            Plaintext pt;
            parm->encoder->encode(temp, pt);
            plain_data.push_back(pt);
            curPtr += slot_count;
        }
    }
}

BFVLongPlaintext::BFVLongPlaintext(BFVParm* parm, uint64_t* data, size_t len) {
    this->len         = len;
    size_t slot_count = parm->slot_count;  // TODO:: this slot_count use SEALcontext? BFVLongPlaintext contain it.
    size_t count      = len / slot_count;

    if (len % slot_count) {
        count++;
    }
    size_t i, j;
    if (slot_count >= len) {
        Plaintext pt;
        parm->encoder->encode(vector<uint64_t>(data, data + len), pt);
        plain_data.push_back(pt);
    }
    else {
        uint64_t *curPtr = data, *endPtr = data + len, *end;
        while (curPtr < endPtr) {
            end        = endPtr - curPtr > slot_count ? slot_count + curPtr : endPtr;
            slot_count = endPtr - curPtr > slot_count ? slot_count : endPtr - curPtr;
            vector<uint64_t> temp(curPtr, end);
            Plaintext pt;
            parm->encoder->encode(temp, pt);
            plain_data.push_back(pt);
            curPtr += slot_count;
        }
    }
}

vector<uint64_t> BFVLongPlaintext::decode_uint(BFVParm* parm) const {
    vector<uint64_t> data(len);  // 5

    size_t size = plain_data.size();  // 1

    size_t solut_cout = parm->slot_count;  // 8192

    for (size_t i = 0; i < size; i++) {
        vector<uint64_t> temp;
        parm->encoder->decode(plain_data[i], temp);

        if (i < size - 1) {
            copy(temp.begin(), temp.end(), data.begin() + i * solut_cout);
        }
        else {
            size_t tail_len = len % solut_cout;
            tail_len        = tail_len ? tail_len : solut_cout;
            copy(temp.begin(), temp.begin() + tail_len, data.begin() + i * solut_cout);
        }
    }
    return data;
}

/*
     BFVLongPlaintext: for field [-p,p]
*/
BFVLongPlaintext::BFVLongPlaintext(BFVParm* parm, int64_t data) {
    // TODO value len =1, use the BFV batchencoder to encode the palaintext
    len = 1;
    Plaintext pt;
    vector<int64_t> temp(parm->slot_count, data);
    parm->encoder->encode(temp, pt);
    plain_data.push_back(pt);
}

BFVLongPlaintext::BFVLongPlaintext(BFVParm* parm, vector<int64_t> data) {
    len               = data.size();
    size_t slot_count = parm->slot_count;  // TODO:: this slot_count use SEALcontext? BFVLongPlaintext contain it.
    size_t count      = len / slot_count;

    if (len % slot_count) {
        count++;
    }
    size_t i, j;
    if (slot_count >= len) {
        Plaintext pt;
        parm->encoder->encode(data, pt);
        plain_data.push_back(pt);
    }
    else {
        vector<int64_t>::iterator curPtr = data.begin(), endPtr = data.end(), end;
        while (curPtr < endPtr) {
            end        = endPtr - curPtr > slot_count ? slot_count + curPtr : endPtr;
            slot_count = endPtr - curPtr > slot_count ? slot_count : endPtr - curPtr;
            vector<int64_t> temp(curPtr, end);
            Plaintext pt;
            parm->encoder->encode(temp, pt);
            plain_data.push_back(pt);
            curPtr += slot_count;
        }
    }
}

BFVLongPlaintext::BFVLongPlaintext(BFVParm* parm, int64_t* data, size_t len) {
    this->len         = len;
    size_t slot_count = parm->slot_count;  // TODO:: this slot_count use SEALcontext? BFVLongPlaintext contain it.
    size_t count      = len / slot_count;

    if (len % slot_count) {
        count++;
    }
    size_t i, j;
    if (slot_count >= len) {
        Plaintext pt;
        parm->encoder->encode(vector<int64_t>(data, data + len), pt);
        plain_data.push_back(pt);
    }
    else {
        int64_t *curPtr = data, *endPtr = data + len, *end;
        while (curPtr < endPtr) {
            end        = endPtr - curPtr > slot_count ? slot_count + curPtr : endPtr;
            slot_count = endPtr - curPtr > slot_count ? slot_count : endPtr - curPtr;
            vector<int64_t> temp(curPtr, end);
            Plaintext pt;
            parm->encoder->encode(temp, pt);
            plain_data.push_back(pt);
            curPtr += slot_count;
        }
    }
}

vector<int64_t> BFVLongPlaintext::decode_int(BFVParm* parm) const {
    vector<int64_t> data(len);  // 5

    size_t size = plain_data.size();  // 1

    size_t solut_cout = parm->slot_count;  // 8192

    for (size_t i = 0; i < size; i++) {
        vector<int64_t> temp;
        parm->encoder->decode(plain_data[i], temp);

        if (i < size - 1) {
            copy(temp.begin(), temp.end(), data.begin() + i * solut_cout);
        }
        else {
            size_t tail_len = len % solut_cout;
            tail_len        = tail_len ? tail_len : solut_cout;
            copy(temp.begin(), temp.begin() + tail_len, data.begin() + i * solut_cout);
        }
    }
    return data;
}

/*
     BFVLongPlaintext:END
*/

BFVLongCiphertext::BFVLongCiphertext(const Ciphertext& ct) {
    len = 1;
    cipher_data.push_back(ct);
}

/*
     BFVLongCiphertext: for field [0,2p]
*/

BFVLongCiphertext::BFVLongCiphertext(uint64_t data, const BFVKey* party) {
    // TODO:
    len = 1;
    Plaintext pt;
    vector<uint64_t> temp(party->parm->slot_count, data);
    party->parm->encoder->encode(temp, pt);
    Ciphertext ct;
    party->encryptor->encrypt(pt, ct);
    cipher_data.push_back(ct);
}

BFVLongCiphertext::BFVLongCiphertext(uint64_t* data, size_t len, const BFVKey* party) {
    this->len = len;
    size_t slot_count =
        party->parm->slot_count;  // TODO:: this slot_count use SEALcontext? BFVLongPlaintext contain it.
    size_t count = len / slot_count;

    if (len % slot_count) {
        count++;
    }
    size_t i, j;
    if (slot_count >= len) {
        Plaintext pt;
        Ciphertext ct;
        party->parm->encoder->encode(vector<uint64_t>(data, data + len), pt);
        party->encryptor->encrypt(pt, ct);
        cipher_data.push_back(ct);
    }
    else {
        uint64_t *curPtr = data, *endPtr = data + len, *end;
        while (curPtr < endPtr) {
            end        = endPtr - curPtr > slot_count ? slot_count + curPtr : endPtr;
            slot_count = endPtr - curPtr > slot_count ? slot_count : endPtr - curPtr;
            vector<uint64_t> temp(curPtr, end);
            Plaintext pt;
            Ciphertext ct;
            party->parm->encoder->encode(temp, pt);
            party->encryptor->encrypt(pt, ct);
            cipher_data.push_back(ct);
            curPtr += slot_count;
        }
    }
}

/*
     BFVLongCiphertext: for field [-p,p]
*/

BFVLongCiphertext::BFVLongCiphertext(int64_t data, const BFVKey* party) {
    // TODO:
    len = 1;
    Plaintext pt;
    vector<int64_t> temp(party->parm->slot_count, data);
    party->parm->encoder->encode(temp, pt);
    Ciphertext ct;
    party->encryptor->encrypt(pt, ct);
    cipher_data.push_back(ct);
}

BFVLongCiphertext::BFVLongCiphertext(int64_t* data, size_t len, const BFVKey* party) {
    this->len = len;
    size_t slot_count =
        party->parm->slot_count;  // TODO:: this slot_count use SEALcontext? BFVLongPlaintext contain it.
    size_t count = len / slot_count;

    if (len % slot_count) {
        count++;
    }
    size_t i, j;
    if (slot_count >= len) {
        Plaintext pt;
        Ciphertext ct;
        party->parm->encoder->encode(vector<int64_t>(data, data + len), pt);
        party->encryptor->encrypt(pt, ct);
        cipher_data.push_back(ct);
    }
    else {
        int64_t *curPtr = data, *endPtr = data + len, *end;
        while (curPtr < endPtr) {
            end        = endPtr - curPtr > slot_count ? slot_count + curPtr : endPtr;
            slot_count = endPtr - curPtr > slot_count ? slot_count : endPtr - curPtr;
            vector<int64_t> temp(curPtr, end);
            Plaintext pt;
            Ciphertext ct;
            party->parm->encoder->encode(temp, pt);
            party->encryptor->encrypt(pt, ct);
            cipher_data.push_back(ct);
            curPtr += slot_count;
        }
    }
}

/*
     BFVLongCiphertext:END
*/

BFVLongCiphertext::BFVLongCiphertext(const BFVLongPlaintext& lpt, const BFVKey* party) {
    len         = lpt.len;
    size_t size = lpt.plain_data.size();
    cipher_data.resize(size);
#pragma omp parallel for if (size > 1)
    for (size_t i = 0; i < size; i++) {
        party->encryptor->encrypt(lpt.plain_data[i], cipher_data[i]);
    }
}

BFVLongPlaintext BFVLongCiphertext::decrypt(const BFVKey* party) const {
    BFVLongPlaintext lpt;
    lpt.len     = len;
    size_t size = cipher_data.size();
    lpt.plain_data.resize(size);
#pragma omp parallel for if (size > 1)
    for (size_t i = 0; i < size; i++) {
        party->decryptor->decrypt(cipher_data[i], lpt.plain_data[i]);
    }
    return lpt;
}

void BFVLongCiphertext::add_plain_inplace(BFVLongPlaintext& lpt, Evaluator* evaluator) {
    if (len == 1)  // cipher text len =1
    {
        len = lpt.len;
        Ciphertext ct(cipher_data[0]);
        cipher_data.pop_back();
        size_t size = lpt.plain_data.size();
        cipher_data = vector<Ciphertext>(size);
        for (size_t i = 0; i < size; i++) {
            Ciphertext ctemp;
            evaluator->add_plain(ct, lpt.plain_data[i], ctemp);
            cipher_data[i] = ctemp;
        }
    }
    else if (lpt.len == 1) {
        for (size_t i = 0; i < cipher_data.size(); i++) {
            evaluator->add_plain_inplace(cipher_data[i], lpt.plain_data[0]);
        }
    }
    else if (len == lpt.len) {
        for (size_t i = 0; i < cipher_data.size(); i++) {
            evaluator->add_plain_inplace(cipher_data[i], lpt.plain_data[i]);
        }
    }
    else {
        char buf[100];
        sprintf(buf, "In add_plain_inplace: Length of BFVLongCiphertext(%ld) and BFVLongPlaintext(%ld) mismatch", len,
                lpt.len);
        throw bfv_lenth_error(buf);
    }
}

BFVLongCiphertext BFVLongCiphertext::add_plain(BFVLongPlaintext& lpt, Evaluator* evaluator) const {
    BFVLongCiphertext lct;
    lct.len = 0;
    if (len == 1) {
        lct.len = lpt.len;
        for (size_t i = 0; i < lpt.plain_data.size(); i++) {
            Ciphertext ct;
            evaluator->add_plain(cipher_data[0], lpt.plain_data[i], ct);
            lct.cipher_data.push_back(ct);
        }
    }
    else if (lpt.len == 1) {
        lct.len = len;
        for (size_t i = 0; i < cipher_data.size(); i++) {
            Ciphertext ct;
            evaluator->add_plain(cipher_data[i], lpt.plain_data[0], ct);
            lct.cipher_data.push_back(ct);
        }
    }
    else if (len == lpt.len) {
        lct.len = len;
        for (size_t i = 0; i < cipher_data.size(); i++) {
            Ciphertext ct;
            evaluator->add_plain(cipher_data[i], lpt.plain_data[i], ct);
            lct.cipher_data.push_back(ct);
        }
    }
    else {
        char buf[100];
        sprintf(buf, "In add_plain: Length of BFVLongCiphertext(%ld) and LongPlaintext(%ld) mismatch", len, lpt.len);
        throw bfv_lenth_error(buf);
    }
    return lct;
}

void BFVLongCiphertext::add_inplace(BFVLongCiphertext& lct, Evaluator* evaluator) {
    if (len == 1) {
        len = lct.len;
        Ciphertext ct(cipher_data[0]);
        cipher_data.pop_back();
        for (Ciphertext cct : lct.cipher_data) {
            Ciphertext ctemp;
            evaluator->add(ct, cct, ctemp);
            cipher_data.push_back(ctemp);
        }
    }
    else if (lct.len == 1) {
        for (size_t i = 0; i < cipher_data.size(); i++)
            evaluator->add_inplace(cipher_data[i], lct.cipher_data[0]);
    }
    else if (len == lct.len) {
        for (size_t i = 0; i < cipher_data.size(); i++)
            evaluator->add_inplace(cipher_data[i], lct.cipher_data[i]);
    }
    else {
        char buf[100];
        sprintf(buf, "In add_inplace: Length of BFVLongCiphertext(%ld) and BFVLongCiphertext(%ld) mismatch", len,
                lct.len);
        throw bfv_lenth_error(buf);
    }
}

BFVLongCiphertext BFVLongCiphertext::add(BFVLongCiphertext& lct, Evaluator* evaluator) const {
    BFVLongCiphertext lcct;
    lcct.len = 0;
    if (len == 1) {
        lcct.len = lct.len;
        for (size_t i = 0; i < lct.cipher_data.size(); i++) {
            Ciphertext ct;
            evaluator->add(cipher_data[0], lct.cipher_data[i], ct);
            lcct.cipher_data.push_back(ct);
        }
    }
    else if (lct.len == 1) {
        lcct.len = len;
        for (size_t i = 0; i < cipher_data.size(); i++) {
            Ciphertext ct;
            evaluator->add(cipher_data[i], lct.cipher_data[0], ct);
            lcct.cipher_data.push_back(ct);
        }
    }
    else if (len == lct.len) {
        lcct.len = len;
        for (size_t i = 0; i < cipher_data.size(); i++) {
            Ciphertext ct;
            evaluator->add(cipher_data[i], lct.cipher_data[i], ct);
            lcct.cipher_data.push_back(ct);
        }
    }
    else {
        char buf[100];
        sprintf(buf, "In add: Length of BFVLongCiphertext(%ld) and BFVLongCiphertext(%ld) mismatch", len, lct.len);
        throw bfv_lenth_error(buf);
    }
    return lcct;
}

void BFVLongCiphertext::sub_plain_inplace(BFVLongPlaintext& lpt, Evaluator* evaluator) {
    if (len == 1)  // cipher text len =1
    {
        len = lpt.len;
        Ciphertext ct(cipher_data[0]);
        cipher_data.pop_back();
        for (Plaintext pt : lpt.plain_data) {
            Ciphertext ctemp;
            evaluator->sub_plain(ct, pt, ctemp);
            cipher_data.push_back(ctemp);
        }
    }
    else if (lpt.len == 1) {
        for (size_t i = 0; i < cipher_data.size(); i++)
            evaluator->sub_plain_inplace(cipher_data[i], lpt.plain_data[0]);
    }
    else if (len == lpt.len) {
        for (size_t i = 0; i < cipher_data.size(); i++)
            evaluator->sub_plain_inplace(cipher_data[i], lpt.plain_data[i]);
    }
    else {
        char buf[100];
        sprintf(buf, "In sub_plain_inplace: Length of BFVLongCiphertext(%ld) and BFVLongPlaintext(%ld) mismatch", len,
                lpt.len);
        throw bfv_lenth_error(buf);
    }
}

BFVLongCiphertext BFVLongCiphertext::sub_plain(BFVLongPlaintext& lpt, Evaluator* evaluator) const {
    {
        BFVLongCiphertext lct;
        lct.len = 0;
        if (len == 1) {
            lct.len = lpt.len;
            for (size_t i = 0; i < lpt.plain_data.size(); i++) {
                Ciphertext ct;
                evaluator->sub_plain(cipher_data[0], lpt.plain_data[i], ct);
                lct.cipher_data.push_back(ct);
            }
        }
        else if (lpt.len == 1) {
            lct.len = len;
            for (size_t i = 0; i < cipher_data.size(); i++) {
                Ciphertext ct;
                evaluator->sub_plain(cipher_data[i], lpt.plain_data[0], ct);
                lct.cipher_data.push_back(ct);
            }
        }
        else if (len == lpt.len) {
            lct.len = len;
            for (size_t i = 0; i < cipher_data.size(); i++) {
                Ciphertext ct;
                evaluator->sub_plain(cipher_data[i], lpt.plain_data[i], ct);
                lct.cipher_data.push_back(ct);
            }
        }
        else {
            char buf[100];
            sprintf(buf, "In sub_plain: Length of BFVLongCiphertext(%ld) and LongPlaintext(%ld) mismatch", len,
                    lpt.len);
            throw bfv_lenth_error(buf);
        }
        return lct;
    }
}

void BFVLongCiphertext::sub_inplace(BFVLongCiphertext& lct, Evaluator* evaluator) {
    if (len == 1) {
        len = lct.len;
        Ciphertext ct(cipher_data[0]);
        cipher_data.pop_back();
        for (Ciphertext cct : lct.cipher_data) {
            Ciphertext ctemp;
            evaluator->sub(ct, cct, ctemp);
            cipher_data.push_back(ctemp);
        }
    }
    else if (lct.len == 1) {
        for (size_t i = 0; i < cipher_data.size(); i++)
            evaluator->sub_inplace(cipher_data[i], lct.cipher_data[0]);
    }
    else if (len == lct.len) {
        for (size_t i = 0; i < cipher_data.size(); i++)
            evaluator->sub_inplace(cipher_data[i], lct.cipher_data[i]);
    }
    else {
        char buf[100];
        sprintf(buf, "In sub_inplace: Length of BFVLongCiphertext(%ld) and BFVLongCiphertext(%ld) mismatch", len,
                lct.len);
        throw bfv_lenth_error(buf);
    }
}

BFVLongCiphertext BFVLongCiphertext::sub(BFVLongCiphertext& lct, Evaluator* evaluator) const {
    BFVLongCiphertext lcct;
    lcct.len = 0;
    if (len == 1) {
        lcct.len = lct.len;
        for (size_t i = 0; i < lct.cipher_data.size(); i++) {
            Ciphertext ct;
            evaluator->sub(cipher_data[0], lct.cipher_data[i], ct);
            lcct.cipher_data.push_back(ct);
        }
    }
    else if (lct.len == 1) {
        lcct.len = len;
        for (size_t i = 0; i < cipher_data.size(); i++) {
            Ciphertext ct;
            evaluator->sub(cipher_data[i], lct.cipher_data[0], ct);
            lcct.cipher_data.push_back(ct);
        }
    }
    else if (len == lct.len) {
        lcct.len = len;
        for (size_t i = 0; i < cipher_data.size(); i++) {
            Ciphertext ct;
            evaluator->sub(cipher_data[i], lct.cipher_data[i], ct);
            lcct.cipher_data.push_back(ct);
        }
    }
    else {
        char buf[100];
        sprintf(buf, "In sub: Length of BFVLongCiphertext(%ld) and BFVLongCiphertext(%ld) mismatch", len, lct.len);
        throw bfv_lenth_error(buf);
    }
    return lcct;
}

void BFVLongCiphertext::multiply_plain_inplace(BFVLongPlaintext& lpt, Evaluator* evaluator) {
    size_t size = cipher_data.size();
    if (len == 1) {
        len  = lpt.len;
        size = lpt.plain_data.size();
        Ciphertext ct(cipher_data[0]);
        cipher_data.resize(size);
#pragma omp parallel for if (size > MAX_SZ)
        for (size_t i = 0; i < size; i++) {
            Ciphertext ctemp;
            evaluator->multiply_plain(ct, lpt.plain_data[i], ctemp);
            cipher_data[i] = ctemp;
        }
    }
    else if (lpt.len == 1) {
#pragma omp parallel for if (size > MAX_SZ)
        for (size_t i = 0; i < size; i++) {
            evaluator->multiply_plain_inplace(cipher_data[i], lpt.plain_data[0]);
        }
    }
    else if (len == lpt.len) {
#pragma omp parallel for if (size > MAX_SZ)
        for (size_t i = 0; i < size; i++) {
            evaluator->multiply_plain_inplace(cipher_data[i], lpt.plain_data[i]);
        }
    }
    else {
        char buf[100];
        sprintf(buf, "In multiply_plain_inplace: Length of LongCiphertext(%ld) and LongPlaintext(%ld) mismatch", len,
                lpt.len);
        throw bfv_lenth_error(buf);
    }
}

BFVLongCiphertext BFVLongCiphertext::multiply_plain(BFVLongPlaintext& lpt, Evaluator* evaluator) const {
    BFVLongCiphertext lct;
    lct.len     = 0;
    size_t size = cipher_data.size();
    if (len == 1) {
        lct.len = lpt.len;
        size    = lpt.plain_data.size();
        lct.cipher_data.resize(size);
        // #pragma omp parallel for if (size > MAX_SZ)
        for (size_t i = 0; i < size; i++) {
            evaluator->multiply_plain(cipher_data[0], lpt.plain_data[i], lct.cipher_data[i]);
        }
    }
    else if (lpt.len == 1) {
        lct.len = len;
        lct.cipher_data.resize(size);
        // #pragma omp parallel for if (size > MAX_SZ)
        for (size_t i = 0; i < size; i++) {
            evaluator->multiply_plain(cipher_data[i], lpt.plain_data[0], lct.cipher_data[i]);
        }
    }
    else if (len == lpt.len) {
        lct.len = len;
        lct.cipher_data.resize(size);
        // #pragma omp parallel for if (size > MAX_SZ)
        for (size_t i = 0; i < size; i++) {
            evaluator->multiply_plain(cipher_data[i], lpt.plain_data[i], lct.cipher_data[i]);
        }
    }
    else {
        char buf[100];
        sprintf(buf, "In multiply_plain: Length of BFVLongCiphertext(%ld) and BFVLongPlaintext(%ld) mismatch", len,
                lpt.len);
        throw bfv_lenth_error(buf);
    }
    return lct;
}

void BFVLongCiphertext::multiply_inplace(BFVLongCiphertext& lct, Evaluator* evaluator) {
    if (len == 1) {
        len = lct.len;
        Ciphertext ct(cipher_data[0]);
        cipher_data.pop_back();
        for (Ciphertext cct : lct.cipher_data) {
            Ciphertext ctemp;
            evaluator->multiply(ct, cct, ctemp);
            cipher_data.push_back(ctemp);
        }
    }
    else if (lct.len == 1) {
        for (size_t i = 0; i < cipher_data.size(); i++)
            evaluator->multiply_inplace(cipher_data[i], lct.cipher_data[0]);
    }
    else if (len == lct.len) {
        for (size_t i = 0; i < cipher_data.size(); i++)
            evaluator->multiply_inplace(cipher_data[i], lct.cipher_data[i]);
    }
    else {
        char buf[100];
        sprintf(buf, "In multiply_inplace: Length of BFVLongCiphertext(%ld) and BFVLongCiphertext(%ld) mismatch", len,
                lct.len);
        throw bfv_lenth_error(buf);
    }
}

BFVLongCiphertext BFVLongCiphertext::multiply(BFVLongCiphertext& lct, Evaluator* evaluator) const {
    BFVLongCiphertext ret;
    ret.len     = len;
    size_t size = cipher_data.size();
    if (len == 1) {
        ret.len = lct.len;
        size    = lct.cipher_data.size();
        ret.cipher_data.resize(size);
#pragma omp parallel for if (size > 1)
        for (size_t i = 0; i < size; i++) {
            Ciphertext ct;
            evaluator->multiply(cipher_data[0], lct.cipher_data[i], ct);
            ret.cipher_data[i] = ct;
        }
    }
    else if (lct.len == 1) {
        ret.cipher_data.resize(size);
#pragma omp parallel for if (size > 1)
        for (size_t i = 0; i < size; i++) {
            Ciphertext ct;
            evaluator->multiply(cipher_data[i], lct.cipher_data[0], ct);
            ret.cipher_data[i] = ct;
        }
    }
    else if (len == lct.len) {
        ret.cipher_data.resize(size);
#pragma omp parallel for if (size > 1)
        for (size_t i = 0; i < size; i++) {
            Ciphertext ct;
            evaluator->multiply(cipher_data[i], lct.cipher_data[i], ct);
            ret.cipher_data[i] = ct;
        }
    }
    else {
        char buf[100];
        sprintf(buf, "In multiply: Length of BFVLongCiphertext(%ld) and BFVLongCiphertext(%ld) mismatch", len, lct.len);
        throw bfv_lenth_error(buf);
    }
    return ret;
}

void BFVLongCiphertext::send(NetIO* io, BFVLongCiphertext* lct) {
    assert(lct->len > 0);
    io->send_data(&(lct->len), sizeof(size_t));
    size_t size = lct->cipher_data.size();
    io->send_data(&size, sizeof(size_t));

    vector<uint64_t> ct_sizes(size);
    vector<string> ct_sers(size);
#pragma omp parallel for if (size > MAX_SZ)
    for (size_t ct = 0; ct < size; ct++) {
        std::stringstream os;
        uint64_t ct_size;
        lct->cipher_data[ct].save(os);
        ct_sizes[ct] = os.tellp();
        ct_sers[ct]  = os.str();
    }

    io->send_data(ct_sizes.data(), sizeof(uint64_t) * size);
    for (size_t i = 0; i < size; i++) {
        io->send_data(ct_sers[i].c_str(), ct_sizes[i]);
    }
    ct_sers.clear();
    ct_sers.shrink_to_fit();
}

void BFVLongCiphertext::recv(NetIO* io, BFVLongCiphertext* lct, SEALContext* context) {
    io->recv_data(&(lct->len), sizeof(size_t));
    size_t size;
    io->recv_data(&size, sizeof(size_t));
    lct->cipher_data.resize(size);

    vector<uint64_t> ct_sizes(size);
    char** ct_sers = new char*[size];
    io->recv_data(ct_sizes.data(), sizeof(uint64_t) * size);
    for (size_t i = 0; i < size; i++) {
        ct_sers[i] = new char[ct_sizes[i]];
        io->recv_data(ct_sers[i], ct_sizes[i]);
    }

#pragma omp parallel for if (size > MAX_SZ)
    for (size_t ct = 0; ct < size; ct++) {
        Ciphertext cct;
        std::stringstream is;
        is.write(ct_sers[ct], ct_sizes[ct]);
        cct.unsafe_load(*context, is);
        lct->cipher_data[ct] = cct;
        delete[] ct_sers[ct];
    }
    delete[] ct_sers;
}

// void BFVLongCiphertext::recv(NetIO *io, BFVLongCiphertext *lct, SEALContext *context) {
//     io->recv_data(&(lct->len), sizeof(size_t));
//     size_t size;
//     io->recv_data(&size, sizeof(size_t));

//     vector<uint64_t> ct_sizes(size);
//     io->recv_data(ct_sizes.data(), sizeof(uint64_t) * size);

//     for (size_t ct = 0; ct < size; ct++) {
//         Ciphertext cct;
//         std::stringstream is;
//         char *c_enc_result = new char[ct_sizes[ct]];
//         io->recv_data(c_enc_result, ct_sizes[ct]);
//         is.write(c_enc_result, ct_sizes[ct]);
//         cct.unsafe_load(*context, is);
//         lct->cipher_data.push_back(cct);
//         delete[] c_enc_result;
//     }
// }

void set_poly_coeffs_uniform(uint64_t* poly, uint32_t bitlen, std::shared_ptr<UniformRandomGenerator> random,
                             std::shared_ptr<const SEALContext::ContextData>& context_data) {
    assert(bitlen < 128 && bitlen > 0);
    auto& parms            = context_data->parms();
    auto& coeff_modulus    = parms.coeff_modulus();
    size_t coeff_count     = parms.poly_modulus_degree();
    size_t coeff_mod_count = coeff_modulus.size();
    uint64_t bitlen_mask   = (1ULL << (bitlen % 64)) - 1;

    RandomToStandardAdapter engine(random);
    for (size_t i = 0; i < coeff_count; i++) {
        if (bitlen < 64) {
            uint64_t noise = (uint64_t(engine()) << 32) | engine();
            noise &= bitlen_mask;
            for (size_t j = 0; j < coeff_mod_count; j++) {
                poly[i + (j * coeff_count)] = seal::util::barrett_reduce_64(noise, coeff_modulus[j]);
            }
        }
        else {
            uint64_t noise[2];  // LSB || MSB
            for (int j = 0; j < 2; j++) {
                noise[0] = (uint64_t(engine()) << 32) | engine();
                noise[1] = (uint64_t(engine()) << 32) | engine();
            }
            noise[1] &= bitlen_mask;
            for (size_t j = 0; j < coeff_mod_count; j++) {
                poly[i + (j * coeff_count)] = seal::util::barrett_reduce_128(noise, coeff_modulus[j]);
            }
        }
    }
}

void flood_ciphertext(Ciphertext& ct, std::shared_ptr<const SEALContext::ContextData>& context_data, uint32_t noise_len,
                      MemoryPoolHandle pool) {
    auto& parms            = context_data->parms();
    auto& coeff_modulus    = parms.coeff_modulus();
    size_t coeff_count     = parms.poly_modulus_degree();
    size_t coeff_mod_count = coeff_modulus.size();

    auto noise(seal::util::allocate_poly(coeff_count, coeff_mod_count, pool));
    std::shared_ptr<UniformRandomGenerator> random(parms.random_generator()->create());

    set_poly_coeffs_uniform(noise.get(), noise_len, random, context_data);
    for (size_t i = 0; i < coeff_mod_count; i++) {
        add_poly_poly_coeffmod(noise.get() + (i * coeff_count), ct.data() + (i * coeff_count), coeff_count,
                               coeff_modulus[i], ct.data() + (i * coeff_count));
    }

    set_poly_coeffs_uniform(noise.get(), noise_len, random, context_data);
    for (size_t i = 0; i < coeff_mod_count; i++) {
        add_poly_poly_coeffmod(noise.get() + (i * coeff_count), ct.data(1) + (i * coeff_count), coeff_count,
                               coeff_modulus[i], ct.data(1) + (i * coeff_count));
    }
}