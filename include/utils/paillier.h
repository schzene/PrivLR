#ifndef PRIV_LR_PAILLIER_H__
#define PRIV_LR_PAILLIER_H__

#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <NTL/ZZ.h>
#include <NTL/ZZ_pXFactoring.h>

#include "io.h"

using NTL::ZZ;
using std::vector;

namespace paillier {
    ZZ string2ZZ(const string &str);
    string ZZ2string(ZZ num);

    class Ciphertext {
    public:
        ZZ data;

        Ciphertext() {}
        Ciphertext(const ZZ &_data) : data(_data) {}
    };

    class PublicKey {
    public:
        ZZ n;

        PublicKey() {}
        PublicKey(const ZZ &_n) : n(_n) {}
        PublicKey(const PublicKey &other) : n(other.n) {}
    };

    class PrivateKey {
    public:
        ZZ lambda;
        ZZ mu;
        PublicKey public_key;

        PrivateKey() {}
        PrivateKey(const ZZ &_lambda, const ZZ &_mu, const PublicKey &_public_key) : lambda(_lambda),
                                                                                     mu(_mu),
                                                                                     public_key(_public_key) {}
    };

    void keygen(PublicKey &public_key, PrivateKey &private_key, size_t key_length = 1024);

    inline Ciphertext encrypt(const ZZ &message, const PublicKey &public_key) {
        ZZ n = public_key.n, r = NTL::RandomBnd(n);
        return Ciphertext(((n * message + 1) * NTL::PowerMod(r, n, n * n)) % (n * n));
    }

    inline ZZ decrypt(const Ciphertext &ciphertext, const PrivateKey &private_key) {
        ZZ n = private_key.public_key.n;
        auto L = [&n](const ZZ &x) { return (x - 1) / n; };
        ZZ message = (L(NTL::PowerMod(ciphertext.data, private_key.lambda, n * n)) * private_key.mu) % n;
        return message;
    }

    inline Ciphertext add(const Ciphertext &ciphertext1, const Ciphertext &ciphertext2, const PublicKey &public_key) {
        return Ciphertext((ciphertext1.data * ciphertext2.data) % (public_key.n * public_key.n));
    }

    inline void add_inplace(Ciphertext &destination, const Ciphertext &source, const PublicKey &public_key) {
        destination.data = (destination.data * source.data) % (public_key.n * public_key.n);
    }

    inline Ciphertext add_plain(const Ciphertext &ciphertext, const ZZ &plaintext, const PublicKey &public_key) {
        Ciphertext ct = encrypt(plaintext, public_key);
        return Ciphertext((ciphertext.data * ct.data) % (public_key.n * public_key.n));
    }

    inline void add_plain_inplace(Ciphertext &destination, const ZZ &source, const PublicKey &public_key) {
        Ciphertext ct = encrypt(source, public_key);
        destination.data = (destination.data * ct.data) % (public_key.n * public_key.n);
    }

    inline Ciphertext mul_plain(const Ciphertext &ciphertext, const ZZ &plaintext, const PublicKey &public_key) {
        return Ciphertext(NTL::PowerMod(ciphertext.data, plaintext, public_key.n * public_key.n));
    }

    inline void mul_plain_inplace(Ciphertext &destination, const ZZ &source, const PublicKey &public_key) {
        destination.data = NTL::PowerMod(destination.data, source, public_key.n * public_key.n);
    }

    class lenth_error : public std::exception {
        const char *message;

    public:
        lenth_error(const char *msg) : message(msg) {}
        const char *what() const throw() override {
            return message;
        }
    };

    class CipherVector {
    public:
        vector<ZZ> data;
        PublicKey public_key;

        CipherVector() {}
        CipherVector(const PublicKey &_public_key) : public_key(_public_key) {}
        CipherVector(const vector<ZZ> &messages, const PublicKey &public_key);
        vector<ZZ> decrypt(const PrivateKey &private_key) const;
        CipherVector add(const CipherVector &other) const;
        void add_inplace(const CipherVector &other);
        CipherVector add_plain(const vector<ZZ> &other) const;
        void add_plain_inplace(const vector<ZZ> &other);
        CipherVector mul_plain(const vector<ZZ> &other) const;
        void mul_plain_inplace(const vector<ZZ> &other);

        static void send(IOPack *io_pack, CipherVector *cv);
        static void recv(IOPack *io_pack, CipherVector *cv);
    };
}

#endif