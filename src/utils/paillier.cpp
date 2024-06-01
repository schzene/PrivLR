#include "paillier.h"

namespace paillier {
    ZZ string2ZZ(const string &str) {
        ZZ num = NTL::conv<ZZ>(str[0]);
        long len = str.length();
        for (long i = 1; i < len; i++) {
            num *= 128;
            num += NTL::conv<ZZ>(str[i]);
        }
        return num;
    }

    string ZZ2string(ZZ num) {
        long len = ceil(NTL::log(num) / log(128));
        string str;
        str.resize(len);
        for (long i = len - 1; i >= 0; i--) {
            str[i] = NTL::conv<int>(num % 128);
            num /= 128;
        }
        return str;
    }

    void keygen(PublicKey &public_key, PrivateKey &private_key, size_t key_length) {
        ZZ p, q, phi, n, lambda, mu;
        do {
            p = NTL::GenPrime_ZZ(key_length);
            q = NTL::GenPrime_ZZ(key_length);
            while (p == q) {
                q = NTL::GenPrime_ZZ(key_length);
            }
            n = p * q;
            phi = (p - 1) * (q - 1);

        } while (NTL::GCD(n, phi) != 1);

        lambda = phi / NTL::GCD(p - 1, q - 1);
        mu = NTL::InvMod(lambda, n);
        public_key = PublicKey(n);
        private_key = PrivateKey(lambda, mu, public_key);
    }

    CipherVector::CipherVector(const vector<ZZ> &messages, const PublicKey &_public_key) : public_key(_public_key) {
        const size_t size = messages.size();
        data.resize(size);
        const ZZ n = public_key.n;
        const ZZ n_square = n * n;
        const ZZ second_part = NTL::PowerMod(NTL::RandomBnd(n), n, n_square);
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (size_t i = 0; i < size; i++) {
            data[i] = ((n * messages[i] + 1) * second_part) % (n_square);
        }
    }

    vector<ZZ> CipherVector::decrypt(const PrivateKey &private_key) const {
        const ZZ n = private_key.public_key.n;
        const ZZ n_square = n * n;
        auto L = [&n](const ZZ &x) { return (x - 1) / n; };
        const size_t size = data.size();
        vector<ZZ> messages(size);
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (size_t i = 0; i < size; i++) {
            messages[i] = (L(NTL::PowerMod(data[i], private_key.lambda, n_square)) * private_key.mu) % n;
        }
        return messages;
    }

    CipherVector CipherVector::add(const CipherVector &other) const {
        CipherVector result(public_key);
        const size_t this_size = data.size();
        const size_t other_size = other.data.size();
        const ZZ n = public_key.n;
        const ZZ n_square = n * n;
        if (this_size == 1) {
            result.data.resize(other_size);
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < other_size; i++) {
                result.data[i] = (data[0] * other.data[i]) % n_square;
            }
        } else if (other_size == 1) {
            result.data.resize(this_size);
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < this_size; i++) {
                result.data[i] = (data[i] * other.data[0]) % n_square;
            }
        } else if (this_size == other_size) {
            result.data.resize(this_size);
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < this_size; i++) {
                result.data[i] = (data[i] * other.data[i]) % n_square;
            }
        } else {
            char buf[100];
            sprintf(buf, "Length of CipherVector(%ld) and CipherVector(%ld) mismatch", this_size, other_size);
            throw lenth_error(buf);
        }
        return result;
    }

    void CipherVector::add_inplace(const CipherVector &other) {
        const size_t this_size = data.size();
        const size_t other_size = other.data.size();
        const ZZ n = public_key.n;
        const ZZ n_square = n * n;
        if (this_size == 1) {
            ZZ ct = data[0];
            data.resize(other_size);
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < other_size; i++) {
                data[i] = (ct * other.data[i]) % n_square;
            }
        } else if (other_size == 1) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < this_size; i++) {
                data[i] = (data[i] * other.data[0]) % n_square;
            }
        } else if (this_size == other_size) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < this_size; i++) {
                data[i] = (data[i] * other.data[i]) % n_square;
            }
        } else {
            char buf[100];
            sprintf(buf, "Length of CipherVector(%ld) and CipherVector(%ld) mismatch", this_size, other_size);
            throw lenth_error(buf);
        }
    }

    CipherVector CipherVector::add_plain(const vector<ZZ> &other) const {
        CipherVector result(public_key);
        const size_t this_size = data.size();
        const size_t other_size = other.size();
        const ZZ n = public_key.n;
        const ZZ n_square = n * n;
        ZZ r = NTL::RandomBnd(n);
        if (this_size == 1) {
            result.data.resize(other_size);
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < other_size; i++) {
                ZZ ct = ((n * other[i] + 1) * NTL::PowerMod(r, n, n_square)) % n_square;
                result.data[i] = (data[0] * ct) % n_square;
            }
        } else if (other_size == 1) {
            result.data.resize(this_size);
            ZZ ct = ((n * other[0] + 1) * NTL::PowerMod(r, n, n_square)) % n_square;
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < this_size; i++) {
                result.data[i] = (data[i] * ct) % n_square;
            }
        } else if (this_size == other_size) {
            result.data.resize(this_size);
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < other_size; i++) {
                ZZ ct = ((n * other[i] + 1) * NTL::PowerMod(r, n, n_square)) % n_square;
                result.data[i] = (data[i] * ct) % n_square;
            }
        } else {
            char buf[100];
            sprintf(buf, "Length of CipherVector(%ld) and message(%ld) mismatch", this_size, other_size);
            throw lenth_error(buf);
        }
        return result;
    }

    void CipherVector::add_plain_inplace(const vector<ZZ> &other) {
        const size_t this_size = data.size();
        const size_t other_size = other.size();
        const ZZ n = public_key.n;
        const ZZ n_square = n * n;
        ZZ r = NTL::RandomBnd(n);
        if (this_size == 1) {
            ZZ this_data = data[0];
            data.resize(other_size);
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < other_size; i++) {
                ZZ ct = ((n * other[i] + 1) * NTL::PowerMod(r, n, n_square)) % n_square;
                data[i] = (this_data * ct) % n_square;
            }
        } else if (other_size == 1) {
            ZZ ct = ((n * other[0] + 1) * NTL::PowerMod(r, n, n_square)) % n_square;
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < this_size; i++) {
                data[i] = (data[i] * ct) % n_square;
            }
        } else if (this_size == other_size) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < other_size; i++) {
                ZZ ct = ((n * other[i] + 1) * NTL::PowerMod(r, n, n_square)) % n_square;
                data[i] = (data[i] * ct) % n_square;
            }
        } else {
            char buf[100];
            sprintf(buf, "Length of CipherVector(%ld) and message(%ld) mismatch", this_size, other_size);
            throw lenth_error(buf);
        }
    }

    CipherVector CipherVector::mul_plain(const vector<ZZ> &other) const {
        CipherVector result(public_key);
        const size_t this_size = data.size();
        const size_t other_size = other.size();
        const ZZ n = public_key.n;
        const ZZ n_square = n * n;
        if (this_size == 1) {
            result.data.resize(other_size);
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < other_size; i++) {
                result.data[i] = NTL::PowerMod(data[0], other[i], n_square) % n_square;
            }
        } else if (other_size == 1) {
            result.data.resize(this_size);
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < other_size; i++) {
                result.data[i] = NTL::PowerMod(data[i], other[0], n_square) % n_square;
            }
        } else if (this_size == other_size) {
            result.data.resize(this_size);
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < other_size; i++) {
                result.data[i] = NTL::PowerMod(data[i], other[i], n_square) % n_square;
            }
        } else {
            char buf[100];
            sprintf(buf, "Length of CipherVector(%ld) and message(%ld) mismatch", this_size, other_size);
            throw lenth_error(buf);
        }
        return result;
    }

    void CipherVector::mul_plain_inplace(const vector<ZZ> &other) {
        const size_t this_size = data.size();
        const size_t other_size = other.size();
        const ZZ n = public_key.n;
        const ZZ n_square = n * n;
        if (this_size == 1) {
            ZZ this_data = data[0];
            data.resize(other_size);
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < other_size; i++) {
                data[i] = NTL::PowerMod(this_data, other[i], n_square) % n_square;
            }
        } else if (other_size == 1) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < other_size; i++) {
                data[i] = NTL::PowerMod(data[i], other[0], n_square) % n_square;
            }
        } else if (this_size == other_size) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (size_t i = 0; i < other_size; i++) {
                data[i] = NTL::PowerMod(data[i], other[i], n_square) % n_square;
            }
        } else {
            char buf[100];
            sprintf(buf, "Length of CipherVector(%ld) and message(%ld) mismatch", this_size, other_size);
            throw lenth_error(buf);
        }
    }

    void CipherVector::send(IOPack *io_pack, CipherVector *cv) {
        string pk_str = ZZ2string(cv->public_key.n);
        size_t pk_str_size = pk_str.size();
        io_pack->send_data(&pk_str_size, sizeof(size_t));
        io_pack->send_data(pk_str.data(), sizeof(char) * pk_str_size);
        size_t data_size = cv->data.size();
        io_pack->send_data(&data_size, sizeof(size_t));
        for (size_t i = 0; i < data_size; i++) {
            string data_str = ZZ2string(cv->data[i]);
            size_t data_str_size = data_str.size();
            io_pack->send_data(&data_str_size, sizeof(size_t));
            io_pack->send_data(data_str.data(), sizeof(char) * data_str_size);
        }
    }

    void CipherVector::recv(IOPack *io_pack, CipherVector *cv) {
        size_t pk_str_size;
        string pk_str;
        io_pack->recv_data(&pk_str_size, sizeof(size_t));
        pk_str.resize(pk_str_size);
        io_pack->recv_data(pk_str.data(), sizeof(char) * pk_str_size);
        cv->public_key.n = string2ZZ(pk_str);
        size_t data_size;
        io_pack->recv_data(&data_size, sizeof(size_t));
        cv->data.resize(data_size);
        for (size_t i = 0; i < data_size; i++) {
            size_t data_str_size;
            string data_str;
            io_pack->recv_data(&data_str_size, sizeof(size_t));
            data_str.resize(data_str_size);
            io_pack->recv_data(data_str.data(), sizeof(char) * data_str_size);
            cv->data[i] = string2ZZ(data_str);
        }
    }
}
