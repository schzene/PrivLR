#include <cassert>
#include <party.h>
#include <random>
#include <utils.h>
#include <vector>
using namespace paillier;

ZZ stringToNumber(string str) {
    ZZ number = NTL::conv<ZZ>(str[0]);
    long len = str.length();
    for (long i = 1; i < len; i++) {
        number *= 128;
        number += NTL::conv<ZZ>(str[i]);
    }
    return number;
}

string numberToString(ZZ num) {
    long len = ceil(NTL::log(num) / log(128));
    string str;
    str.resize(len);
    for (long i = len - 1; i >= 0; i--) {
        str[i] = NTL::conv<int>(num % 128);
        num /= 128;
    }
    return str;
}

int main(int argc, const char **argv) {
    PublicKey pk;
    PrivateKey sk;
    keygen(pk, sk);
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // std::uniform_real_distribution<> dist(0, 1);
    // size_t size = 10;

    // vector<double> a(size, dist(gen));
    // vector<double> b(size, dist(gen));
    // vector<ZZ> a_fixed(size), b_fixed(size);
    // for (size_t i = 0; i < size; i++) {
    //     a_fixed[i] = uint64_t(a[i] * (1ull << 31));
    //     b_fixed[i] = uint64_t(b[i] * (1ull << 31));
    // }

    // CipherVector a_secret(a_fixed, pk), b_secret(b_fixed, pk);
    // a_secret.add_inplace(b_secret);
    // auto res_fixed = a_secret.decrypt(sk);
    // vector<double> res(size);
    // for (size_t i = 0; i < size; i++) {
    //     res[i] = NTL::to_double(res_fixed[i]) / (1ull << 31) - a[i] - b[i];
    //     // std::cout << res[i] << "\n";
    // }
    // ZZ test = a_secret.data[0];
    // std::cout << test << "\n\n";
    // std::cout << stringToNumber(numberToString(test)) << "\n";

    int party_ = argv[1][0] - '0';
    if (party_ == ALICE) {
        std::cout << "Party: ALICE"
                  << "\n";
    } else {
        party_ = BOB;
        std::cout << "Party: BOB"
                  << "\n";
    }
    IOPack *io_pack = new IOPack(party_);
    if (party_ == ALICE) {
        vector<ZZ> a(10, ZZ(10));
        CipherVector ct_a(a, pk);
        CipherVector::send(io_pack, &ct_a);
        CipherVector ct_aa;
        CipherVector::recv(io_pack, &ct_aa);
        vector<ZZ> aa = ct_aa.decrypt(sk);
        for (size_t i = 0; i < 10; i++) {
            std::cout << a[i] - aa[i] << "\n";
        }
    } else {
        CipherVector ct_a;
        CipherVector::recv(io_pack, &ct_a);
        vector<ZZ> a(10, ZZ(5));
        ct_a.add_plain_inplace(a);
        CipherVector::send(io_pack, &ct_a);
    }
}