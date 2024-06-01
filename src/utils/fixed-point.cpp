/*
Authors: Deevashwer Rathee
Copyright:
Copyright (c) 2021 Microsoft Research
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "fixed-point.h"

FixArray FixArray::subset(int i, int j) {
    assert(i >= 0 && j <= size && i < j);
    int sz = j - i;
    FixArray ret(this->party, sz, this->signed_, this->ell, this->s);
    memcpy(ret.data, this->data + i, sz * sizeof(uint64_t));
    return ret;
}

FixArray concat(const vector<FixArray> &x) {
    int N = x.size();
    int sz = x[0].size;
    bool signed_ = x[0].signed_;
    int ell = x[0].ell;
    int s = x[0].s;
    int party = x[0].party;
    for (int i = 1; i < N; i++) {
        sz += x[i].size;
        assert(signed_ == x[i].signed_);
        assert(ell == x[i].ell);
        assert(s == x[i].s);
        assert(party == x[i].party);
    }
    FixArray ret(party, sz, signed_, ell, s);
    int offset = 0;
    for (int i = 0; i < N; i++) {
        int n = x[i].size;
        memcpy(ret.data + offset, x[i].data, n * sizeof(uint64_t));
        offset += n;
    }
    return ret;
}

template <class T>
std::vector<T> FixArray::get_native_type() {
    assert(this->party == PUBLIC);
    if constexpr (is_same_v<T, uint32_t> || is_same_v<T, uint64_t>) {
        assert(this->signed_ == false);
    }
    vector<T> ret(this->size);
    double den = pow(2.0, this->s);
    for (int i = 0; i < this->size; i++) {
        int64_t data_ = (this->signed_ ? signed_val(this->data[i], this->ell) : this->data[i]);
        ret[i] = T(data_ / den);
    }
    return ret;
}

std::ostream &operator<<(std::ostream &os, FixArray &other) {
    assert(other.party == PUBLIC);
    vector<double> dbl_other = other.get_native_type<double>();
    os << "(ell: " << other.ell << ", s: " << other.s << ") \t[";
    for (int i = 0; i < other.size; i++) {
        int64_t data_ =
            (other.signed_ ? signed_val(other.data[i], other.ell) : other.data[i]);
        std::string tmp_data = std::bitset<64>(data_).to_string();
        os << dbl_other[i] << " int=" << data_ << " ("
           << tmp_data.substr(64 - other.ell, 64) << ")\t";
    }
    os << "]";
    return os;
}

template vector<uint32_t> FixArray::get_native_type();
template vector<int32_t> FixArray::get_native_type();
template vector<uint64_t> FixArray::get_native_type();
template vector<int64_t> FixArray::get_native_type();
template vector<float> FixArray::get_native_type();
template vector<double> FixArray::get_native_type();

FixArray FixOp::input(int party_, int sz, uint64_t *data_, bool signed__, int ell_, int s_) {
    FixArray ret((party_ == PUBLIC ? party_ : this->party), sz, signed__, ell_, s_);
    uint64_t ell_mask_ = ret.ell_mask();
    if ((this->party == party_) || (party_ == PUBLIC)) {
        memcpy(ret.data, data_, sz * sizeof(uint64_t));
        for (int i = 0; i < sz; i++) {
            ret.data[i] &= ell_mask_;
        }
    } else {
        for (int i = 0; i < sz; i++) {
            ret.data[i] = 0;
        }
    }
    return ret;
}

FixArray FixOp::input(int party_, int sz, uint64_t data_, bool signed__, int ell_, int s_) {
    FixArray ret((party_ == PUBLIC ? party_ : this->party), sz, signed__, ell_, s_);
    uint64_t ell_mask_ = ret.ell_mask();
    if ((this->party == party_) || (party_ == PUBLIC)) {
        for (int i = 0; i < sz; i++) {
            ret.data[i] = data_ & ell_mask_;
        }
    } else {
        for (int i = 0; i < sz; i++) {
            ret.data[i] = 0;
        }
    }
    return ret;
}

FixArray FixOp::output(int party_, const FixArray &x) {
    if (x.party == PUBLIC) {
        return x;
    }
    int sz = x.size;
    int ret_party = (party_ == PUBLIC || party_ == x.party ? PUBLIC : x.party);
    FixArray ret(ret_party, sz, x.signed_, x.ell, x.s);
#ifdef _OPENMP
#pragma omp parallel num_threads(2)
    {
        if (omp_get_thread_num() == 1 && party_ != BOB) {
#endif
            if (party == ALICE) {
                iopack->io_rev->recv_data(ret.data, sz * sizeof(uint64_t));
            } else { // party == BOB
                iopack->io_rev->send_data(x.data, sz * sizeof(uint64_t));
            }
#ifdef _OPENMP
        } else if (omp_get_thread_num() == 0 && party_ != ALICE) {
#endif
            if (party == ALICE) {
                iopack->io->send_data(x.data, sz * sizeof(uint64_t));
            } else { // party == sci::BOB
                iopack->io->recv_data(ret.data, sz * sizeof(uint64_t));
            }
#ifdef _OPENMP
        }
    }
#endif
    uint64_t ell_mask_ = x.ell_mask();
    for (int i = 0; i < sz; i++) {
        ret.data[i] = (ret.data[i] + x.data[i]) & ell_mask_;
    }
    return ret;
}

FixArray FixOp::add(const FixArray &x, const FixArray &y) {
    assert(x.size == y.size);
    assert(x.signed_ == y.signed_);
    assert(x.ell == y.ell);
    assert(x.s == y.s);

    bool x_cond, y_cond;
    int party_;
    if (x.party == PUBLIC && y.party == PUBLIC) {
        x_cond = false;
        y_cond = false;
        party_ = PUBLIC;
    } else {
        x_cond = (x.party == PUBLIC) && (this->party == BOB);
        y_cond = (y.party == PUBLIC) && (this->party == BOB);
        party_ = this->party;
    }
    FixArray ret(party_, x.size, x.signed_, x.ell, x.s);
    uint64_t ell_mask_ = x.ell_mask();
    for (int i = 0; i < x.size; i++) {
        ret.data[i] = ((x_cond ? 0 : x.data[i]) + (y_cond ? 0 : y.data[i])) & ell_mask_;
    }
    return ret;
}

FixArray FixOp::add(const FixArray &x, uint64_t y) {
    FixArray y_fix = this->input(PUBLIC, x.size, y, x.signed_, x.ell, x.s);
    return this->add(x, y_fix);
}

FixArray FixOp::sub(const FixArray &x, const FixArray &y) {
    FixArray neg_y = this->mul(y, uint64_t(-1));
    return this->add(x, neg_y);
}

FixArray FixOp::sub(const FixArray &x, uint64_t y) {
    FixArray y_fix = this->input(PUBLIC, x.size, y, x.signed_, x.ell, x.s);
    return this->sub(x, y_fix);
}

FixArray FixOp::sub(uint64_t x, const FixArray &y) {
    FixArray x_fix = this->input(PUBLIC, y.size, x, y.signed_, y.ell, y.s);
    return this->sub(x_fix, y);
}

FixArray FixOp::mul(const FixArray &x, const FixArray &y, int ell,
                    uint8_t *msb_x, uint8_t *msb_y) {
    assert(x.party != PUBLIC || y.party != PUBLIC);
    assert(x.size == y.size);
    assert(x.signed_ || (x.signed_ == y.signed_));
    assert(ell >= x.ell && ell >= y.ell && ell <= x.ell + y.ell);
    assert(ell < 64);
    FixArray ret(this->party, x.size, x.signed_, ell, x.s + y.s);
    if (x.party == PUBLIC || y.party == PUBLIC) {
        FixArray x_ext = this->extend(x, ell, msb_x);
        FixArray y_ext = this->extend(y, ell, msb_y);
        uint64_t ret_mask = ret.ell_mask();
        for (int i = 0; i < x.size; i++) {
            ret.data[i] = (x_ext.data[i] * y_ext.data[i]) & ret_mask;
        }
    } else {
        mult->hadamard_product(x.size, x.data, y.data, ret.data, x.ell, y.ell, ell,
                               x.signed_, y.signed_, MultMode::None, msb_x, msb_y);
    }
    return ret;
}

FixArray FixOp::mul(const FixArray &x, uint64_t y, int ell, uint8_t *msb_x) {
    assert(ell >= x.ell);
    FixArray ret;
    if (ell > x.ell) {
        ret = this->extend(x, ell, msb_x);
    } else {
        ret = x;
    }
    uint64_t ell_mask_ = ret.ell_mask();
    for (int i = 0; i < x.size; i++) {
        ret.data[i] = (y * ret.data[i]) & ell_mask_;
    }
    return ret;
}

FixArray FixOp::scale_up(const FixArray &x, int ell, int s) {
    assert(ell - x.ell <= s - x.s);
    assert(s >= x.s);
    FixArray ret(x.party, x.size, x.signed_, ell, s);
    uint64_t ell_mask_ = ret.ell_mask();
    for (int i = 0; i < x.size; i++) {
        ret.data[i] = (x.data[i] << (s - x.s)) & ell_mask_;
    }
    return ret;
}

FixArray FixOp::reduce(const FixArray &x, int ell) {
    assert(ell <= x.ell && ell > 0);
    FixArray ret(x.party, x.size, x.signed_, ell, x.s);
    uint64_t ell_mask_ = ret.ell_mask();
    for (int i = 0; i < x.size; i++) {
        ret.data[i] = x.data[i] & ell_mask_;
    }
    return ret;
}

// A0 \in (1/4, 1)
inline uint64_t recp_lookup_c0(uint64_t index, int m) {
    uint64_t k = 1ULL << m;
    double p = 1 + (double(index) / double(k));
    double A1 = 1.0 / (p * (p + 1.0 / double(k)));
    int32_t scale = m + 3;
    uint64_t mask = (1ULL << scale) - 1;
    uint64_t val = uint64_t(A1 * (1ULL << scale)) & mask;
    return val;
}

// A1 \in (1/2, 1)
inline uint64_t recp_lookup_c1(uint64_t index, int m) {
    uint64_t k = 1ULL << m;
    double p = 1 + (double(index) / double(k));
    double z = (p * (p + (1.0 / double(k))));
    double A1 = ((1.0 / double(k * 2)) + sqrt(z)) / z;
    int32_t scale = 2 * m + 2;
    uint64_t mask = (1ULL << scale) - 1;
    uint64_t val = uint64_t(A1 * (1ULL << scale)) & mask;
    return val;
}