#include <protocols/linear.h>
using namespace PrivLR;

int main(int argc, const char **argv) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(-1, 1);
    std::uniform_real_distribution<> positive_dist(0, 1);

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

    size_t size = 10;
    size_t data_size = 20;
    vector<vector<double>> in_a(data_size, vector<double>(size));
    vector<double> in_b(size);
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < data_size; j++) {
            in_a[j][i] = dist(gen);
        }
        in_b[i] = dist(gen);
    }
    Linear *linear = new Linear(party_, io_pack);
    auto res = linear->dot_product(in_a, in_b);

    if (party_ == ALICE) {
        for (size_t i = 0; i < data_size; i++) {
            io_pack->send_data(in_a[i].data(), sizeof(double) * size);
        }
        io_pack->send_data(in_b.data(), sizeof(double) * size);
        io_pack->send_data(res.data(), sizeof(double) * data_size);
    } else {
        vector<double> res_a(data_size), true_res(data_size);
        vector<vector<double>> in_a_a(data_size, vector<double>(size));
        vector<double> in_b_a(size);
        for (size_t i = 0; i < data_size; i++) {
            io_pack->recv_data(in_a_a[i].data(), sizeof(double) * size);
        }
        io_pack->recv_data(in_b_a.data(), sizeof(double) * size);
        io_pack->recv_data(res_a.data(), sizeof(double) * data_size);
        for (size_t j = 0; j < data_size; j++) {
            for (size_t i = 0; i < size; i++) {
                true_res[j] += (in_a[j][i] + in_a_a[j][i]) * (in_b[i] + in_b_a[i]);
                // std::cout << in_a_b[j][i] * in_b[i] << "\n";
            }
        }
        std::cout << "error: " << "\n";
        for (size_t j = 0; j < data_size; j++) {
            cout  << true_res[j] - res[j] - res_a[j] << "\n";
        }
    }
}