#include "logistic.h"

namespace PrivLR {
    Logistic::Logistic(int party, IOPack *io_pack) {
        assert(party == ALICE || party == BOB);
        assert(io_pack != nullptr);
        this->io_pack = io_pack;
        this->linear = new Linear(party, io_pack);
        this->non_linear = new NonLinear(party, io_pack);
    }

    Logistic::~Logistic() {
        if (linear != nullptr) {
            delete linear;
            linear = nullptr;
        }
        if (non_linear != nullptr) {
            delete non_linear;
            non_linear = nullptr;
        }
    }

    void Logistic::gradAscent(vector<vector<double>> &datas, vector<int> &label,
                              int max_cycles, const double alpha) {
        const size_t data_size = datas.size(), size = datas[0].size();
        weight = vector<double>(size, .5);
        while (max_cycles > 0) {
            vector<double> classified = linear->dot_product(datas, weight);
            vector<double> h = non_linear->sigmoid(classified);
            vector<double> error(data_size);
            // START_TIMER
            for (int i = 0; i < data_size; i++) {
                double dist = label[i] - h[i];
                if (abs(dist) < 1e-10) {
                    dist = 0;
                }
                error[i] = dist;
            }
            vector<double> delta_weight = linear->dot_product(datas, error, true);
            for (int i = 0; i < size; i++) {
                weight[i] += alpha * delta_weight[i];
            }
            printf("Cycle remain: %3d", max_cycles);
            max_cycles--;

            // Not a protocol content, only for statistical purposes
            {
                double sum_error = 0.;
                vector<double> h_remote(data_size);
                vector<int> label_remote(data_size);
                io_pack->send_data(h.data(), sizeof(double) * data_size, false);
                io_pack->send_data(label.data(), sizeof(int) * data_size, false);
                io_pack->recv_data(h_remote.data(), sizeof(double) * data_size, false);
                io_pack->recv_data(label_remote.data(), sizeof(int) * data_size, false);
                for (size_t i = 0; i < data_size; i++) {
                    sum_error += -1 * (label[i] + label_remote[i]) * log(h[i] + h_remote[i]) - (1 - label[i] - label_remote[i]) * log(1 - h[i] - h_remote[i]);
                }
                printf("loss: %.10lf\n", sum_error / data_size);
            }
        }
    }
}