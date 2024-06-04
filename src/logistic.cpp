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
            printf("Cycle remain: %3d", max_cycles);
            vector<double> classified = linear->dot_product(datas, weight);
            vector<double> h = non_linear->sigmoid(classified);
            vector<double> error(data_size);
#ifdef _OPENMP
#pragma omp parallel for
#endif
            // START_TIMER
            for (int i = 0; i < data_size; i++) {
                double dist = label[i] - h[i];
                if (abs(dist) < 1e-10) {
                    dist = 0;
                }
                error[i] = dist;
            }

            // Not a protocol content, only for statistical purposes
            {
                double sum_error = 0;
                vector<double> error_remote(data_size);
                io_pack->send_data(error.data(), sizeof(double) * data_size, false);
                io_pack->recv_data(error_remote.data(), sizeof(double) * data_size, false);
                for (size_t i = 0; i < data_size; i++) {
                    sum_error += abs(error_remote[i] + error[i]);
                }
                printf(", Error: %8.6lf\n", sum_error / data_size);
            }
            vector<double> delta_weight = linear->dot_product(datas, error, true);
            for (int i = 0; i < size; i++) {
                weight[i] += alpha * delta_weight[i];
            }
            max_cycles--;
        }
    }

    void Logistic::stocGradAscent(vector<vector<double>> &datas, vector<int> &label,
                                  int numIter, double alpha) {
        double h = 0.0;
        int i = 0;
        int j = 0;
        double error = 0.0;
        const size_t data_size = datas.size(), size = datas[0].size();
        vector<int> rand_idx(data_size);
        weight = vector<double>(size, 1);
        for (i = 0; i < data_size; i++) {
            rand_idx[i] = i;
        }

        for (int k = 0; k < numIter; k++) {
            std::random_shuffle(rand_idx.begin(), rand_idx.end());

            for (i = 0; i < data_size; i++) {
                alpha = 4 / (1 + k + i) + alpha;
                h = non_linear->sigmoid(linear->dot_product(datas[rand_idx[i]], weight));
                error = label[rand_idx[i]] - h;
                for (j = 0; j < weight.size(); j++) {
                    weight[j] += alpha * error * datas[rand_idx[i]][j];
                }
            }
        }
    }
}