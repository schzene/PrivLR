#include <logistic-paillier.h>

using namespace PrivLR_Paillier;
int party = 1, num_iter = 25, dataset = 1;

INIT_TIMER

void load_dataset_base(vector<vector<double>>& datas, vector<int>& label, const string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        exit(1);
    }
    string line;
    while (std::getline(file, line)) {
        std::istringstream record(line);
        vector<double> data;
        data.push_back(1);
        double temp;
        while (record >> temp) {
            data.push_back(temp);
        }
        label.push_back(int(temp));
        data.pop_back();
        datas.push_back(data);
    }
}

double scalarProduct(const vector<double>& w, const vector<double>& x) {
    double ret = 0.0;
    for (int i = 0; i < w.size(); i++) {
        ret += w[i] * x[i];
    }
    return ret;
}

inline double sigmoid(const double& z) {
    return 1 / (1 + exp(-z));
}

vector<vector<double>> matTranspose(vector<vector<double>>& dataMat) {
    vector<vector<double>> ret(dataMat[0].size(), vector<double>(dataMat.size(), 0));
    for (int i = 0; i < ret.size(); i++)
        for (int j = 0; j < ret[0].size(); j++)
            ret[i][j] = dataMat[j][i];
    return ret;
}

void gradAscent(vector<double>& weight, vector<vector<double>>& dataMat, vector<int>& labelMat, int maxCycles = 1000,
                double alpha = 0.01) {
    const size_t data_size          = dataMat.size();
    vector<vector<double>> dataMatT = matTranspose(dataMat);
    while (maxCycles > 0) {
        vector<double> h;
        vector<double> error;
        double sum_err = 0;
        for (auto& data : dataMat)
            h.push_back(sigmoid(scalarProduct(data, weight)));
        for (int i = 0; i < labelMat.size(); i++) {
            double dist = labelMat[i] - h[i];
            if (abs(dist) < 1e-10)
                dist = 0;
            error.push_back(dist);
        }
        for (int i = 0; i < weight.size(); i++)
            weight[i] += alpha * scalarProduct(dataMatT[i], error);
        double sum_error = 0.;
        for (int i = 0; i < data_size; ++i) {
            sum_error += -1 * labelMat[i] * log(h[i]) - (1 - labelMat[i]) * log(1 - h[i]);
        }
        printf("loss: %.10lf\n", sum_error / data_size);
        maxCycles--;
    }
}

inline int classify(vector<double>& data, vector<double>& weights) {
    return sigmoid(scalarProduct(data, weights)) > 0.5 ? 1 : 0;
}

double testResult(vector<vector<double>>& testDataMat, vector<int>& testDataLabel, vector<double>& weight) {
    double errCount = 0.0;
    double dataSize = testDataMat.size();
    for (int i = 0; i < dataSize; i++)
        if (classify(testDataMat[i], weight) != testDataLabel[i])
            errCount += 1.0;
    return errCount / dataSize;
}

void testPlaintext() {
    INIT_TIMER
    std::cout << "**************************************************\n"
              << "testPlaintext(Base)\n"
              << "**************************************************\n";
    vector<vector<double>> base_train_mat;
    vector<int> base_train_label;
    string base_train_file;
    if (dataset == 1) {
        base_train_file = "/data/PrivLR/ACAD_train";
    } else if (dataset == 2) {
        base_train_file = "/data/PrivLR/HFCR_train";
    } else {
        base_train_file = "/data/PrivLR/WIBC_train";
    }
    load_dataset_base(base_train_mat, base_train_label, base_train_file);

    vector<vector<double>> base_test_mat;
    vector<int> base_test_label;
    string base_test_file;
    if (dataset == 1) {
        base_test_file = "/data/PrivLR/ACAD_test";
    } else if (dataset == 2) {
        base_test_file = "/data/PrivLR/HFCR_test";
    } else {
        base_test_file = "/data/PrivLR/WIBC_test";
    }
    load_dataset_base(base_test_mat, base_test_label, base_test_file);

    vector<double> base_weight(base_train_mat[0].size(), 1);

    START_TIMER
    gradAscent(base_weight, base_train_mat, base_train_label, num_iter, 0.008);
    STOP_TIMER("base")
    auto err = testResult(base_test_mat, base_test_label, base_weight);
    std::cout << "accuracy: " << (1 - err) * 100 << " %\n";
    std::cout << "**************************************************\n\n";
}

void PrivLR_test(int& _party) {
    std::cout << "**************************************************\n"
              << "PrivLR_test:\n"
              << "**************************************************\n";
    string train_file, test_file;
    if (_party == ALICE) {
        std::cout << "Party: ALICE"
                  << "\n";
        if (dataset == 1) {
            train_file = "/data/PrivLR/ACAD_alice";
            test_file = "/data/PrivLR/ACAD_alice_test";
        } else if (dataset == 2) {
            train_file = "/data/PrivLR/HFCR_alice";
            test_file = "/data/PrivLR/HFCR_alice_test";
        } else {
            train_file = "/data/PrivLR/WIBC_alice";
            test_file = "/data/PrivLR/WIBC_alice_test";
        }
    }
    else {
        _party = BOB;
        std::cout << "Party: BOB"
                  << "\n";
        if (dataset == 1) {
            train_file = "/data/PrivLR/ACAD_bob";
            test_file = "/data/PrivLR/ACAD_bob_test";
        } else if (dataset == 2) {
            train_file = "/data/PrivLR/HFCR_bob";
            test_file = "/data/PrivLR/HFCR_bob_test";
        } else {
            train_file = "/data/PrivLR/WIBC_bob";
            test_file = "/data/PrivLR/WIBC_bob_test";
        }
    }
    IOPack* io_pack = new IOPack(_party);

    vector<vector<double>> train_mat, test_mat;
    vector<int> train_label, test_label;
    load_dataset(train_mat, train_label, train_file);
    load_dataset(test_mat, test_label, test_file);
    const size_t size = train_mat[0].size();

    Logistic* logistic = new Logistic(_party, io_pack);
    START_TIMER
    logistic->gradAscent(train_mat, train_label, num_iter, 0.008);
    STOP_TIMER("train")

    size_t comm   = io_pack->get_comm();
    size_t rounds = io_pack->get_rounds();
    if (comm < 1024) {
        printf("data size of communication: %ld B\n", comm);
    }
    else if (comm < 1024 * 1024) {
        printf("data size of communication: %.2lf KB\n", comm / 1024.);
    }
    else if (comm < 1024 * 1024 * 1024) {
        printf("data size of communication: %.2lf MB\n", comm / (1024. * 1024.));
    }
    else {
        printf("data size of communication: %.2lf MB\n", comm / (1024. * 1024. * 1024.));
    }
    std::cout << "rounds of communication: " << rounds << "\n";

    size_t test_size      = test_mat.size();
    vector<double> result = logistic->classify(test_mat);
    if (_party == BOB) {
        io_pack->send_data(result.data(), sizeof(double) * test_size);
        io_pack->send_data(test_label.data(), sizeof(int) * test_size);
    }
    else {
        vector<double> result_remote(test_size);
        vector<int> test_label_remote(test_size);
        io_pack->recv_data(result_remote.data(), sizeof(double) * test_size);
        io_pack->recv_data(test_label_remote.data(), sizeof(int) * test_size);
        double acc = 0.;
        for (size_t i = 0; i < test_size; i++) {
            int res = result[i] + result_remote[i] > 0.5 ? 1 : 0;
            if (res == (test_label[i] + test_label_remote[i])) {
                acc += 1.;
            }
        }
        acc /= test_size;
        std::cout << "accuracy: " << acc * 100 << " %\n";
    }

    delete io_pack;
    delete logistic;
    std::cout << "**************************************************\n\n";
}

int main(int argc, const char** argv) {
    assert(argc == 4);
    party = atoi(argv[1]);
    dataset = atoi(argv[2]);
    num_iter = atoi(argv[3]);
    PrivLR_test(party);
    if (party == BOB) {
        testPlaintext();
    }
}