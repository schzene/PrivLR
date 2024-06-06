#include <logistic.h>

using namespace PrivLR;

INIT_TIMER

int main(int argc, const char **argv) {
    string train_file, test_file;
    int party_ = argv[1][0] - '0';
    if (party_ == ALICE) {
        std::cout << "Party: ALICE"
                  << "\n";
        train_file = "/data/PrivLR/WIBC_alice";
        test_file = "/data/PrivLR/WIBC_alice_test";
    } else {
        party_ = BOB;
        std::cout << "Party: BOB"
                  << "\n";
        train_file = "/data/PrivLR/WIBC_bob";
        test_file = "/data/PrivLR/WIBC_bob_test";
    }
    IOPack *io_pack = new IOPack(party_);

    vector<vector<double>> train_mat, test_mat;
    vector<int> train_label, test_label;
    load_dataset(train_mat, train_label, train_file);
    load_dataset(test_mat, test_label, test_file);
    const size_t size = train_mat[0].size();

    Logistic *logistic = new Logistic(party_, io_pack);
    START_TIMER
    logistic->gradAscent(train_mat, train_label, 10, 0.008);
    STOP_TIMER("train")

    size_t comm = io_pack->get_comm();
    size_t rounds = io_pack->get_rounds();
    if (comm < 1024) {
        printf("data size of communication: %ld B\n", comm);
    } else if (comm < 1024 * 1024) {
        printf("data size of communication: %.2lf KB\n", comm / 1024.);
    } else if (comm < 1024 * 1024 * 1024) {
        printf("data size of communication: %.2lf MB\n", comm / (1024. * 1024.));
    } else {
        printf("data size of communication: %.2lf MB\n", comm / (1024. * 1024. * 1024.));
    }
    std::cout << "rounds of communication: " << rounds << "\n";

    size_t test_size = test_mat.size();
    vector<double> result = logistic->classify(test_mat), result_remote(test_size);
    vector<int> test_label_remote(test_size);
    io_pack->send_data(result.data(), sizeof(double) * test_size);
    io_pack->recv_data(result_remote.data(), sizeof(double) * test_size);
    io_pack->send_data(test_label.data(), sizeof(int) * test_size);
    io_pack->recv_data(test_label_remote.data(), sizeof(int) * test_size);
    double acc = 0.;
    for (size_t i = 0; i < test_size; i++) {
        int res = result[i] + result_remote[i] > 0.5? 1: 0;
        if (res == (test_label[i] + test_label_remote[i])) {
            acc += 1.;
        }
    }
    acc /= test_size;
    std::cout << "accuracy: " << acc * 100 << " %\n";

    delete io_pack;
    delete logistic;
}