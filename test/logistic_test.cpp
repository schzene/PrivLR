#include <logistic.h>

using namespace PrivLR;

int main(int argc, const char **argv)
{
    string train_file, test_file;
    int party_ = argv[1][0] - '0';
    if (party_ == ALICE)
    {
        std::cout << "Party: ALICE"
                  << "\n";
        train_file = "/root/PrivLR/traindata_alice.txt";
        test_file = "/root/PrivLR/traindata_alice.txt";
    }
    else
    {
        party_ = BOB;
        std::cout << "Party: BOB"
                  << "\n";
        train_file = "/root/PrivLR/traindata_bob.txt";
        test_file = "/root/PrivLR/traindata_bob.txt";
    }
    IOPack *io_pack = new IOPack(party_);

    vector<vector<double>> train_mat;
    vector<int> train_label;
    load_dataset(train_mat, train_label, train_file);

    vector<vector<double>> test_mat;
    vector<int> test_label;
    load_dataset(test_mat, test_label, test_file);

    // double err = testResult(test_mat, test_label, weight);

    // std::cout << "the error rate is: " << err << std::endl;
    Logistic *logistic = new Logistic(party_, io_pack);
    logistic->gradAscent(train_mat, train_label, 100, 0.001);
    std::cout << "weight = " << logistic->weight[0] << ", " << logistic->weight[1] << ", " << logistic->weight[2] << "\n";
    // vector<double> res = logistic->classify(test_mat);
    // vector<double> res_remote(res.size());
    // if (party_ == ALICE) {
    //     io_pack->recv_data(res_remote.data(), sizeof(double) * res.size());
    // } else {
    //     io_pack->send_data(res.data(), sizeof(double) * res.size());
    // }

    delete logistic;
    delete io_pack;
}