#include <logistic.h>

using namespace PrivLR;

int main(int argc, const char **argv) {
    vector<vector<double>> train_mat;
    vector<int> train_label;
    string train_file("traindata.txt");
    load_dataset(train_mat, train_label, train_file);

    vector<vector<double>> test_mat;
    vector<int> test_label;
    string testFile("traindata.txt");
    load_dataset(test_mat, test_label, testFile);

    // double err = testResult(test_mat, test_label, weight);

    // std::cout << "the error rate is: " << err << std::endl;
    IOPack *io_pack;
    if (argc > 1) {
        Logistic *logistic = new Logistic(argv[1][0] - '0', io_pack);
        logistic->gradAscent(train_mat, train_label);
        auto res = logistic->classify(test_mat);
        delete logistic;
    }
    delete io_pack;
}