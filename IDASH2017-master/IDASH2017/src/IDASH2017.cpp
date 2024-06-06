#include <NTL/BasicThreadPool.h>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include "GD.h"
#include "TestGD.h"

using namespace std;
using namespace NTL;

/* root@vultr:/home/john/Music/C16/IDASH2017-master/IDASH2017/Debug# ./MyIDASH2017 "../data/data103x1579.txt" 1 7 5 1 -1 1 5 1
 *
 * run: ./IDASH2017 trainfile isYfirst numIter k gammaUp gammaDown isInitZero fold isEncrypted testfile
 * ./HEML string bool long long double double bool long bool string
 * example: ./IDASH2017 "../data/data103x1579.txt" 1 7 5 1 -1 1 5 1
 * example: ./IDASH2017 "../data/1_training_data_csv" 1 7 5 1 -1 1 0 1 "../data/1_testing_data_csv"
 *
 * parameters:
 * trainfile - path to train file
 * isYfirst - {0,1} y parameter first OR last
 * numIter - number of iterations
 * kdeg - degree of sigmoid approximation function k in {3,5,7}
 * gammaUp - corresponds to learning rate
 * gammaDown - corresponds to learning rate
 * isInitZero - is initial weights zero or average
 * fold - folding method if arguments <= 8 we use folding method
 * isEncrypted - encrypted or plain
 * testfile - path to test file (checks if number of arguments > 8 then we use standard method
 *
 * current files that in data folder (filename isYfirst):
 * "../data/data5x500.txt" false
 * "../data/data9x1253.txt" false
 * "../data/data15x1500.txt" false
 * "../data/data16x101.txt" false
 * "../data/data27x148.txt" false
 * "../data/data51x653.txt" false
 * "../data/data67x216.txt" false
 * "../data/data103x1579.txt" true
 * "../data/1_training_data.csv" true
 * "../data/1_testing_data.csv" true
 *
 * FYI: approx 3 suggested iter: 4, 9, 18, 36, ...
 * FYI: approx 5 suggested iter: 3, 7, 14, 28, ...
 * FYI: approx 7 suggested iter: 3, 7, 14, 28, ...
 */


int main(int argc, char **argv) {
	//	size_t currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	//	size_t peakAfterSchemeSize = getPeakRSS() >> 20;
	//	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	//	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;
	SetNumThreads(36);

	//./MyIDASH2017 \"../data/data103x1579.txt\" 1 7 5 1 -1 1 5 1
	//string trainfile("../data/data103x1579.txt"); //string trainfile(argv[1]);

	string trainfile("../data/idash18x1579.txt");
	//string trainfile("../data/edin.txt");
	//string trainfile("../data/lbw.txt");
	//string trainfile("../data/nhanes3.txt");
	//string trainfile("../data/pcs.txt"); 
	//string trainfile("../data/uis.txt"); 
	bool isYfirst = 1;                            //atoi(argv[2]);
	long numIter = 1;                             //atol(argv[3]);
	long kdeg = 5;                                //atol(argv[4]);
	double gammaUp = 1;                           //atof(argv[5]);
	double gammaDown = -1;                        //atof(argv[6]);
	bool isInitZero = 1;                          //atoi(argv[7]);

	long fold = 10;                                //atol(argv[8]);
	
	bool isEncrypted = 1;                         //atoi(argv[9]);
	//string testfile = argc > 10? string(argv[10]) : trainfile;
	//-----------------------------------------

	if(argc > 10) {
			cout << "THIS SHOULD NOT HAPPEN!" << endl;
			exit(0);
	} else {
		long sampleDim = 0, factorDim = 0;
		cout << "double** zData = GD::zDataFromFile(trainfile, factorDim, sampleDim, 1);" << endl;
		double** zData = GD::zDataFromFile(trainfile, factorDim, sampleDim, 1); //isYfirst);
		cout << "printData(zData, , ,) " << endl;
		cout << factorDim << " " << sampleDim << endl;
		GD::printData(zData,factorDim, sampleDim);

		//double** zInvB = GD::zInvBFromFile(zData, factorDim, sampleDim, 1); //isYfirst);
		//cout << "printData(zInvB, , ,) " << endl;
		//GD::printData(zInvB,factorDim, sampleDim);

		//GD::shuffleZData(zData, factorDim, sampleDim);
		if(isEncrypted) {
			// # ./MyIDASH2017 "../data/data103x1579.txt" 1 5 5 1 -1 1 5 1
			cout << endl << "# ./MyIDASH2017 \"../data/data103x1579.txt\" 1 5 5 1 -1 1 5 1" << endl << endl << endl;
			TestGD::testEncNLGDFOLD(fold, zData, factorDim, sampleDim, isYfirst, numIter, kdeg, gammaUp, gammaDown, isInitZero);
		} else {
			cout << "THIS SHOULD NOT HAPPEN!" << endl;
			exit(0);
		}
	}


	return 0;
}
