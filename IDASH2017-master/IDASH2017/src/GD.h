
////////////////////////////////////////////////////////////////////////////////////
#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif
////////////////////////////////////////////////////////////////////////////////////

#ifndef IDASH2017_GD_H_
#define IDASH2017_GD_H_

#include <iostream>
#include <assert.h>
#include <cmath>
#include <set>

#include "Ciphertext.h"
#include "NTL/ZZX.h"
#include "Scheme.h"
#include "TestScheme.h"
#include "SecretKey.h"
#include "TimeUtils.h"


static double degree3[3] = {-0.5,0.15012,-0.001593};
//static double degree5[4] = {-0.5,0.19131,-0.0045963, 0.0000412332};
static double degree5[4] = {+0.5, -0.19131, +0.0045963,   -0.0000412332};
static double degree7[5] = {-0.5,0.216884,-0.00819276,0.000165861,-0.00000119581};

//static double degree3[3] = {-0.5,0.19,-0.0035};
//static double degree5[4] = {-0.5,0.2166,-0.0077,0.00011};
//static double degree7[5] = {-0.5,0.216884,-0.00819276,0.000165861,-0.00000119581};

using namespace std;

class GD {

public:

	static double** zDataFromFile(string& path, long& factorDim, long& sampleDim, bool isfirst = true);
	static double** zInvBFromFile(double **zDataTrain, long& factorDim, long& sampleDim, bool isfirst = true, double epsilon = 1e-8);
	static void printData(double** zData, long factorDim, long sampleDim);

	static void shuffleZData(double** zData, long factorDim, long sampleDim);

	static void normalizeZData(double** zData, long factorDim, long sampleDim);
	static void normalizezData2(double** zDataLearn, double** zDataTest, long factorDim, long sampleDimLearn, long sampleDimTest);

	static void initialWDataVDataAverage(double* wData, double* vData, double** zData, long factorDim, long sampleDim);
	static void initialWDataVDataZero(double* wData, double* vData, long factorDim);

	static double* plainIP(double** a, double* b, long factorDim, long sampleDim);
	static double* plainSigmoid(long approxDeg, double** zData, double* ip, long factorDim, long sampleDim, double gamma);

	static void plainLGDstep(double* wData, double* grad, long factorDim);
	static void plainMLGDstep(double* wData, double* vData, double* grad, long factorDim, double eta);
	static void plainNLGDstep(double* wData, double* vData, double* grad, long factorDim, double eta);

	static void plainLGDL2step(double* wData, double* grad, long factorDim, double lambda);
	static void plainMLGDL2step(double* wData, double* vData, double* grad, long factorDim, double eta, double lambda);
	static void plainNLGDL2step(double* wData, double* vData, double* grad, long factorDim, double eta, double lambda);

	static void plainLGDiteration(long approxDeg, double** zData, double* wData, long factorDim, long sampleDim, double gamma);
	static void plainMLGDiteration(long approxDeg, double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta);
	static void plainNLGDiteration(long approxDeg, double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta);

	static void plainLGDL2iteration(long approxDeg, double** zData, double* wData, long factorDim, long sampleDim, double gamma, double lambda);
	static void plainMLGDL2iteration(long approxDeg, double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta, double lambda);
	static void plainNLGDL2iteration(long approxDeg, double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta, double lambda);

	//-----------------------------------------

	static double trueIP(double* a, double* b, long size);

	static void trueLGDiteration(double** zData, double* wData, long factorDim, long sampleDim, double gamma);
	static void trueMLGDiteration(double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta);
	static void trueNLGDiteration(double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta);

	static void trueLGDL2iteration(double** zData, double* wData, long factorDim, long sampleDim, double gamma, double lambda);
	static void trueMLGDL2iteration(double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta, double lambda);
	static void trueNLGDL2iteration(double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta, double lambda);

	static void calculateAUC(double** zData, double* wData, long factorDim, long sampleDim, double& correctness, double& AUC);
	static double calculateMSE(double* wData1, double* wData2, long factorDim);
	static double calculateNMSE(double* wData1, double* wData2, long factorDim);
        static double calculateACC(double** zData, double* wData, long factorDim, long sampleDim, double& correctness, double& auc);
	
	static size_t getPeakRSS();
	static size_t getCurrentRSS();

};

#endif /* SGD_SGD_H_ */
