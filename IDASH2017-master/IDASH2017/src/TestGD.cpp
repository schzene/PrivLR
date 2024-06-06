#include "TestGD.h"

#include "Ciphertext.h"
#include "NTL/ZZX.h"
#include "Scheme.h"
#include "SecretKey.h"
#include "TimeUtils.h"
#include <cmath>
#include <NTL/BasicThreadPool.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>

#include "SerializationUtils.h"  // error: ‘SerializationUtils’ has not been declared

#include "CipherGD.h"
#include "GD.h"



long TestGD::suggestLogN(long lambda, long logQ) {
    long NBnd = ceil(logQ * (lambda + 110) / 3.6);
    double logNBnd = log2((double)NBnd);
    return (long)ceil(logNBnd);
}


void TestGD::testPlainNLGD(double** zDataTrain, double** zDataTest, long factorDim, long sampleDimTrain, long sampleDimTest,
        bool isYfirst, long numIter, long kdeg, double gammaUp, double gammaDown, bool isInitZero) {
    cout << "isYfirst = " << isYfirst << ", number of iterations = " << numIter << ", g_k = " << kdeg << endl;
    cout << "factors = " << factorDim << ", train samples = " << sampleDimTrain << ", test samples = " << sampleDimTest << endl;
    cout << "gammaUp = " << gammaUp << ", gammaDown = " << gammaDown << ", isInitZero = " << isInitZero << endl;

    TimeUtils timeutils;

    double* pwData = new double[factorDim];
    double* pvData = new double[factorDim];

    double* twData = new double[factorDim];
    double* tvData = new double[factorDim];

    GD::normalizezData2(zDataTrain, zDataTest, factorDim, sampleDimTrain, sampleDimTest);

    if(isInitZero) {
        GD::initialWDataVDataZero(pwData, pvData, factorDim);
        GD::initialWDataVDataZero(twData, tvData, factorDim);
    } else {
        GD::initialWDataVDataAverage(pwData, pvData, zDataTrain, factorDim, sampleDimTrain);
        GD::initialWDataVDataAverage(twData, tvData, zDataTrain, factorDim, sampleDimTrain);
    }

    double alpha0, alpha1, eta, gamma;
    double plaincor, plainauc, truecor, trueauc;

    alpha0 = 0.01;
    alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

    for (long iter = 0; iter < numIter; ++iter) {
        eta = (1 - alpha0) / alpha1;
        gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDimTrain : gammaUp / (iter - gammaDown) / sampleDimTrain;

        GD::plainNLGDiteration(kdeg, zDataTrain, pwData, pvData, factorDim, sampleDimTrain, gamma, eta);
        GD::trueNLGDiteration(zDataTrain, twData, tvData, factorDim, sampleDimTrain, gamma, eta);

        alpha0 = alpha1;
        alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
    }

    cout << "------PLAIN-------" << endl;
    GD::calculateAUC(zDataTest, pwData, factorDim, sampleDimTest, plaincor, plainauc);
    cout << "------------------" << endl;

    cout << "-------TRUE-------" << endl;
    GD::calculateAUC(zDataTest, twData, factorDim, sampleDimTest, truecor, trueauc);
    cout << "------------------" << endl;

//  GD::calculateMSE(twData, pwData, factorDim);
//  GD::calculateNMSE(twData, pwData, factorDim);

}

void TestGD::testEncNLGDFOLD(long fold, double** zData, long factorDim, long sampleDim,
            bool isYfirst, long numIter, long kdeg, double gammaUp, double gammaDown, bool isInitZero) {
    cout << endl << "testEncNLGDFOLD( ) " << endl << endl;
    cout << "isYfirst = " << isYfirst << ", number of iterations = " << numIter << ", g_k = " << kdeg << endl;
    cout << "factors = " << factorDim << ", samples = " << sampleDim << ", fold = " << fold << endl;
    cout << "gammaUp = " << gammaUp << ", gammaDown = " << gammaDown << ", isInitZero = " << isInitZero << endl;

    long sampleDimTest = sampleDim / fold;
    long sampleDimTrain = sampleDim - sampleDimTest;

    long fdimBits = (long)ceil(log2(factorDim));
    long sdimBits = (long)ceil(log2(sampleDimTrain));

    long wBits = 30;
    long pBits = 20;
    long lBits = 5;
    long aBits = 3;
    long kBits = (long)ceil(log2(kdeg));

//  long logQ = isInitZero ? (wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits) :
//          (sdimBits + wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits);

//  long logN = TestGD::suggestLogN(80, logQ);
    //long logN = 16;
        //long logQ = 990; // 991.300840336 > logQ  to ensure 128-bit security level. Security Parameter λ
    long bBits = min(logN - 1 - sdimBits, fdimBits);
    long batch = 1 << bBits;
    long sBits = sdimBits + bBits;
    long slots =  1 << sBits;
    long cnum = (long)ceil((double)factorDim / batch);

    cout << "batch = " << batch << ", slots = " << slots << ", cnum = " << cnum << endl;

    
    cout << "logQ = " << logQ << ", logN = " << logN << endl;
    
//  cout << "HEAAN PARAMETER logQ: " << logQ << endl;
//  cout << "HEAAN PARAMETER logN: " << logN << endl;

    TimeUtils timeutils;
    timeutils.start("Scheme generating...");
    Ring ring;
    SecretKey secretKey(ring);
    Scheme scheme(secretKey, ring);
    scheme.addLeftRotKeys(secretKey);
    scheme.addRightRotKeys(secretKey);
    timeutils.stop("Scheme generation");
    CipherGD cipherGD(cnum, batch, factorDim, scheme, secretKey);

    long logT=3, logI=4;
    timeutils.start("Bootstrap Key generating");
    long bootlogq = 30  +10;
    long lognslots = (long)ceil(log2(batch));  //batch = factorDim / cnum;
    //scheme.addBootKey(secretKey, logn, logq+logI);
    scheme.addBootKey(secretKey, lognslots, bootlogq +logI);
    timeutils.stop("Bootstrap Key generated");


    timeutils.start("Polynomial generating...");
    long np = ceil((pBits + logQ + logN + 2)/59.);
    uint64_t* rpoly = new uint64_t[np << logN];
    cipherGD.generateAuxPoly(rpoly, slots, batch, pBits);
    timeutils.stop("Polynomial generation");

    double* pwData = new double[factorDim];
    double* pvData = new double[factorDim];

    double* twData = new double[factorDim];
    double* tvData = new double[factorDim];

    double* cwData = new double[factorDim];
    double* cvData = new double[factorDim];

    Ciphertext* encZData = new Ciphertext[cnum];
    Ciphertext* encZInvB = new Ciphertext[cnum];
    Ciphertext* encWData = new Ciphertext[cnum];
    Ciphertext* encVData = new Ciphertext[cnum];

    double **zDataTrain, **zDataTest;

    zDataTrain = new double*[sampleDimTrain];
    zDataTest = new double*[sampleDimTest];

    GD::shuffleZData(zData, factorDim, sampleDim);

    double enccor, encauc, truecor, trueauc;
    double averenccor = 0, averencauc = 0, avertruecor = 0, avertrueauc = 0;
    double averevalutime = 0;

    for (long fnum = 0; fnum < fold; ++fnum) {
        cout << " !!! START " << fnum + 1 << " FOLD !!! " << endl;

        for (long i = 0; i < sampleDimTest; ++i) {
            zDataTest[i] = zData[fnum * sampleDimTest + i];
        }
        for (long j = 0; j < fnum; ++j) {
            for (long i = 0; i < sampleDimTest; ++i) {
                zDataTrain[j * sampleDimTest + i] = zData[j * sampleDimTest + i];
            }
        }
        for (long i = (fnum + 1) * sampleDimTest; i < sampleDim; ++i) {
            zDataTrain[i - sampleDimTest] = zData[i];
        }

        double** zInvB = GD::zInvBFromFile(zDataTrain, factorDim, sampleDimTrain, 1); //isYfirst);
        cout << endl << "zInvB : " << endl;
        GD::printData(zInvB, factorDim, sampleDimTrain);

        
        timeutils.start("Encrypting zInvB...");
        cipherGD.encZData(encZInvB, zInvB, slots, factorDim, sampleDimTrain, batch, cnum, wBits+10, logQ);
        timeutils.stop("zInvB encryption");
        for (long i = 0; i < cnum; ++i) {
            SerializationUtils::writeCiphertext(encZInvB[i], "encZInvB["+ std::to_string(i) +"].txt");
        }

        timeutils.start("Encrypting zData...");
        cipherGD.encZData(encZData, zDataTrain, slots, factorDim, sampleDimTrain, batch, cnum, wBits, logQ);
        timeutils.stop("zData encryption");
        for (long i = 0; i < cnum; ++i) {
            SerializationUtils::writeCiphertext(encZData[i], "encZData["+ std::to_string(i) +"].txt");
        }

        timeutils.start("Encrypting wData and vData...");
        if(isInitZero) {
            cipherGD.encWVDataZero(encWData, encVData, cnum, slots, wBits, logQ);
            for (long i = 0; i < cnum; ++i) {
                SerializationUtils::writeCiphertext(encWData[i], "encWData["+ std::to_string(i) +"].txt");
            }
        } else {
            cout << "THIS SHOULD NOT HAPPEN!" << endl;
            exit(0);
        }
        timeutils.stop("wData and vData encryption");

        if(isInitZero) {
            GD::initialWDataVDataZero(pwData, pvData, factorDim);
            GD::initialWDataVDataZero(twData, tvData, factorDim);
        } else {
            cout << "THIS SHOULD NOT HAPPEN!" << endl;
            exit(0);
        }

        //-----------------------------------------

        double alpha0, alpha1, eta, gamma;

        alpha0 = 0.01;
        alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
        double alliterationtime = 0;
        for (long iter = 0; iter < numIter; ++iter) {
            cout << " !!! START " << iter + 1 << " ITERATION !!! " << endl;
            eta = (1 - alpha0) / alpha1;
            //gamma = 10.0 / (iter + 1) / sampleDimTrain;
	    gamma = pow(0.9, iter);

            cout << "encWData.logq before: " << encWData[0].logq << endl;
            timeutils.start("Enc NLGD");
            cipherGD.encNLGDiteration(kdeg, encZData, encZInvB, encWData, encVData, rpoly, cnum, 1.0 + gamma, eta, sBits, bBits, wBits, pBits, aBits);
            timeutils.stop("Enc NLGD");
            cout << "encWData.logq after: " << encWData[0].logq << endl;
            alliterationtime += timeutils.timeElapsed;

            cout << "----ENCRYPTED--encVData---" << endl;
            cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits);
            GD::calculateAUC(zDataTest, cwData, factorDim, sampleDimTest, enccor, encauc);
            cout << "------------------" << endl;

            GD::plainNLGDiteration(kdeg, zDataTrain, pwData, pvData, factorDim, sampleDimTrain, gamma, eta);
            GD::trueNLGDiteration(zDataTrain, twData, tvData, factorDim, sampleDimTrain, gamma, eta);

        //if ( encVData[0].logq <= 310 + wBits && iter < numIter-1 || encVData[0].logq < wBits && iter == numIter-1 ) {
        if ( 0 ) {
            
            timeutils.start("Use Bootstrap To Recrypt Ciphertext");         
                NTL_EXEC_RANGE(cnum, first, last);
                for(long i = first; i < last; ++i){
                    scheme.modDownToAndEqual(encWData[i], bootlogq);
                    encWData[i].n = batch;
                    scheme.bootstrapAndEqual(encWData[i], bootlogq, logQ, logT, logI);
                    encWData[i].n = slots;
                 }
                NTL_EXEC_RANGE_END
                NTL_EXEC_RANGE(cnum, first, last);
                for(long i = first; i < last; ++i){
                    scheme.modDownToAndEqual(encVData[i], bootlogq);
                    encVData[i].n = batch;
                    scheme.bootstrapAndEqual(encVData[i], bootlogq, logQ, logT, logI);
                    encVData[i].n = slots;
                 }
                NTL_EXEC_RANGE_END              

            timeutils.stop("Use Bootstrap To Recrypt Ciphertext");

        }


            alpha0 = alpha1;
            alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
            cout << " !!! STOP " << iter + 1 << " ITERATION !!! " << endl;
            cout << "CurrentRSS (MB): " << ( GD::getCurrentRSS() /1024.0/1024.0 ) << endl;
                cout << "PeakRSS    (MB): " << ( GD::getPeakRSS() /1024.0/1024.0 )    << endl;
        }
        averevalutime += alliterationtime;

        averenccor += enccor;
        averencauc += encauc;
        avertruecor += truecor;
        avertrueauc += trueauc;

        cout << endl << endl << endl << endl;
        cout << "----ENCRYPTED-----" << endl;
        cipherGD.decWData(cvData, encVData, factorDim, batch, cnum, wBits);
        GD::calculateAUC(zDataTest, cvData, factorDim, sampleDimTest, enccor, encauc);
        cout << "------------------" << endl;

        cout << "-------TRUE-------" << endl;
        GD::calculateAUC(zDataTest, tvData, factorDim, sampleDimTest, truecor, trueauc);
        cout << "------------------" << endl;


        cout << " !!! STOP " << fnum + 1 << " FOLD !!! " << endl;
        cout << "------------------" << endl;
        cout << "CurrentRSS (MB): " << ( GD::getCurrentRSS() /1024.0/1024.0 ) << endl;
            cout << "PeakRSS    (MB): " << ( GD::getPeakRSS() /1024.0/1024.0 )    << endl;
    }

    cout << "Average Encrypted correctness: " << averenccor/fold  << "%" << endl;
    cout << "Average Encrypted AUC: " << averencauc/fold << endl;
    cout << "Average True correctness: " << avertruecor/fold << "%" << endl;
    cout << "Average True AUC: " << avertrueauc/fold << endl;
    cout << "Average Evaluation Time: " << averevalutime/fold << endl;
    
    cout << "CurrentRSS (MB): " << ( GD::getCurrentRSS() /1024.0/1024.0 ) << endl;
        cout << "PeakRSS    (MB): " << ( GD::getPeakRSS() /1024.0/1024.0 )    << endl;
}

void TestGD::testPlainNLGDFOLD(long fold, double** zData, long factorDim, long sampleDim,
            bool isYfirst, long numIter, long kdeg, double gammaUp, double gammaDown, bool isInitZero) {
    cout << "isYfirst = " << isYfirst << ", number of iterations = " << numIter << ", g_k = " << kdeg << endl;
    cout << "factors = " << factorDim << ", samples = " << sampleDim << ", fold = " << fold << endl;
    cout << "gammaUp = " << gammaUp << ", gammaDown = " << gammaDown << ", isInitZero = " << isInitZero << endl;

    long sampleDimTest = sampleDim / fold;
    long sampleDimTrain = sampleDim - sampleDimTest;

    TimeUtils timeutils;

    double* pwData = new double[factorDim];
    double* pvData = new double[factorDim];

    double* twData = new double[factorDim];
    double* tvData = new double[factorDim];

    double **zDataTrain, **zDataTest;

    zDataTrain = new double*[sampleDimTrain];
    zDataTest = new double*[sampleDimTest];

    GD::normalizeZData(zData, factorDim, sampleDim);
    GD::shuffleZData(zData, factorDim, sampleDim);

    double plaincor, plainauc, truecor, trueauc;
    double averplaincor = 0, averplainauc = 0, avertruecor = 0, avertrueauc = 0;

    for (long fnum = 0; fnum < fold; ++fnum) {
        cout << " !!! START " << fnum + 1 << " FOLD !!! " << endl;

        for (long i = 0; i < sampleDimTest; ++i) {
            zDataTest[i] = zData[fnum * sampleDimTest + i];
        }
        for (long j = 0; j < fnum; ++j) {
            for (long i = 0; i < sampleDimTest; ++i) {
                zDataTrain[j * sampleDimTest + i] = zData[j * sampleDimTest + i];
            }
        }
        for (long i = (fnum + 1) * sampleDimTest; i < sampleDim; ++i) {
            zDataTrain[i - sampleDimTest] = zData[i];
        }

        if(isInitZero) {
            GD::initialWDataVDataZero(pwData, pvData, factorDim);
            GD::initialWDataVDataZero(twData, tvData, factorDim);
        } else {
            GD::initialWDataVDataAverage(pwData, pvData, zDataTrain, factorDim, sampleDimTrain);
            GD::initialWDataVDataAverage(twData, tvData, zDataTrain, factorDim, sampleDimTrain);
        }

        //-----------------------------------------

        double alpha0, alpha1, eta, gamma;

        alpha0 = 0.01;
        alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

        for (long iter = 0; iter < numIter; ++iter) {
            eta = (1 - alpha0) / alpha1;
            cout << eta << endl;
            gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDimTrain : gammaUp / (iter - gammaDown) / sampleDimTrain;
            GD::plainNLGDiteration(kdeg, zDataTrain, pwData, pvData, factorDim, sampleDimTrain, gamma, eta);
            GD::trueNLGDiteration(zDataTrain, twData, tvData, factorDim, sampleDimTrain, gamma, eta);

            alpha0 = alpha1;
            alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
        }
        cout << "------PLAIN-------" << endl;
        GD::calculateAUC(zDataTest, pwData, factorDim, sampleDimTest, plaincor, plainauc);
        cout << "------------------" << endl;

        cout << "-------TRUE-------" << endl;
        GD::calculateAUC(zDataTest, twData, factorDim, sampleDimTest, truecor, trueauc);
        cout << "------------------" << endl;

        averplaincor += plaincor;
        averplainauc += plainauc;
        avertruecor += truecor;
        avertrueauc += trueauc;

        cout << " !!! STOP " << fnum + 1 << " FOLD !!! " << endl;
        cout << "------------------" << endl;
    }

    cout << "Average Plain correctness: " << averplaincor << "%" << endl;
    cout << "Average Plain AUC: " << averplainauc << endl;
    cout << "Average True correctness: " << avertruecor << "%" << endl;
    cout << "Average True AUC: " << avertrueauc << endl;
}

