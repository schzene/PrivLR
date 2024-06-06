#include "GD.h"

#include <NTL/ZZ.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>     //setw
#include <string.h>
#include <string>
#include <vector>
#include <algorithm>    // std::shuffle + std::copy
#include <array>        // std::array
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock


using namespace NTL;

double** GD::zDataFromFile(string& path, long& factorDim, long& sampleDim, bool isfirst) {
	vector<vector<double>> zline;
	factorDim = 1; 	// dimension of x
	sampleDim = 0;	// number of samples
	ifstream openFile(path.data());
	if(openFile.is_open()) {
		string line, temp;
		getline(openFile, line);
		long i;
		for(i = 0; i < line.length(); ++i) if(line[i] == ',' ) factorDim++;

		while(getline(openFile, line)){
			size_t start=0, end=0;
			vector<double> vecline;
			do {
				end = line.find_first_of (',', start);
				temp = line.substr(start,end);
				vecline.push_back(atof(temp.c_str()));
				start = end + 1;
			} while(start);
			zline.push_back(vecline);
			sampleDim++;
		}
	} else {
		cout << "Error: cannot read file" << endl;
	}


	for (long i = 1; i < factorDim; ++i) {
		double rowmin = zline[0][i], rowmax = zline[0][i];
		for (long j = 0; j < sampleDim; ++j) {
			rowmin = std::min(rowmin, zline[j][i]);  //std::abs(
			rowmax = std::max(rowmax, zline[j][i]);
		}

		if(rowmax - rowmin < 1e-10) {
			for (long j = 0; j < sampleDim; ++j) 
				zline[j][i] = 1.0;
		}else{
			for (long j = 0; j < sampleDim; ++j) 
				zline[j][i] = (zline[j][i]-rowmin)/(rowmax - rowmin);
		}

	}

	double** zData = new double*[sampleDim];
	if(isfirst) {
		for(long j = 0; j < sampleDim; ++j){
			double* zj = new double[factorDim];
			zj[0] = 2 * zline[j][0] - 1;
			for(long i = 1; i < factorDim; ++i){
				zj[i] = zj[0] * zline[j][i];
			}
			zData[j] = zj;
		}
	} else {
		cout << "THIS SHOULD NOT HAPPEN! " << endl;
		exit(0);

	}
	return zData;
}

double** GD::zInvBFromFile(double **zDataTrain, long& factorDim, long& sampleDim, bool isfirst, double epsilon) {
	vector<vector<double>> zline;

	double** zData =  new double*[sampleDim];;
	double** zInvB =  new double*[sampleDim];;
	if(isfirst) {
		for(long j = 0; j < sampleDim; ++j){
			double* zj = new double[factorDim];
			double* Bj = new double[factorDim];
			for(long i = 0; i < factorDim; ++i){
				zj[i] = std::abs(zDataTrain[j][i]);
				Bj[i] = std::abs(zDataTrain[j][i]);
			}
			zData[j] = zj;
			zInvB[j] = Bj;
		}

		for(long j = 0; j < sampleDim; ++j){
			for(long i = 1; i < factorDim; ++i){
				zInvB[j][0] += zInvB[j][i];
			}
		}

		for(long j = 0; j < sampleDim; ++j){
			for(long i = 1; i < factorDim; ++i){
				zInvB[j][i] = zInvB[j][0];
			}
		}

		for(long j = 0; j < sampleDim; ++j){
			for(long i = 0; i < factorDim; ++i){
				zInvB[j][i] *= zData[j][i];
			}
		}

		for(long i = 0; i < factorDim; ++i){
			for(long j = 1; j < sampleDim; ++j){
				zInvB[0][i] += zInvB[j][i];
			}
		}

		for(long i = 0; i < factorDim; ++i){
			for(long j = 1; j < sampleDim; ++j){
				zInvB[j][i] = zInvB[0][i];
			}
		}

		for(long i = 0; i < factorDim; ++i){
			for(long j = 0; j < sampleDim; ++j){
				zInvB[j][i] = 1.0 / (epsilon + .25 * zInvB[j][i]);
			}
		}

		return zInvB;

	} else {
		cout << "THIS SHOULD NOT HAPPEN !" << endl;
		for(long j = 0; j < sampleDim; ++j){
			double* zj = new double[factorDim];
			zj[0] = 2 * zline[j][factorDim - 1] - 1;
			for(long i = 1; i < factorDim; ++i){
				zj[i] = zj[0] * zline[j][i-1];
			}
			zData[j] = zj;
		}
		return zData;
	}
	
}

void GD::printData(double** zData, long factorDim, long sampleDim) {
	for (long i = 0; i < (sampleDim>2?2:sampleDim); ++i) {
		cout << "The " << std::setw(3) << i << "-th Row:";
		for (long j = 0; j < factorDim; ++j) {
			cout << std::showpos << std::fixed << std::setw(12)  << zData[i][j] ;
		}
		cout << endl << endl;
	}
}

void GD::shuffleZData(double** zData, long factorDim, long sampleDim) {
	srand(time(NULL));
	double* tmp = new double[factorDim];
	for (long i = 0; i < sampleDim; ++i) {
		long idx = i + rand() / (RAND_MAX / (sampleDim - i) + 1);
		copy(zData[i], zData[i] + factorDim, tmp);
		copy(zData[idx], zData[idx] + factorDim, zData[i]);
		copy(tmp, tmp + factorDim, zData[idx]);
	}
}

void GD::normalizeZData(double** zData, long factorDim, long sampleDim) {
	long i, j;
	double m;
	for (i = 0; i < factorDim; ++i) {
		m = 0.0;
		for (j = 0; j < sampleDim; ++j) {
			m = max(m, abs(zData[j][i]));
		}

		if(m < 1e-10) continue;

		for (j = 0; j < sampleDim; ++j) {
			zData[j][i] /= m;
		}
	}
}

void GD::normalizezData2(double** zDataLearn, double** zDataTest, long factorDim, long sampleDimLearn, long sampleDimTest) {
	long i, j;
	double m;
	for (i = 0; i < factorDim; ++i) {
		m = 0.0;
		for (j = 0; j < sampleDimLearn; ++j) {
			m = max(m, abs(zDataLearn[j][i]));
		}
		for (j = 0; j < sampleDimTest; ++j) {
			m = max(m, abs(zDataTest[j][i]));
		}
		if(m < 1e-10) continue;
		for (j = 0; j < sampleDimLearn; ++j) {
			zDataLearn[j][i] /= m;
		}
		for (j = 0; j < sampleDimTest; ++j) {
			zDataTest[j][i] /= m;
		}
	}
}

void GD::initialWDataVDataAverage(double* wData, double* vData, double** zData, long factorDim, long sampleDim) {
	long sdimBits = (long)ceil(log2(sampleDim));
	long sdimPow = 1 << sdimBits;
	for (long i = 0; i < factorDim; ++i) {
		double tmp = 0.0;
		for (long j = 0; j < sampleDim; ++j) {
			tmp += zData[j][i];
		}
		tmp /= sdimPow;
		wData[i] = tmp;
		vData[i] = tmp;
	}
}

void GD::initialWDataVDataZero(double* wData, double* vData, long factorDim) {
	for (long i = 0; i < factorDim; ++i) {
		wData[i] = 0.0;
		vData[i] = 0.0;
	}
}

double* GD::plainIP(double** a, double* b, long factorDim, long sampleDim) {
	double* res = new double[sampleDim]();
	for (long j = 0; j < sampleDim; ++j) {
		for(long i = 0; i < factorDim; ++i) {
			res[j] += a[j][i] * b[i];
		}
	}
	return res;
}

double* GD::plainSigmoid(long approxDeg, double** zData, double* ip, long factorDim, long sampleDim, double gamma) {
	double* grad = new double[factorDim]();
	if(approxDeg == 3) {
		for (long i = 0; i < factorDim; ++i) {
			for (long j = 0; j < sampleDim; ++j) {
				grad[i] += (degree3[0] + ip[j] * degree3[1] + pow(ip[j], 3) * degree3[2]) * zData[j][i];
			}
			grad[i] *= gamma;
		}
	} else if(approxDeg == 5) {
		for (long i = 0; i < factorDim; ++i) {
			for (long j = 0; j < sampleDim; ++j) {
				grad[i] += (degree5[0] + ip[j] * degree5[1] + pow(ip[j], 3) * degree5[2] + pow(ip[j], 5) * degree5[3]) * zData[j][i];
			}
			grad[i] *= gamma;
		}
	} else {
		for (long i = 0; i < factorDim; ++i) {
			for (long j = 0; j < sampleDim; ++j) {
				grad[i] += (degree7[0] + ip[j] * degree7[1] + pow(ip[j], 3) * degree7[2] + pow(ip[j], 5) * degree7[3] + pow(ip[j], 7) * degree7[4]) * zData[j][i];
			}
			grad[i] *= gamma;
		}
	}
	return grad;
}

void GD::plainLGDstep(double* wData, double* grad, long factorDim) {
	for (long i = 0; i < factorDim; ++i) {
		wData[i] -= grad[i];
	}
}

void GD::plainMLGDstep(double* wData, double* vData, double* grad, long factorDim, double eta) {
	for (long i = 0; i < factorDim; ++i) {
		vData[i] = eta * vData[i] + grad[i];
		wData[i] -= vData[i];
	}
}

void GD::plainNLGDstep(double* wData, double* vData, double* grad, long factorDim, double eta) {
	for (long i = 0; i < factorDim; ++i) {
		double tmpw = vData[i] - grad[i];
		vData[i] = (1.0 - eta) * tmpw + eta * wData[i];
		wData[i] = tmpw;
	}
}

void GD::plainLGDL2step(double* wData, double* grad, long factorDim, double lambda) {
	//TODO: implement method
}

void GD::plainMLGDL2step(double* wData, double* vData, double* grad, long factorDim, double eta, double lambda) {
	//TODO: implement method
}

void GD::plainNLGDL2step(double* wData, double* vData, double* grad, long factorDim, double eta, double lambda) {
	//TODO: implement method
}

void GD::plainLGDiteration(long approxDeg, double** zData, double* wData, long factorDim, long sampleDim, double gamma) {
	double* ip = plainIP(zData, wData, factorDim, sampleDim);
	double* grad = plainSigmoid(approxDeg, zData, ip, factorDim, sampleDim, gamma);
	plainLGDstep(wData, grad, factorDim);
	delete[] ip;
	delete[] grad;
}

void GD::plainMLGDiteration(long approxDeg, double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta) {
	double* ip = plainIP(zData, wData, factorDim, sampleDim);
	double* grad = plainSigmoid(approxDeg, zData, ip, factorDim, sampleDim, gamma);
	plainMLGDstep(wData, vData, grad, factorDim, eta);
	delete[] ip;
	delete[] grad;
}


void GD::plainNLGDiteration(long approxDeg, double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta) {
	double* ip = plainIP(zData, vData, factorDim, sampleDim);
	double* grad = plainSigmoid(approxDeg, zData, ip, factorDim, sampleDim, gamma);
	plainNLGDstep(wData, vData, grad, factorDim, eta);

	delete[] ip;
	delete[] grad;
}

void GD::plainLGDL2iteration(long approxDeg, double** zData, double* wData, long factorDim, long sampleDim, double gamma, double lambda) {
	//TODO: implement method
}

void GD::plainMLGDL2iteration(long approxDeg, double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta, double lambda) {
	//TODO: implement method
}

void GD::plainNLGDL2iteration(long approxDeg, double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta, double lambda) {
	//TODO: implement method
}

//-----------------------------------------

double GD::trueIP(double* a, double* b, long size) {
	double res = 0.0;
	for(long i = 0; i < size; ++i) {
		res += a[i] * b[i];
	}
	return res;
}

void GD::trueLGDiteration(double** zData, double* wData, long factorDim, long sampleDim, double gamma) {
	double* grad = new double[factorDim]();

	for(long j = 0; j < sampleDim; ++j) {
		double ip = trueIP(wData, zData[j], factorDim);
		double tmp = (ip > 15.0) ? 0 : (ip < -15.0) ? -1.0 : - 1. / (1. + exp(ip));
		for(int i = 0; i < factorDim; ++i) {
			grad[i] += tmp * (double) zData[j][i];
		}
	}

	for (int i = 0; i < factorDim; ++i) {
		wData[i] -= gamma * grad[i];
	}
	delete[] grad;
}

void GD::trueMLGDiteration(double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta) {
	double* grad = new double[factorDim]();

	for(long j = 0; j < sampleDim; ++j) {
		double ip = trueIP(wData, zData[j], factorDim);
		double tmp = (ip > 15.0) ? 0 : (ip < -15.0) ? -1.0 : - 1. / (1. + exp(ip));
		for(long i = 0; i < factorDim; ++i) {
			grad[i] += tmp * (double) zData[j][i];
		}
	}

	for (int i = 0; i < factorDim; ++i) {
		vData[i] = eta * vData[i] + gamma * grad[i];
		wData[i] -= vData[i];
	}
	delete[] grad;
}

void GD::trueNLGDiteration(double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta) {
	double* grad = new double[factorDim]();

	for(long j = 0; j < sampleDim; ++j) {
		double ip = trueIP(vData, zData[j], factorDim);
		double tmp = (ip > 15.0) ? 0 : (ip < -15.0) ? -1.0 : - 1. / (1. + exp(ip));
		for(long i = 0; i < factorDim; ++i) {
			grad[i] += tmp * (double) zData[j][i];
		}
	}
	for (long i = 0; i < factorDim; ++i) {
		double tmpw = vData[i] - gamma * grad[i];
		vData[i] = (1.0 - eta) * tmpw + eta * wData[i];
		wData[i] = tmpw;
	}
	delete[] grad;
}

void GD::trueLGDL2iteration(double** zData, double* wData, long factorDim, long sampleDim, double gamma, double lambda) {
	//TODO: implement method
}

void GD::trueMLGDL2iteration(double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta, double lambda) {
	//TODO: implement method
}

void GD::trueNLGDL2iteration(double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta, double lambda) {
	//TODO: implement method
}

void GD::calculateAUC(double** zData, double* wData, long factorDim, long sampleDim, double& correctness, double& auc) {
	cout << "w:";
	for (long i = 0; i < factorDim; ++i) {
		cout << wData[i] << ",";
	}
	cout << endl;

	long TN = 0, FP = 0;

    vector<double> thetaTN;
    vector<double> thetaFP;

    for(int i = 0; i < sampleDim; ++i){
        if(zData[i][0] > 0){
            if(GD::trueIP(zData[i], wData, factorDim) < 0) TN++;
            thetaTN.push_back(zData[i][0] * GD::trueIP(zData[i] + 1, wData + 1, factorDim - 1));
        } else{
            if(GD::trueIP(zData[i], wData, factorDim) < 0) FP++;
            thetaFP.push_back(zData[i][0] * GD::trueIP(zData[i] + 1, wData + 1, factorDim - 1));
        }
    }

    correctness = 100.0 - (100.0 * (FP + TN) / sampleDim);
//    cout << "Failure rate: (y = 1) " << TN << "/" << thetaTN.size() << " + (y = 0) " << FP << "/" ;
//    cout << thetaFP.size() << " = " <<  (100.0 * (FP + TN) / sampleDim) << " %." << endl;
    cout << "ACC: " << correctness  << " %." << endl;

    if(thetaFP.size() == 0 || thetaTN.size() == 0) {
        cout << "n_test_yi = 0 : cannot compute AUC" << endl;
        auc = 0.0;
    } else{
        auc = 0.0;
        for(long i = 0; i < thetaTN.size(); ++i){
            for(long j = 0; j < thetaFP.size(); ++j){
                if(thetaFP[j] <= thetaTN[i]) auc++;
            }
        }
        auc /= thetaTN.size() * thetaFP.size();
        cout << "AUC: " << auc << endl << endl;
    }
}


double GD::calculateMSE(double* wData1, double* wData2, long factorDim) {
    double res= 0.0;

    for(long i = 0; i < factorDim; ++i) {
        res += pow(wData1[i] - wData2[i], 2.0);
    }
    res /= factorDim;
    cout << "MSE = " << res << endl;
    return res;
}


double GD::calculateNMSE(double* wData1, double* wData2, long factorDim) {
    double res= 0.0;

    for(long i = 0; i < factorDim; ++i) {
        res += pow(wData1[i], 2.0);
    }
    res /= factorDim;

    double mse = GD::calculateMSE(wData1, wData2, factorDim);
    res = mse / res;

    cout << "NMSE = " << res << endl;
    return res;

}


double GD::calculateACC(double** zData, double* wData, long factorDim, long sampleDim, double& correctness, double& auc) {
//	cout << "w:";
//	for (long i = 0; i < factorDim; ++i) {
//		cout << wData[i] << ",";
//	}
//	cout << endl;

	long TN = 0, FP = 0;

    vector<double> thetaTN;
    vector<double> thetaFP;

    for(int i = 0; i < sampleDim; ++i){
        if(zData[i][0] > 0){
            if(trueIP(zData[i], wData, factorDim) < 0) TN++;
            thetaTN.push_back(zData[i][0] * trueIP(zData[i] + 1, wData + 1, factorDim - 1));
        } else{
            if(trueIP(zData[i], wData, factorDim) < 0) FP++;
            thetaFP.push_back(zData[i][0] * trueIP(zData[i] + 1, wData + 1, factorDim - 1));
        }
    }

    correctness = 100.0 - (100.0 * (FP + TN) / sampleDim);

    if(thetaFP.size() == 0 || thetaTN.size() == 0) {
        cout << "n_test_yi = 0 : cannot compute ACC" << endl;
        auc = 0.0;
    } else{
        auc = 0.0;
        for(long i = 0; i < thetaTN.size(); ++i){
            for(long j = 0; j < thetaFP.size(); ++j){
                if(thetaFP[j] <= thetaTN[i]) auc++;
            }
        }
        auc /= thetaTN.size() * thetaFP.size();
        //cout << "AUC: " << auc << endl;
    }

    cout << "ACC: " << correctness << endl;
    return correctness;
}



/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t GD::getPeakRSS() {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
        return (size_t)0L;      /* Can't open? */
    if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
    {
        close( fd );
        return (size_t)0L;      /* Can't read? */
    }
    close( fd );
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
	/* BSD, Linux, and OSX -------------------------------------- */
	struct rusage rusage;
	getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
	return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;          /* Unsupported. */
#endif
}

/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t GD::getCurrentRSS() {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
        (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
        return (size_t)0L;      /* Can't access? */
    return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
	/* Linux ---------------------------------------------------- */
	long rss = 0L;
	FILE *fp = NULL;
	if ((fp = fopen("/proc/self/statm", "r")) == NULL)
		return (size_t) 0L; /* Can't open? */
	if (fscanf(fp, "%*s%ld", &rss) != 1) {
		fclose(fp);
		return (size_t) 0L; /* Can't read? */
	}
	fclose(fp);
	return (size_t) rss * (size_t) sysconf(_SC_PAGESIZE);

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t)0L;          /* Unsupported. */
#endif
}
