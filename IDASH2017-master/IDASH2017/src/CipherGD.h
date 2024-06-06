#ifndef IDASH2017_CIPHERGD_H_
#define IDASH2017_CIPHERGD_H_

#include <Scheme.h>
#include <SecretKey.h>

#include <complex>

#define variablename( val )    # val 

using namespace std;
using namespace NTL;

class CipherGD {
public:
	Scheme& scheme;
	SecretKey& secretKey;
	long& factorDim;
	long& batch;
	long& cnum;

	CipherGD(long& cnum, long& batch, long& factorDim, Scheme& scheme, SecretKey& secretKey) : cnum(cnum), batch(batch), factorDim(factorDim), scheme(scheme), secretKey(secretKey) {}

	void encZData(Ciphertext* encZData, double** zData, long slots, long factorDim, long sampleDim, long batch, long cnum, long wBits, long logQ);

	void encWDataAverage(Ciphertext* encWData, Ciphertext* encZData, long cnum, long sBits, long bBits);

	void encWVDataAverage(Ciphertext* encWData, Ciphertext* encVData, Ciphertext* encZData, long cnum, long sBits, long bBits);

	void encWDataZero(Ciphertext* encWData, long cnum, long slots, long wBits, long logQ);

	void encWVDataZero(Ciphertext* encWData, Ciphertext* encVData, long cnum, long slots, long wBits, long logQ);

	void generateAuxPoly(uint64_t* poly, long slots, long batch, long pBits);

	void encInnerProduct(Ciphertext& encIP, Ciphertext* encZData, Ciphertext* encWData, uint64_t* rpoly, long cnum, long bBits, long wBits, long pBits);

	void encSigmoid(long kdeg, Ciphertext* encZData, Ciphertext* encGrad, Ciphertext* encBinv, Ciphertext& encIP, long cnum, double gamma, long sBits, long bBits, long pBits, long wBits, long aBits);

	void encLGDstep(Ciphertext* encWData, Ciphertext* encGrad, long cnum);
	void encMLGDstep(Ciphertext* encWData, Ciphertext* encVData, Ciphertext* encGrad, double eta, long cnum, long pBits);
	void encNLGDstep(Ciphertext* encWData, Ciphertext* encVData, Ciphertext* encGrad, double eta, long cnum, long wBits, long pBits);

	void encLGDiteration(long kdeg, Ciphertext* encZData, Ciphertext* encWData, uint64_t* rpoly, long cnum, double gamma, long sBits, long bBits, long wBits, long pBits, long aBits);
	void encMLGDiteration(long kdeg, Ciphertext* encZData, Ciphertext* encWData, Ciphertext* encVData, uint64_t* rpoly, long cnum, double gamma, double eta, long sBits, long bBits, long wBits, long pBits, long aBits);



	void encNLGDiteration(long kdeg, Ciphertext* encZData, Ciphertext* encZInvB, Ciphertext* encWData, Ciphertext* encVData, uint64_t* rpoly, long cnum, double gamma, double eta, long sBits, long bBits, long wBits, long pBits, long aBits);

	void decWData(double* wData, Ciphertext* encWData, long factorDim, long batch, long cnum, long wBits);
};

#endif /* CIPHERGD_H_ */
