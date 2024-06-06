#include "CipherGD.h"
#include "GD.h"

#include <Ciphertext.h>
#include <EvaluatorUtils.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>

#include <assert.h>

void CipherGD::encZData(Ciphertext* encZData, double** zData, long slots, long factorDim, long sampleDim, long batch, long cnum, long wBits, long logQ) {
	complex<double>* pzData = new complex<double>[slots];
	for (long i = 0; i < cnum - 1; ++i) {
		for (long j = 0; j < sampleDim; ++j) {
			for (long l = 0; l < batch; ++l) {
				pzData[batch * j + l].real(zData[j][batch * i + l]);
			}
		}
		scheme.encrypt(encZData[i], pzData, slots, wBits, logQ);
	}

	long rest = factorDim - batch * (cnum - 1);
	for (long j = 0; j < sampleDim; ++j) {
		for (long l = 0; l < rest; ++l) {
			pzData[batch * j + l].real(zData[j][batch * (cnum - 1) + l]);
		}
		for (long l = rest; l < batch; ++l) {
			pzData[batch * j + l] = 0;
		}
	}
	scheme.encrypt(encZData[cnum - 1], pzData, slots, wBits, logQ);

	delete[] pzData;
}

void CipherGD::encWDataAverage(Ciphertext* encWData, Ciphertext* encZData, long cnum, long sBits, long bBits) {
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		encWData[i] = encZData[i];
		Ciphertext rot;
		for (long l = bBits; l < sBits; ++l) {
			scheme.leftRotateFast(rot, encWData[i], (1 << l));
			scheme.addAndEqual(encWData[i], rot);
		}
		scheme.divByPo2AndEqual(encWData[i], sBits - bBits);
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encWVDataAverage(Ciphertext* encWData, Ciphertext* encVData, Ciphertext* encZData, long cnum, long sBits, long bBits) {
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		encWData[i] = encZData[i];
		Ciphertext rot;
		for (long l = bBits; l < sBits; ++l) {
			scheme.leftRotateFast(rot, encWData[i], (1 << l));
			scheme.addAndEqual(encWData[i], rot);
		}
		scheme.divByPo2AndEqual(encWData[i], sBits - bBits);
		encVData[i] = encWData[i];
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encWDataZero(Ciphertext* encWData, long cnum, long slots, long wBits, long logQ) {
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		 scheme.encryptZeros(encWData[i], slots, wBits, logQ);
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encWVDataZero(Ciphertext* encWData, Ciphertext* encVData, long cnum, long slots, long wBits, long logQ) {
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.encryptSingle(encWData[i], 0.0, wBits, logQ);
		encWData[i].n = slots;
		encVData[i].copy(encWData[i]);
		encVData[i].n = slots;
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::generateAuxPoly(uint64_t* poly, long slots, long batch, long pBits) {
	complex<double>* pvals = new complex<double>[slots];
	for (long j = 0; j < slots; j += batch) {
		pvals[j].real(1.0);
	}
	ZZ* msg = new ZZ[N];
	scheme.ring.encode(msg, pvals, slots, pBits);
	long np = ceil((pBits + logQ + logN + 2)/59.);
	scheme.ring.CRT(poly, msg, np);
	delete[] pvals;
	delete[] msg;
}

void CipherGD::encInnerProduct(Ciphertext& encIP, Ciphertext* encZData, Ciphertext* encVData, uint64_t* rpoly, long cnum, long bBits, long wBits, long pBits) {
	Ciphertext* encIPvec = new Ciphertext[cnum];

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		assert(encZData[i].logq >=  encVData[i].logq);
        scheme.modDownTo(encIPvec[i], encZData[i], encVData[i].logq);
		scheme.multAndEqual(encIPvec[i], encVData[i]);

		Ciphertext rot;
		for (long l = 0; l < bBits; ++l) {
			scheme.leftRotateFast(rot, encIPvec[i], (1 << l));
			scheme.addAndEqual(encIPvec[i], rot);
		}
		rot.kill();
	}
	NTL_EXEC_RANGE_END


	encIP.copy(encIPvec[0]);
	for (long i = 1; i < cnum; ++i) {
		scheme.addAndEqual(encIP, encIPvec[i]);
	}



	scheme.multByPolyNTTAndEqual(encIP, rpoly, pBits, pBits); 
	Ciphertext tmp;
	for (long l = 0; l < bBits; ++l) {
		scheme.rightRotateFast(tmp, encIP, (1 << l));
		scheme.addAndEqual(encIP, tmp);
	}
	tmp.kill();
	scheme.reScaleByAndEqual(encIP, pBits);

	for (long i=0; i < cnum; ++i) encIPvec[i].kill();
	delete[] encIPvec;

}

void CipherGD::encSigmoid(long kdeg, Ciphertext* encZData, Ciphertext* encGrad, Ciphertext* encBinv, Ciphertext& encIP, long cnum, double gamma, long sBits, long bBits, long pBits, long wBits, long aBits) {

	Ciphertext encIP2(encIP);

	scheme.multAndEqual(encIP2, encIP);

	scheme.reScaleByAndEqual(encIP2, encIP.logp);          

	if(kdeg == 3) {
		cout << "SFESFFE" << endl;
		exit(0);
	} else if (kdeg == 5) {
		Ciphertext encIP4;
		scheme.square(encIP4, encIP2);
		scheme.reScaleByAndEqual(encIP4, encIP2.logp);

		scheme.multByConstAndEqual(encIP2, degree5[2] / degree5[3], wBits);
		scheme.reScaleByAndEqual(encIP2, wBits);


        if(encIP4.logq > encIP2.logq) scheme.modDownToAndEqual(encIP4, encIP2.logq);
		if(encIP4.logq < encIP2.logq) scheme.modDownToAndEqual(encIP2, encIP4.logq);		
		if(encIP4.logp != encIP2.logp) {cout << "ESEF"; exit(0);}
		scheme.addAndEqual(encIP4, encIP2);

		scheme.addConstAndEqual(encIP4, degree5[1] / degree5[3], encIP4.logp);


		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
		//for (long i = 0; i < cnum; ++i) {	
			scheme.multByConst(encGrad[i], encZData[i], gamma * degree5[3], wBits + pBits);
		    scheme.reScaleByAndEqual(encGrad[i], pBits);

		    Ciphertext ctIP(encIP);	
			if(encGrad[i].logq > ctIP.logq)
				scheme.modDownToAndEqual(encGrad[i], ctIP.logq);
			if(encGrad[i].logq < ctIP.logq)
				scheme.modDownToAndEqual(ctIP, encGrad[i].logq);		
			scheme.multAndEqual(encGrad[i], ctIP);
			scheme.reScaleByAndEqual(encGrad[i], ctIP.logp);

			Ciphertext ctIP4(encIP4);
			if(encGrad[i].logq > ctIP4.logq)
				scheme.modDownToAndEqual(encGrad[i], ctIP4.logq);
			if(encGrad[i].logq < ctIP4.logq)
				scheme.modDownToAndEqual(ctIP4, encGrad[i].logq);
			scheme.multAndEqual(encGrad[i], ctIP4);
			scheme.reScaleByAndEqual(encGrad[i], ctIP4.logp);
		

			Ciphertext tmp;
			scheme.multByConst(tmp, encZData[i], gamma * degree5[0], wBits);
			scheme.modDownToAndEqual(tmp, encGrad[i].logq);

			scheme.addAndEqual(encGrad[i], tmp);

			tmp.kill();
			ctIP4.kill();
			ctIP.kill();

		}
		NTL_EXEC_RANGE_END;

		encIP4.kill();

	} else {
		cout << "EFDSFEE" << endl;
		exit(0);

	}


	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		Ciphertext tmp;
		for (long l = bBits; l < sBits; ++l) {
			scheme.leftRotateFast(tmp, encGrad[i], (1 << l));
			scheme.addAndEqual(encGrad[i], tmp);
		}
		tmp.kill();

				Ciphertext ctBinv(encBinv[i]);
					if (encGrad[i].logq > ctBinv.logq)
						scheme.modDownToAndEqual(encGrad[i], ctBinv.logq);
					if (encGrad[i].logq < ctBinv.logq)
						scheme.modDownToAndEqual(ctBinv, encGrad[i].logq);

					scheme.multAndEqual(encGrad[i], encBinv[i]);
					scheme.reScaleByAndEqual(encGrad[i], encBinv[i].logp);
					ctBinv.kill();

	}
	NTL_EXEC_RANGE_END;

}


void CipherGD::encNLGDstep(Ciphertext* encWData, Ciphertext* encVData, Ciphertext* encGrad, double eta, long cnum, long wBits, long pBits) {

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
	//for (long i = 0; i < cnum; ++i) {
		if(encGrad[i].logp > encVData[i].logp) 
			scheme.reScaleByAndEqual(encGrad[i], encGrad[i].logp-encVData[i].logp);
		if(encGrad[i].logp < encVData[i].logp) 
			scheme.reScaleByAndEqual(encVData[i], encVData[i].logp-encGrad[i].logp);
		scheme.modDownToAndEqual(encVData[i], encGrad[i].logq);

		Ciphertext ctmpw;
		scheme.add(ctmpw, encVData[i], encGrad[i]);

		scheme.multByConst(encVData[i], ctmpw, 1. - eta, pBits);

		scheme.multByConstAndEqual(encWData[i], eta, pBits);

		if (encWData[i].logq > encVData[i].logq) 
			scheme.modDownToAndEqual(encWData[i], encVData[i].logq);
		if (encWData[i].logq < encVData[i].logq) 
			scheme.modDownToAndEqual(encVData[i], encWData[i].logq);
		if (encWData[i].logp != encVData[i].logp) 
			{ cout << "logp != logp" ;exit(0); }


		scheme.addAndEqual(encVData[i], encWData[i]);

		scheme.reScaleByAndEqual(encVData[i], pBits);
		encWData[i].copy(ctmpw);

		ctmpw.kill();
	}
	NTL_EXEC_RANGE_END;
}





void CipherGD::encNLGDiteration(long kdeg, Ciphertext* encZData, Ciphertext* encZInvB, Ciphertext* encWData, Ciphertext* encVData, uint64_t* rpoly, long cnum, double gamma, double eta, long sBits, long bBits, long wBits, long pBits, long aBits) {
 	Ciphertext* encGrad = new Ciphertext[cnum];
	Ciphertext encIP;

	encInnerProduct(encIP, encZData, encVData, rpoly, cnum, bBits, wBits, pBits);


	encSigmoid(kdeg, encZData, encGrad, encZInvB, encIP, cnum, gamma, sBits, bBits, pBits, wBits, aBits);

	cout << endl << "encGrad[i].logp = " << encGrad[0].logp << endl;

         
	encNLGDstep(encWData, encVData, encGrad, eta, cnum, wBits, pBits);

	for( long i=0; i < cnum; ++i) encGrad[i].kill();
	delete[] encGrad;
}

void CipherGD::decWData(double* wData, Ciphertext* encWData, long factorDim, long batch, long cnum, long wBits) {
	for (long i = 0; i < (cnum - 1); ++i) {
		complex<double>* dcw = scheme.decrypt(secretKey, encWData[i]);
		for (long j = 0; j < batch; ++j) {
			wData[batch * i + j] = dcw[j].real();
		}
		delete[] dcw;
	}
	complex<double>* dcw = scheme.decrypt(secretKey, encWData[cnum-1]);
	long rest = factorDim - batch * (cnum - 1);
	for (long j = 0; j < rest; ++j) {
		wData[batch * (cnum - 1) + j] = dcw[j].real();
	}
	delete[] dcw;
}
