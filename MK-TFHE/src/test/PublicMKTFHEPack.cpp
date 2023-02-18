#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "tfhe.h"
#include "polynomials.h"
#include "lwesamples.h"
#include "lwekey.h"
#include "lweparams.h"
#include "tlwe.h"
#include "tgsw.h"
#include "secretshare.hpp"


#include "mkTFHEparams.h"
#include "mkTFHEkeys.h"
#include "mkTFHEkeygen.h"
#include "mkTFHEsamples.h"
#include "mkTFHEfunctions.h"



using namespace std;

#define Psize 16 // number of plaintext bits in packed ciphertext
#define lweDim 560
#define num_parties 2

// MK LWE sample Packed (a_1, ..., a_k, b_1, ..., b_l)
struct MKLweSampleP {
	Torus32* a; //-- the parties*n coefs of the mask
    Torus32* b;  //
   	double current_variance; //-- average noise of the sample
   	const int32_t parties;
   	const int32_t n;

#ifdef __cplusplus
   MKLweSampleP(const LweParams* LWEparams, const MKTFHEParams* MKparams);
   ~MKLweSampleP();
   MKLweSampleP(const MKLweSample&)=delete;
   MKLweSampleP& operator=(const MKLweSampleP&)=delete;
#endif
};

MKLweSampleP* new_MKLweSampleP(const LweParams* LWEparams, const MKTFHEParams* MKparams) {
    return new MKLweSampleP(LWEparams, MKparams);
}

// MK LWE sample Packed (a_1, ..., a_k, b_1, ..., b_l)
MKLweSampleP::MKLweSampleP(const LweParams* LWEparams, const MKTFHEParams* MKparams) :
		parties(MKparams->parties), n(LWEparams->n)
{
	this->a = new Torus32[parties*n];
    this->b = new Torus32[Psize];
    this->current_variance = 0.0;
}

MKLweSampleP::~MKLweSampleP() {
    delete[] a;
}

void MKlweCopyP(MKLweSampleP* result, const MKLweSampleP* sample, const MKTFHEParams* params){
    const int32_t n = params->n;
    const int32_t parties = params->parties;

    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result->a[i*n+j] = sample->a[i*n+j];
        }
    }
    
    for (int i=0; i < Psize; i++)
    {
        result->b[i] = sample->b[i];
    }
    

    result->current_variance = sample->current_variance;
}

void MKlweNoiselessTrivialP(MKLweSampleP* result, Torus32 mu, const MKTFHEParams* params){
    const int32_t parties = params->parties;
    const int32_t n = params->n;

    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result->a[i*n +j] = 0;
        }
    }

    for (int l=0; l < Psize; l++)
    {
      result->b[l] = mu;
    }

    result->current_variance = 0.0;
}

void MKlweSubToP(MKLweSampleP* result, const MKLweSampleP* sample, const MKTFHEParams* MKparams){
    const int32_t n = MKparams->n;
    const int32_t parties = MKparams->parties;

    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result->a[i*n+j] -= sample->a[i*n+j];
        }
    }

    for (int l=0; l < Psize; l++)
      result->b[l] -= sample->b[l];

    result->current_variance += sample->current_variance;
}

void MKbootsNANDP(MKLweSample *result, const MKLweSampleP *ca, const MKLweSampleP *cb, 
        const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) {
    
    int parties = result->parties;
    int n = result->n;

    static const Torus32 MU = modSwitchToTorus32(1, 8);
    MKLweSampleP *temp_result = new_MKLweSampleP(LWEparams, MKparams);

    //compute: (0,1/8) - ca - cb
    static const Torus32 NandConst = modSwitchToTorus32(1, 8);
    MKlweNoiselessTrivialP(temp_result, NandConst, MKparams);
    MKlweSubToP(temp_result, ca, MKparams);
    MKlweSubToP(temp_result, cb, MKparams);

    //if the phase is positive, the result is 1/8
    //if the phase is positive, else the result is -1/8

    MKLweSample *temp_result_boots = new_MKLweSample(LWEparams, MKparams);
    MKLweSample *result_boots = new_MKLweSample(LWEparams, MKparams);


    for (int j=0; j < parties * n; j++){
        temp_result_boots->a[j] = temp_result->a[j];
    }

    for (int i=0; i < Psize; i++)
    {
        temp_result_boots->b = temp_result->b[i];
        MKtfhe_bootstrapFFT_v2m2(result_boots, bkFFT, MU, temp_result_boots, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);

        MKlweCopy(&result[i], result_boots, MKparams);

    }

    delete_MKLweSample(temp_result_boots);
    delete_MKLweSample(result_boots);

}

void MkNANDP(MKLweSampleP *result, const MKLweSampleP *ca, const MKLweSampleP *cb, 
        const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) {

    //compute: (0,1/8) - ca - cb
    static const Torus32 NandConst = modSwitchToTorus32(1, 8);
    MKlweNoiselessTrivialP(result, NandConst, MKparams);
    MKlweSubToP(result, ca, MKparams);
    MKlweSubToP(result, cb, MKparams);

}


// **********************************************************************************
// ********************************* MAIN *******************************************
// **********************************************************************************


void dieDramatically(string message) {
    cerr << message << endl;
    abort();
}


void MKbootsAsymEncryptP(MKLweSampleP *result, int32_t msg[Psize], Torus32* crs, Torus32 publicMatrix[Psize])
{

    Torus32 _1s8 = modSwitchToTorus32(1, 8);
    Torus32 mu[Psize];

    double alpha = result->current_variance;
    int parties = result->parties;
    int n = result->n;

    for(int i=0; i < Psize; i++)
        mu[i] = msg[i] ? _1s8 : -_1s8;

    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result->a[i*n +j] = crs[i*n +j];
        }
    }

    for (int l = 0; l < Psize; ++l)
    {
        result->b[l] = gaussian32(mu[l], alpha) + publicMatrix[l];
    }

    result->current_variance = alpha*alpha;

}

int* MKbootsAsymDecryptP(MKLweSampleP *result, Torus32 extSecretMatrix[Psize][num_parties*lweDim])
{

    int parties = result->parties;
    int n = result->n;
    Torus32 ExtA[Psize];
    memset(ExtA, 0, sizeof(ExtA));

    for (int l=0; l < Psize; l++){
        for (int i=0; i < parties * n; i++){
            ExtA[l] += ( extSecretMatrix[l][i] * result->a[i]);
        }
    }

    Torus32 temp[Psize]; // stores b - <a,s>
    static int dec_msg[Psize];

    // decryption of ciphertext
    for (int l = 0; l < Psize; ++l)
    {
        temp[l] = result->b[l] - ExtA[l];
        // cout << temp[l] << " ";
        dec_msg[l] = (temp[l] > 0) ? 1 : 0;

    }

    cout << "\n\n";

    return dec_msg;

}

void MKGenExtCipherP(MKLweSampleP* result, MKLweSampleP* sample, int partyID, Torus32 publicMatrix[num_parties][Psize]){

    double alpha = sample->current_variance;
    int parties = sample->parties;
    int n = sample->n;

    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result->a[i*n +j] = sample->a[i*n +j];
        } 
    }

    for (int l=0; l < Psize; l++)
        result->b[l] = sample->b[l];

    for (int i=0; i < parties; i++)
    {
        if (i != partyID)
        {
            for (int l = 0; l < Psize; ++l)
            {
                result->b[l] += publicMatrix[i][l];
            }
            
        }
    }

    result->current_variance = alpha*alpha;


}


int32_t main(int32_t argc, char **argv) {

    // generate params
    static const int32_t k = 1;
    static const double ks_stdev = 3.05e-5; // 2.44e-5; //standard deviation
    static const double bk_stdev = 3.72e-9; // 3.29e-10; //standard deviation
    static const double max_stdev = 0.012467; // 0.012467; //max standard deviation for a 1/4 msg space
    static const int32_t n = 560; //500;            // LWE modulus
    static const int32_t n_extract = 1024;    // LWE extract modulus (used in bootstrapping)
    static const int32_t hLWE = 0;         // HW secret key LWE --> not used
    static const double stdevLWE = 3.05e-5; // 0.012467;      // LWE ciphertexts standard deviation
    static const int32_t Bksbit = 2;       // Base bit key switching
    static const int32_t dks = 8;          // dimension key switching
    static const double stdevKS = ks_stdev; // 2.44e-5;       // KS key standard deviation
    static const int32_t N = 1024;            // RLWE,RGSW modulus
    static const int32_t hRLWE = 0;        // HW secret key RLWE,RGSW --> not used
    static const double stdevRLWEkey = bk_stdev; // 3.29e-10; // 0; // 0.012467;  // RLWE key standard deviation
    static const double stdevRLWE = bk_stdev; // 3.29e-10; // 0; // 0.012467;     // RLWE ciphertexts standard deviation
    static const double stdevRGSW = bk_stdev; // 3.29e-10;     // RGSW ciphertexts standard deviation
    static const int32_t Bgbit = 8;        // Base bit gadget
    static const int32_t dg = 4;           // dimension gadget
    static const double stdevBK = bk_stdev; // 3.29e-10;       // BK standard deviation
    int32_t parties = num_parties;      // number of parties

    // new parameters
    // 2 parties, B=2^9, d=3 -> works
    // 4 parties, B=2^8, d=4 -> works
    // 8 parties, B=2^6, d=5 -> works



    // params
    LweParams *LWEparams = new_LweParams(n, ks_stdev, max_stdev);
    MKTFHEParams *MKparams = new_MKTFHEParams(n, n_extract, hLWE, stdevLWE, Bksbit, dks, stdevKS, N,
                            hRLWE, stdevRLWEkey, stdevRLWE, stdevRGSW, Bgbit, dg, stdevBK, parties);
    LweParams *extractedLWEparams = new_LweParams(n_extract, ks_stdev, max_stdev);
    TLweParams *RLWEparams = new_TLweParams(N, k, bk_stdev, max_stdev);

    cout << "\nStarting key generation . . .\n" << endl;
    clock_t begin_KG = clock();

    // LWE key
    MKLweKey* MKlwekey = new_MKLweKey(LWEparams, MKparams);
    MKLweKeyGen(MKlwekey);
    cout << "KeyGen MKlwekey: DONE!" << endl;

    // RLWE key
    MKRLweKey* MKrlwekey = new_MKRLweKey(RLWEparams, MKparams);
    MKRLweKeyGen(MKrlwekey);
    cout << "KeyGen MKrlwekey: DONE!" << endl;

    
    // LWE key extracted
    MKLweKey* MKextractedlwekey = new_MKLweKey(extractedLWEparams, MKparams);
    MKtLweExtractKey(MKextractedlwekey, MKrlwekey);
    cout << "KeyGen MKextractedlwekey: DONE!" << endl;

    // bootstrapping + key switching keys
    MKLweBootstrappingKey_v2* MKlweBK = new_MKLweBootstrappingKey_v2(LWEparams, RLWEparams, MKparams);
    MKlweCreateBootstrappingKey_v2(MKlweBK, MKlwekey, MKrlwekey, MKextractedlwekey,
                                extractedLWEparams, LWEparams, RLWEparams, MKparams);
    cout << "KeyGen MKlweBK: DONE!" << endl;

    // bootstrapping FFT + key switching keys
    MKLweBootstrappingKeyFFT_v2* MKlweBK_FFT = new_MKLweBootstrappingKeyFFT_v2(MKlweBK, LWEparams, RLWEparams, MKparams);
    cout << "KeyGen MKlweBK_FFT: DONE!" << endl;

    clock_t end_KG = clock();
    double time_KG = ((double) end_KG - begin_KG)/CLOCKS_PER_SEC;
    cout << "Finished key generation" << endl;
    cout << "Time for key generation: " << time_KG << " seconds" << endl;


    int32_t msg;
    int32_t msg1[Psize];
    int32_t msg2[Psize];

    // use current time as seed for the random generator
    srand(time(0));

    printf("\n1st palintext: ");
    for(int i = 0; i < Psize; i++)
    {
      msg = rand() % 2;
      msg1[i] = msg;
      printf(" %d", msg);
    }

    printf("\n2nd palintext: ");
    for(int i = 0; i < Psize; i++)
    {
      msg = rand() % 2;
      msg2[i] = msg;
      printf(" %d", msg);
    }

    cout << "\nGenerating Public Matrix . . ." << endl;
    auto crs_a = new Torus32[parties*n];
    double alpha = MKlwekey->MKparams->stdevLWE;

    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            crs_a[i*n +j] = uniformTorus32_distrib(generator);
        }
    }

    Torus32 PublicMatrix[parties][Psize];
    for (int i=0; i < parties; i++)
    {
        for (int l=0; l < Psize; l++)
        {
            PublicMatrix[i][l] = gaussian32(0, alpha);
        }
    }


    for (int l=0; l < Psize; l++)
    {
        for (int i=0; i < parties; i++)
        {
            for (int j=0; j < n; j++)
            {
                PublicMatrix[i][l] += ( MKlwekey->key[i].key[j] * crs_a[i*n + j]);
            }
        }
    }

    // generate 2 samples in input
    MKLweSampleP *cipher0 = new_MKLweSampleP(LWEparams, MKparams);
    MKLweSampleP *cipher1 = new_MKLweSampleP(LWEparams, MKparams);

    struct timespec enc_start = {0, 0};
	struct timespec enc_end = {0, 0};

    clock_gettime(CLOCK_MONOTONIC, &enc_start);
    MKbootsAsymEncryptP(cipher0, msg1, crs_a, PublicMatrix[0]);
    clock_gettime(CLOCK_MONOTONIC, &enc_end);

    MKbootsAsymEncryptP(cipher1, msg2, crs_a, PublicMatrix[1]);

    cout << "\nFinished Encryption" << endl;
    cout << "\nEncryption time: " << (((double)enc_end.tv_nsec + 1.0e+9 * enc_end.tv_sec) - ((double)enc_start.tv_nsec + 1.0e+9 * enc_start.tv_sec)) * 1.0e-9 << " sec" << endl;

    MKLweSampleP *ext_cipher0 = new_MKLweSampleP(LWEparams, MKparams);
    MKLweSampleP *ext_cipher1 = new_MKLweSampleP(LWEparams, MKparams);

    MKGenExtCipherP(ext_cipher0, cipher0, 0, PublicMatrix);
    MKGenExtCipherP(ext_cipher1, cipher1, 1, PublicMatrix);

    cout << "\nGenerating Extended Secret Matrix . . ." << endl;
    Torus32 ExtSecretMatrix[Psize][num_parties*n];

    for (int l = 0; l < Psize; ++l)
    {
        for (int i=0; i < parties; i++)
        {
            for (int j = 0; j < n; ++j)
            {
                ExtSecretMatrix[l][i*n+j] = MKlwekey->key[i].key[j];
            }
        }
    }

    struct timespec dec_start = {0, 0};
	struct timespec dec_end = {0, 0};

    clock_gettime(CLOCK_MONOTONIC, &dec_start);
    auto dec_msg0 = MKbootsAsymDecryptP(ext_cipher0, ExtSecretMatrix);
    clock_gettime(CLOCK_MONOTONIC, &dec_end);

    printf("\n1st decrypted msg : ");
    for (int l=0; l < Psize; l++)
    {
        printf(" %d", dec_msg0[l]);
    }

    auto dec_msg1 = MKbootsAsymDecryptP(ext_cipher1, ExtSecretMatrix);
    printf("2nd decrypted msg : ");
    for (int l=0; l < Psize; l++)
    {
        printf(" %d", dec_msg1[l]);
    }

    cout << "\n\nFinished Decryption" << endl;
    cout << "\nDecryption time: " << (((double)dec_end.tv_nsec + 1.0e+9 * dec_end.tv_sec) - ((double)dec_start.tv_nsec + 1.0e+9 * dec_start.tv_sec)) * 1.0e-9 << " sec" << endl;

    // Homomorphic NAND operation
    int outNAND;

    cout << "\nNAND evaluation clear: ";
    for (int l=0; l < Psize; l++)
    {
      outNAND = 1 - (msg1[l] * msg2[l]);
      cout << outNAND << " ";
    }

    cout << "\n";

    // MKLweSample *temp_result = new_MKLweSample_array(Psize, LWEparams, MKparams);
    MKLweSampleP *temp_result = new_MKLweSampleP(LWEparams, MKparams);

    struct timespec eval_start = {0, 0};
	struct timespec eval_end = {0, 0};

    clock_gettime(CLOCK_MONOTONIC, &eval_start);
    MkNANDP(temp_result, ext_cipher0, ext_cipher1, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    clock_gettime(CLOCK_MONOTONIC, &eval_end);

    auto dec_msg_nand = MKbootsAsymDecryptP(temp_result, ExtSecretMatrix);

    cout << "\nDecrypted NAND Msg: ";
    for (int l=0; l < Psize; l++)
    {
        cout << dec_msg_nand[l] << " ";
    }

    cout << "\nFinished Homomrphic NAND Evaluation" << endl;
    cout << "\nNAND Evaluation time: " << (((double)eval_end.tv_nsec + 1.0e+9 * eval_end.tv_sec) - ((double)eval_start.tv_nsec + 1.0e+9 * eval_start.tv_sec)) * 1.0e-9 << " sec" << endl;


    // Threshold Decryption ****************************************************************************

    cout << "\n\nStarting thresholdization on the 1st ciphertext . . ." << endl;

    std::vector<int> subset;

    SecretSharing *ss = new SecretSharing();
    MKKeyShare* share;
    LweKey* tempShare;

    int P = 0;
    cout << "Enter total no. of workers: ";
    cin >> P;

    int th = 0;
    cout << "Enter threshold no. of workers: ";
    cin >> th;

    int tmp;
    cout << "\nEnter party-ids (indices) of participating workers: ";
    for (int i=0; i < th; i++)
    {
        cin >> tmp;
        subset.push_back(tmp);
    }

    // int gropuID = 1;
    int gropuID = findGroupId(subset, th, P);

    Torus32 PartyShareMatrix[P][parties*n];

    struct timespec share_start = {0, 0};
	struct timespec share_end = {0, 0};

    clock_gettime(CLOCK_MONOTONIC, &share_start);
    for (int party=0; party < parties; party++)
    {
      ss->shareSecret(th, P, &(MKlwekey->key[party]), LWEparams);
      share = new MKKeyShare();

      for (int i=0; i < P; i++)
      {
        ss->GetShareSet(i+1, share);
        tempShare = share->GetShare(gropuID);
        for (int j=0; j < n; j++)
        {
          PartyShareMatrix[i][party*n + j] = tempShare->key[j];
        }
      }
    }
    clock_gettime(CLOCK_MONOTONIC, &share_end);
    cout << "\nSecret Sharing time: " << (((double)share_end.tv_nsec + 1.0e+9 * share_end.tv_sec) - ((double)share_start.tv_nsec + 1.0e+9 * share_start.tv_sec)) * 1.0e-9 << " sec" << endl;

    // generate extended secret matrix for partial decryption
    Torus32 ExtSecretMatrixShare[th][Psize][parties*n];

    for (int party=0; party < th; party++)
    {
      for (int l = 0; l < Psize; ++l)
      {
        for (int j=0; j < (parties*n); j++)
        {
          ExtSecretMatrixShare[party][l][j] = PartyShareMatrix[party][j];
        }
      }
    }

    Torus32 PartDecArr[th][Psize];
    Torus32 temp0[Psize]; // stores b - <a,s>
    Torus32 ExtA[Psize];
    memset(ExtA, 0, sizeof(ExtA));

    double sigma = alpha * 2;
    cout << "\nsigma: " << sigma << endl;;

    struct timespec thDec_start = {0, 0};
	struct timespec thDec_end = {0, 0};

    clock_gettime(CLOCK_MONOTONIC, &thDec_start);
    for (int i=0; i < th; i++)
    {
      if (i == 0)
      {
        for (int l=0; l < Psize; l++)
        {
          PartDecArr[i][l] = gaussian32(0, sigma);
          for (int j=0; j < (parties*n); j++)
          {
            PartDecArr[i][l] -= ( ExtSecretMatrixShare[i][l][j] * ext_cipher0->a[j]);
          }
        }
      }
      else
      {
        for (int l=0; l < Psize; l++)
        {
          PartDecArr[i][l] = gaussian32(0, sigma);
          for (int j=0; j < (parties*n); j++)
          {
            PartDecArr[i][l] += ( ExtSecretMatrixShare[i][l][j] * ext_cipher0->a[j]);
          }
        }
      }
    }

    for (int i=0; i < th; i++)
    {
        for (int l=0; l < Psize; l++)
        {
            ExtA[l] += PartDecArr[i][l];
        }
    }

    int dec_msg[Psize];
    for (int l = 0; l < Psize; ++l)
    {
      temp0[l] = ext_cipher0->b[l] - ExtA[l];
      dec_msg[l] = (temp0[l] > 0) ? 1 : 0;

    }

    clock_gettime(CLOCK_MONOTONIC, &thDec_end);
    cout << "\n1st ciphertext -> " << th << "-out-of-" << P << " threshold decryption: ";
    for (int l=0; l < Psize; l++)
    {
      cout << dec_msg[l] << " ";
    }

    cout << "\n";

    cout << "Finished Threshold Decryption" << endl;
    cout << "\nThreshold Decryption time: " << (((double)thDec_end.tv_nsec + 1.0e+9 * thDec_end.tv_sec) - ((double)thDec_start.tv_nsec + 1.0e+9 * thDec_start.tv_sec)) * 1.0e-9 << " sec" << endl;


    while (sigma < 1)
    {
        memset(ExtA, 0, sizeof(ExtA));
        sigma = sigma * 2;
        cout << "\nsigma: " << sigma << " -> ";
        for (int i=0; i < th; i++)
        {
        if (i == 0)
        {
            for (int l=0; l < Psize; l++)
            {
            PartDecArr[i][l] = gaussian32(0, sigma);
            for (int j=0; j < (parties*n); j++)
            {
                PartDecArr[i][l] -= ( ExtSecretMatrixShare[i][l][j] * ext_cipher0->a[j]);
            }
            }
        }
        else
        {
            for (int l=0; l < Psize; l++)
            {
            PartDecArr[i][l] = gaussian32(0, sigma);
            for (int j=0; j < (parties*n); j++)
            {
                PartDecArr[i][l] += ( ExtSecretMatrixShare[i][l][j] * ext_cipher0->a[j]);
            }
            }
        }
        }

        for (int i=0; i < th; i++)
        {
            for (int l=0; l < Psize; l++)
            {
                ExtA[l] += PartDecArr[i][l];
            }
        }

        int dec_msg[Psize];
        for (int l = 0; l < Psize; ++l)
        {
            temp0[l] = ext_cipher0->b[l] - ExtA[l];
            dec_msg[l] = (temp0[l] > 0) ? 1 : 0;
        }

        cout << "1st ciphertext -> " << th << "-out-of-" << P << " threshold decryption: ";
        for (int l=0; l < Psize; l++)
        {
            cout << dec_msg[l] << " ";
        }
    }

    return 0;
}
