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



#include "mkTFHEparams.h"
#include "mkTFHEkeys.h"
#include "mkTFHEkeygen.h"
#include "mkTFHEsamples.h"
#include "mkTFHEfunctions.h"



#define WIDTH 32


using namespace std;



// **********************************************************************************
// ********************************* MAIN *******************************************
// **********************************************************************************


void dieDramatically(string message) {
    cerr << message << endl;
    abort();
}

void MKbootsAsymEncrypt(MKLweSample* result, int32_t message, Torus32* crs, Torus32 pk){

    Torus32 _1s8 = modSwitchToTorus32(1, 8);
    Torus32 mu = message ? _1s8 : -_1s8;

    double alpha = result->current_variance;
    int parties = result->parties;
    int n = result->n;

    result->b = gaussian32(mu, alpha) + pk;

    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result->a[i*n +j] = crs[i*n +j];
        } 
    }
    
    result->current_variance = alpha*alpha;
}

void genExtCipher(MKLweSample* result, MKLweSample* sample, int partyID, vector<Torus32> publicKey){

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

    result->b = sample->b;

    for (int i=0; i < parties; i++)
    {
        if (i != partyID)
        {
            result->b += publicKey[i];
        }
    }

    result->current_variance = alpha*alpha;


}

void genExtCipherArr(MKLweSample* result, MKLweSample* sample, int partyID, vector<Torus32> publicKey)
{
    for (int i = 0; i < WIDTH; i++){
        genExtCipher(&result[i], &sample[i], partyID, publicKey);
    }
}


void MkNAND(MKLweSample *result, const MKLweSample *ca, const MKLweSample *cb, 
        const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) {

    //compute: (0,1/8) - ca - cb
    static const Torus32 NandConst = modSwitchToTorus32(1, 8);
    MKlweNoiselessTrivial(result, NandConst, MKparams);
    MKlweSubTo(result, ca, MKparams);
    MKlweSubTo(result, cb, MKparams);

}

void MKbootsNAND(MKLweSample *result, const MKLweSample *ca, const MKLweSample *cb, 
        const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) {

    static const Torus32 MU = modSwitchToTorus32(1, 8);

    MKLweSample *temp_result = new_MKLweSample(LWEparams, MKparams);

    //compute: (0,1/8) - ca - cb
    static const Torus32 NandConst = modSwitchToTorus32(1, 8);
    MKlweNoiselessTrivial(temp_result, NandConst, MKparams);
    MKlweSubTo(temp_result, ca, MKparams);
    MKlweSubTo(temp_result, cb, MKparams);

    //if the phase is positive, the result is 1/8
    //if the phase is positive, else the result is -1/8
    MKtfhe_bootstrapFFT_v2m2(result, bkFFT, MU, temp_result, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);

    delete_MKLweSample(temp_result);
}

int32_t main(int32_t argc, char **argv) {


    // generate params 
    static const int32_t k = 1;
    static const double ks_stdev = 3.05e-5;// 2.44e-5; //standard deviation
    static const double bk_stdev = 3.72e-9; // 3.29e-10; //standard deviation
    static const double max_stdev = 0.012467; //max standard deviation for a 1/4 msg space
    static const int32_t n = 560; //500;            // LWE modulus
    static const int32_t n_extract = 1024;    // LWE extract modulus (used in bootstrapping)
    static const int32_t hLWE = 0;         // HW secret key LWE --> not used
    static const double stdevLWE = 0.012467;      // LWE ciphertexts standard deviation
    static const int32_t Bksbit = 2;       // Base bit key switching
    static const int32_t dks = 8;          // dimension key switching
    static const double stdevKS = ks_stdev; // 2.44e-5;       // KS key standard deviation
    static const int32_t N = 1024;            // RLWE,RGSW modulus
    static const int32_t hRLWE = 0;        // HW secret key RLWE,RGSW --> not used
    static const double stdevRLWEkey = bk_stdev; // 3.29e-10; // 0; // 0.012467;  // RLWE key standard deviation
    static const double stdevRLWE = bk_stdev; // 3.29e-10; // 0; // 0.012467;     // RLWE ciphertexts standard deviation
    static const double stdevRGSW = bk_stdev; // 3.29e-10;     // RGSW ciphertexts standard deviation 
    static const int32_t Bgbit = 9;        // Base bit gadget
    static const int32_t dg = 3;           // dimension gadget
    static const double stdevBK = bk_stdev; // 3.29e-10;       // BK standard deviation
    static const int32_t parties = 2;      // number of parties

    // new parameters 
    // 2 parties, B=2^9, d=3 -> works
    // 4 parties, B=2^8, d=4 -> works
    // 8 parties, B=2^6, d=5 -> works 
    

    // params
    LweParams *extractedLWEparams = new_LweParams(n_extract, ks_stdev, max_stdev);
    LweParams *LWEparams = new_LweParams(n, ks_stdev, max_stdev);
    TLweParams *RLWEparams = new_TLweParams(N, k, bk_stdev, max_stdev);
    MKTFHEParams *MKparams = new_MKTFHEParams(n, n_extract, hLWE, stdevLWE, Bksbit, dks, stdevKS, N, 
                            hRLWE, stdevRLWEkey, stdevRLWE, stdevRGSW, Bgbit, dg, stdevBK, parties);


    cout << "Params: DONE!" << endl;


   
    // Key generation 
    cout << "Starting KEY GENERATION" << endl;
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
    cout << "Finished KEY GENERATION" << endl;
    cout << "\nTime per KEY GENERATION " << time_KG << " sec" << endl;

    cout << "\nGenerating Public Keys . . .\n" << endl;

    vector<Torus32> PublicKey(parties);
    auto crs_a = new Torus32[parties*n];
    double alpha = MKlwekey->MKparams->stdevLWE;

    for (int i=0; i < parties; i++)
    {
        Torus32 PublicKeyTmp = gaussian32(0, alpha);
        for (int j = 0; j < n; ++j)
        {
            crs_a[i*n +j] = uniformTorus32_distrib(generator);
            PublicKeyTmp += crs_a[i*n +j]*MKlwekey->key[i].key[j];
        } 

        int choice = rand() % 10;
        // cout << "\nChoice: " << choice << endl;

        for (int k=0; k < choice; k++)
        {
            PublicKey[i] += PublicKeyTmp;
        }

    }
    
    int32_t msg1[WIDTH];
    int32_t msg2[WIDTH];
    int32_t outNAND[WIDTH];

    // use current time as seed for the random generator
    srand(time(0));

    printf("\n1st palintext: ");
    for(int i = 0; i < WIDTH; i++)
    {
      msg1[i] = rand() % 2;
      printf(" %d", msg1[i]);
    }

    printf("\n2nd palintext: ");
    for(int i = 0; i < WIDTH; i++)
    {
      msg2[i] = rand() % 2;
      printf(" %d", msg2[i]);
    }

    for (int l=0; l < WIDTH; l++)
    {
        outNAND[l] = 1 - (msg1[l] * msg2[l]);
    }


    // generate 2 samples in input
    MKLweSample *test_in1 = new_MKLweSample_array(WIDTH, LWEparams, MKparams);
    MKLweSample *test_in2 = new_MKLweSample_array(WIDTH, LWEparams, MKparams);

    for (int i=0; i < WIDTH; i++)
    {
        MKbootsAsymEncrypt(&test_in1[i], msg1[i], crs_a, PublicKey[0]);
    }

    for (int i=0; i < WIDTH; i++)
    {
        MKbootsAsymEncrypt(&test_in2[i], msg2[i], crs_a, PublicKey[1]);
    }

    cout << "\nEncryption: DONE!" << endl;

    cout << "\nGenerating Extended Ciphertext . . . \n\n" << endl;

    MKLweSample *ext_test_in1 = new_MKLweSample_array(WIDTH, LWEparams, MKparams);
    MKLweSample *ext_test_in2 = new_MKLweSample_array(WIDTH, LWEparams, MKparams);

    genExtCipherArr(ext_test_in1, test_in1, 0, PublicKey);
    genExtCipherArr(ext_test_in2, test_in2, 1, PublicKey);

    // verify encrypt 
    int32_t mess1_dec[WIDTH];
    int32_t mess2_dec[WIDTH];

    for (int i=0; i < WIDTH; i++)
    {
        mess1_dec[i] = MKbootsSymDecrypt(&ext_test_in1[i], MKlwekey);
    }

    for (int i=0; i < WIDTH; i++)
    {
        mess2_dec[i] = MKbootsSymDecrypt(&ext_test_in2[i], MKlwekey);
    }

    printf("1st Decrypted msg : ");
    for (int l=0; l < WIDTH; l++)
    {
        cout << mess1_dec[l] << " ";
    }
    cout << "\n";

    printf("2nd Decrypted msg : ");
    for (int l=0; l < WIDTH; l++)
    {
        cout << mess2_dec[l] << " ";
    }
    cout << "\n";

    // generate output sample
    MKLweSample *test_out_v2m2 = new_MKLweSample_array(WIDTH, LWEparams, MKparams);

    // evaluate MK bootstrapped NAND 

    struct timespec eval_start = {0, 0};
	struct timespec eval_end = {0, 0};

    clock_gettime(CLOCK_MONOTONIC, &eval_start);
    for (int i=0; i < WIDTH; i++)
    {
        MkNAND(&test_out_v2m2[i], &ext_test_in1[i], &ext_test_in2[i], MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    }
    clock_gettime(CLOCK_MONOTONIC, &eval_end);


    cout << "\n\nFinished MK NAND Gate" << endl;
    cout << "\nNAND Evaluation time: " << (((double)eval_end.tv_nsec + 1.0e+9 * eval_end.tv_sec) - ((double)eval_start.tv_nsec + 1.0e+9 * eval_start.tv_sec)) * 1.0e-9 << " sec" << endl;


    // verify NAND
    cout << "\n\nNAND evaluation clear: ";
    for (int l=0; l < WIDTH; l++)
    {
        cout << outNAND[l] << " ";
    }
    cout << "\n";

    int32_t outNAND_v2m2[WIDTH];
    for (int i=0; i < WIDTH; i++)
    {
        outNAND_v2m2[i] = MKbootsSymDecrypt(&test_out_v2m2[i], MKlwekey);
    }

    cout << "Decrypted NAND evaluation: ";
    for (int l=0; l < WIDTH; l++)
    {
        cout << outNAND_v2m2[l] << " ";
    }
    cout << "\n";


    // delete samples
    delete_MKLweSample(test_out_v2m2);
    delete_MKLweSample(test_in2);
    delete_MKLweSample(test_in1);

    // delete keys
    delete_MKLweBootstrappingKeyFFT_v2(MKlweBK_FFT);
    delete_MKLweBootstrappingKey_v2(MKlweBK);
    delete_MKLweKey(MKextractedlwekey);
    delete_MKRLweKey(MKrlwekey);
    delete_MKLweKey(MKlwekey);
    // delete params
    delete_MKTFHEParams(MKparams);
    delete_TLweParams(RLWEparams);
    delete_LweParams(LWEparams);
    delete_LweParams(extractedLWEparams);


    return 0;
}
