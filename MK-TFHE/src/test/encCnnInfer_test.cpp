#include <iostream>
#include <vector>
#include <cstdio>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <cassert>
#include <functional>
#include <stdio.h>
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
#include <math.h>
#include <string.h>

using namespace std;
#define WIDTH 9

const int INPUT_WIDTH = 28;
const int INPUT_HEIGHT = 28;
const int NUM_CLASSES = 10;
const int WINDOW_SIZE = 4;
const int FILTER_SIZE = 4;
const int NUM_FILTERS = 5;
const int STRIDE = 2;
const int NUM_TEST_IMAGES = 10000;
const int NUM_TEST_LABELS = 10000;

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

void IntAsymEncrypt(MKLweSample *result, int p, Torus32* crs, Torus32 pk)
{
    for (int i = 0; i < WIDTH; i++){
        int bit = p & 1;
        p = (p >> 1);
        MKbootsAsymEncrypt(&result[i], bit, crs, pk);
    }
}

int IntSymDecrypt(MKLweSample *c, MKLweKey* MKlwekey)
{
    int msb = MKbootsSymDecrypt(&c[WIDTH-1], MKlwekey);

    int ans = 0;
    for (int i = 0; i < WIDTH - 1; i++){
        int bit = MKbootsSymDecrypt(&c[i], MKlwekey);
        ans += (bit ^ msb) << i;
    }

    if (msb == 1){
        ans++;
        ans = -ans;
    }

    return ans;
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

// circuit

/** result = result + sample */
void MKlweAddTo(MKLweSample* result, const MKLweSample* sample, const MKTFHEParams* MKparams){
    const int32_t n = MKparams->n;
    const int32_t parties = MKparams->parties;

    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result->a[i*n+j] += sample->a[i*n+j];
        }
    }

    result->b += sample->b;
    result->current_variance += sample->current_variance; 
}

/** result = result + p.sample */
void MKlweAddMulTo(MKLweSample* result, int32_t p, const MKLweSample* sample, const MKTFHEParams* MKparams){
    const int32_t n = MKparams->n;
    const int32_t parties = MKparams->parties;

    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result->a[i*n+j] += p*sample->a[i*n+j];
        }
    }

    result->b += p*sample->b;
    result->current_variance += (p*p)*sample->current_variance; 
}

void MKbootsMUX(MKLweSample *result, const MKLweSample *a, const MKLweSample *b, const MKLweSample *c,
                     const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) {

    static const Torus32 MU = modSwitchToTorus32(1, 8);

    MKLweSample *temp_result = new_MKLweSample(LWEparams, MKparams);
    MKLweSample *temp_result1 = new_MKLweSample(extractedLWEparams, MKparams);
    MKLweSample *u1 = new_MKLweSample(extractedLWEparams, MKparams);
    MKLweSample *u2 = new_MKLweSample(extractedLWEparams, MKparams);


    //compute "AND(a,b)": (0,-1/8) + a + b
    static const Torus32 AndConst = modSwitchToTorus32(-1, 8);
    MKlweNoiselessTrivial(temp_result, AndConst, MKparams);
    MKlweAddTo(temp_result, a, MKparams);
    MKlweAddTo(temp_result, b, MKparams);

    MKtfhe_bootstrapFFT_v2m2(u1, bkFFT, MU, temp_result, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey); 


    //compute "AND(not(a),c)": (0,-1/8) - a + c
    MKlweNoiselessTrivial(temp_result, AndConst, MKparams);
    MKlweSubTo(temp_result, a, MKparams);
    MKlweAddTo(temp_result, c, MKparams);

    MKtfhe_bootstrapFFT_v2m2(u2, bkFFT, MU, temp_result, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);   

    // Add u1=u1+u2
    static const Torus32 MuxConst = modSwitchToTorus32(1, 8);
    MKlweNoiselessTrivial(temp_result1, MuxConst, MKparams);
    MKlweAddTo(temp_result1, u1, MKparams);
    MKlweAddTo(temp_result1, u2, MKparams);

    MKtfhe_bootstrapFFT_v2m2(result, bkFFT, MU, temp_result1, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);


    delete_MKLweSample(u2);
    delete_MKLweSample(u1);
    delete_MKLweSample(temp_result1);
    delete_MKLweSample(temp_result);
}


void MKbootsOR(MKLweSample *result, const MKLweSample *ca, const MKLweSample *cb, 
        const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) {

    static const Torus32 MU = modSwitchToTorus32(1, 8);

    MKLweSample *temp_result = new_MKLweSample(LWEparams, MKparams);

    //compute: (0,1/8) + ca + cb
    static const Torus32 OrConst = modSwitchToTorus32(1, 8);
    MKlweNoiselessTrivial(temp_result, OrConst, MKparams);
    MKlweAddTo(temp_result, ca, MKparams);
    MKlweAddTo(temp_result, cb, MKparams);

    //if the phase is positive, the result is 1/8
    //if the phase is positive, else the result is -1/8
    MKtfhe_bootstrapFFT_v2m2(result, bkFFT, MU, temp_result, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);

    delete_MKLweSample(temp_result);
}


void MKbootsXOR(MKLweSample *result, const MKLweSample *ca, const MKLweSample *cb, 
        const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) {

    static const Torus32 MU = modSwitchToTorus32(1, 8);

    MKLweSample *temp_result = new_MKLweSample(LWEparams, MKparams);

    //compute: (0,1/4) + 2*(ca + cb)
    static const Torus32 XorConst = modSwitchToTorus32(1, 4);
    MKlweNoiselessTrivial(temp_result, XorConst, MKparams);
    MKlweAddMulTo(temp_result, 2, ca, MKparams);
    MKlweAddMulTo(temp_result, 2, cb, MKparams);

    //if the phase is positive, the result is 1/8
    //if the phase is positive, else the result is -1/8
    MKtfhe_bootstrapFFT_v2m2(result, bkFFT, MU, temp_result, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);

    delete_MKLweSample(temp_result);
}


void MKbootsAND(MKLweSample *result, const MKLweSample *ca, const MKLweSample *cb, 
        const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) {

    static const Torus32 MU = modSwitchToTorus32(1, 8);

    MKLweSample *temp_result = new_MKLweSample(LWEparams, MKparams);

    //compute: (0,-1/8) + ca + cb
    static const Torus32 AndConst = modSwitchToTorus32(-1, 8);
    MKlweNoiselessTrivial(temp_result, AndConst, MKparams);
    MKlweAddTo(temp_result, ca, MKparams);
    MKlweAddTo(temp_result, cb, MKparams);

    //if the phase is positive, the result is 1/8
    //if the phase is positive, else the result is -1/8
    MKtfhe_bootstrapFFT_v2m2(result, bkFFT, MU, temp_result, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);

    delete_MKLweSample(temp_result);
}


void CktAddWithCarry(const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey,
            MKLweSample *res, MKLweSample *a, MKLweSample *b, MKLweSample *Cin = NULL)
{
    auto carry = new_MKLweSample(LWEparams, MKparams);
    auto tmp1 = new_MKLweSample(LWEparams, MKparams);
    auto tmp2 = new_MKLweSample(LWEparams, MKparams);
    auto tmp3 = new_MKLweSample(LWEparams, MKparams);

    auto _cin = Cin;

    for (int i = 0; i < WIDTH; i++){

        MKbootsXOR(tmp1, &a[i], &b[i], bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);   // tmp1 = a^b
        MKbootsAND(tmp2, &a[i], &b[i], bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);   // tmp2 = ab
        if (_cin){
            MKbootsXOR(&res[i], tmp1, _cin, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);  // res = Cin^a^b
            MKbootsAND(tmp3, tmp1, _cin, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
            MKbootsOR(carry, tmp2, tmp3, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);    // carry = ab + Cin(a^b)  
        }else{
            MKlweCopy(&res[i], tmp1, MKparams);
            MKlweCopy(carry, tmp2, MKparams);
        }

        _cin = carry;

    }

    MKlweCopy(&res[WIDTH], _cin, MKparams);

    delete_MKLweSample(carry);
    delete_MKLweSample(tmp1);
    delete_MKLweSample(tmp2);
    delete_MKLweSample(tmp3);
}


void CktAdd(const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey,
            MKLweSample *res, MKLweSample *a, MKLweSample *b, MKLweSample *Cin = NULL)
{
    auto carry = new_MKLweSample(LWEparams, MKparams);
    auto tmp1 = new_MKLweSample(LWEparams, MKparams);
    auto tmp2 = new_MKLweSample(LWEparams, MKparams);
    auto tmp3 = new_MKLweSample(LWEparams, MKparams);

    auto _cin = Cin;

    for (int i = 0; i < WIDTH; i++){

        MKbootsXOR(tmp1, &a[i], &b[i], bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);   // tmp1 = a^b
        MKbootsAND(tmp2, &a[i], &b[i], bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);   // tmp2 = ab
        if (_cin){
            MKbootsXOR(&res[i], tmp1, _cin, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);  // res = Cin^a^b
            MKbootsAND(tmp3, tmp1, _cin, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
            MKbootsOR(carry, tmp2, tmp3, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);    // carry = ab + Cin(a^b)  
        }else{
            MKlweCopy(&res[i], tmp1, MKparams);
            MKlweCopy(carry, tmp2, MKparams);
        }

        _cin = carry;

    }

    delete_MKLweSample(carry);
    delete_MKLweSample(tmp1);
    delete_MKLweSample(tmp2);
    delete_MKLweSample(tmp3);
}


void CktInv(const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey,
            MKLweSample *res, MKLweSample *a, MKLweSample *one)
{
    for (int i = 0; i < WIDTH; i++){
        MKbootsXOR(&res[i], &a[i], one, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    }
}


void CktSub(const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey,
            MKLweSample *res, MKLweSample *a, MKLweSample *b, MKLweSample *one)
{
    auto tmp = new_MKLweSample_array(WIDTH, LWEparams, MKparams);

    CktInv(bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, tmp, b, one);


    CktAdd(bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, res, a, tmp, one);
    delete_MKLweSample_array(WIDTH, tmp);
   
}


void nBitMul(const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey,
            MKLweSample *res, MKLweSample *a, MKLweSample *b, MKLweSample *zero, MKLweSample *Cin = NULL)
{
    MKLweSample** BArr = new MKLweSample*[WIDTH];

    auto result = new_MKLweSample_array(2*WIDTH, LWEparams, MKparams);

    for (int i = 0; i < WIDTH; i++) {
        BArr[i] = new_MKLweSample_array(WIDTH, LWEparams, MKparams);
    }

    for (int i=0; i < WIDTH; i++)
    {
        for (int j=0; j < WIDTH; j++)
        {
            MKbootsAND(&BArr[i][j], &a[j], &b[i], bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
        }
    }

    MKlweCopy(&result[0], &BArr[0][0], MKparams);

    auto tmpIn = new_MKLweSample_array(WIDTH, LWEparams, MKparams);

    for (int i=0; i < WIDTH-1; i++)
        MKlweCopy(&tmpIn[i], &BArr[0][i+1], MKparams);

    MKlweCopy(&tmpIn[WIDTH-1], zero, MKparams);

    int ctr=0;

    for (int i=1; i < WIDTH-1; i++)
    {
        auto tmpArr = new_MKLweSample_array(WIDTH+1, LWEparams, MKparams);

        CktAddWithCarry(bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, tmpArr, tmpIn, BArr[i]);

        MKlweCopy(&result[i], &tmpArr[0], MKparams);

        for (int j=0; j < WIDTH; j++)
            MKlweCopy(&tmpIn[j], &tmpArr[j+1], MKparams);

        ctr = i;

    }

    ctr++;

    auto tmpArr = new_MKLweSample_array(WIDTH+1, LWEparams, MKparams);

    CktAddWithCarry(bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, tmpArr, tmpIn, BArr[ctr]);

    for (int i=0; i <= WIDTH; i++)
        MKlweCopy(&result[i+ctr], &tmpArr[i], MKparams);

    // copy WIDTH LSB bits
    for (int i=0; i < WIDTH; i++)
        MKlweCopy(&res[i], &result[i], MKparams);


}


void CktMul(const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey,
            MKLweSample *res, MKLweSample *a, MKLweSample *b)
{

    auto zero = new_MKLweSample(LWEparams, MKparams);
    auto mres = new_MKLweSample(LWEparams, MKparams);

    MKbootsXOR(mres, &a[WIDTH-1], &b[WIDTH-1], bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);

    MKbootsMUX(&a[WIDTH-1], &a[WIDTH-1], zero, &a[WIDTH-1], bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    MKbootsMUX(&b[WIDTH-1], &b[WIDTH-1], zero, &b[WIDTH-1], bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);

    nBitMul(bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, res, a, b, zero);

    MKbootsMUX(&res[WIDTH-1], mres, mres, mres, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);

}


void CktLess(const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey,
            MKLweSample *res, MKLweSample *a, MKLweSample *b, MKLweSample *one)     // res = (a < b)
{
    auto tmp = new_MKLweSample_array(WIDTH, LWEparams, MKparams);

    CktSub(bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, tmp, a, b, one);

    MKlweCopy(res, &tmp[WIDTH-1], MKparams);
    delete_MKLweSample_array(WIDTH, tmp);
}


void CktGeq(const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey,
            MKLweSample *res, MKLweSample *a, MKLweSample *b, MKLweSample *one)     // res = (a >= b)
{
    auto tmp = new_MKLweSample(LWEparams, MKparams);
    CktLess(bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, tmp, a, b, one);
    MKbootsXOR(res, tmp, one, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    delete_MKLweSample(tmp);
}


void CktMax(const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey,
            MKLweSample *res, MKLweSample *a, MKLweSample *b, MKLweSample *one)     // res = max(a, b)
{
    auto tmpRes = new_MKLweSample_array(WIDTH, LWEparams, MKparams);
    CktGeq(bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, tmpRes, a, b, one);
    MKbootsMUX(res, tmpRes, a, b, bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
    
}


vector<vector<vector<MKLweSample*>>> enc_conv2d(vector<vector<MKLweSample*>>& input, vector<vector<vector<MKLweSample*>>>& kernels, vector<MKLweSample*> &bias, 
                        const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
                        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey, 
                        Torus32* crs, vector<Torus32> publicKey, int id) {

    int input_width = input.size();
    int input_height = input[0].size();
    int output_height = (input_height - FILTER_SIZE) / STRIDE + 1;
    int output_width = (input_width - FILTER_SIZE) / STRIDE + 1;

    vector<vector<vector<MKLweSample*>>> outputs(NUM_FILTERS, vector<vector<MKLweSample*>>(output_width, vector<MKLweSample*>(output_height)));

    for (int c = 0; c < NUM_FILTERS; c++)
    {   
        cout << "c: " << c << endl;
        for (int i = 0; i < output_height; i++) {
            for (int j = 0; j < output_width; j++) {
                auto sum = new_MKLweSample_array(WIDTH, LWEparams, MKparams);
                for (int m = 0; m < FILTER_SIZE; m++) {
                    for (int n = 0; n < FILTER_SIZE; n++) {
                        int x = i * STRIDE + m;
                        int y = j * STRIDE + n;
                        cout << "x: " << x << ", y: " << y << endl;
                        // sum += input[x][y] * kernels[c][m][n];
                        auto mulTmp = new_MKLweSample_array(WIDTH, LWEparams, MKparams);
                        CktMul(bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, mulTmp, input[x][y], kernels[c][m][n]);
                        CktAdd(bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, sum, sum, mulTmp);
                    }

                }
                // outputs[c][i][j] = sum + bias[c];
                CktAdd(bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, outputs[c][i][j], sum, bias[c]);
            }

        }
    }

    return outputs;
}


vector<vector<vector<MKLweSample*>>> enc_max_pooling2d(vector<vector<vector<MKLweSample*>>> &input,
                        const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
                        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey, 
                        Torus32* crs, vector<Torus32> publicKey, int id) {

    int input_height = input[0].size();
    int input_width = input[0][0].size();
    int output_height = input_height / WINDOW_SIZE;
    int output_width = input_width / WINDOW_SIZE;

    auto one = new_MKLweSample_array(WIDTH, LWEparams, MKparams);
    IntAsymEncrypt(one, 1, crs, publicKey[id]);
    genExtCipherArr(one, one, id, publicKey);

    auto max_val = new_MKLweSample_array(WIDTH, LWEparams, MKparams);
    IntAsymEncrypt(max_val, -10, crs, publicKey[id]);
    genExtCipherArr(max_val, max_val, id, publicKey);

    vector<vector<vector<MKLweSample*>>> outputs(NUM_FILTERS, vector<vector<MKLweSample*>>(output_height, vector<MKLweSample*>(output_width)));

    for (int c = 0; c < NUM_FILTERS; c++){
        for (int i = 0; i < output_height; i++) {
            for (int j = 0; j < output_width; j++) {
                auto enc_max_val = max_val;
                for (int m = 0; m < WINDOW_SIZE; m++) {
                    for (int n = 0; n < WINDOW_SIZE; n++) {
                        int x = i * WINDOW_SIZE + m;
                        int y = j * WINDOW_SIZE + n;
                        // max_val = max(max_val, input[c][x][y]);
                        CktMax(bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, enc_max_val, enc_max_val, input[c][x][y], one);
                    }
                }
                outputs[c][i][j] = enc_max_val;
            }
        }
    }
    
    return outputs;
}


vector<MKLweSample*> enc_flatten(vector<vector<vector<MKLweSample*>>> &inputs)
{
    vector<MKLweSample*> outputs;

    for (const auto &channel : inputs) {
        for (const auto &row : channel) {
            for (const auto &val : row) {
                outputs.push_back(val);
            }
        }
    }

    return outputs;
}


vector<MKLweSample*> enc_fc_layer(vector<MKLweSample*>& input, vector<vector<MKLweSample*>>& weights, vector<MKLweSample*>& biases, 
                            const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
                            const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey, 
                            Torus32* crs, vector<Torus32> publicKey, int id) {

    int input_size = input.size();
    int output_size = biases.size();

    vector<MKLweSample*> output(output_size);

    for (int i = 0; i < output_size; i++) {
        for (int j = 0; j < input_size; j++) {
        //   output[i] += input[j] * weights[j][i];
            auto mulTmp = new_MKLweSample_array(WIDTH, LWEparams, MKparams);
            CktMul(bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, mulTmp, input[j], weights[j][i]);
            CktAdd(bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, output[i], output[i], mulTmp);
        }
        // output[i] += biases[i];
        CktAdd(bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, output[i], output[i], biases[i]);
    }

    return output;
}


MKLweSample* CktRelu(MKLweSample* x, const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
                            const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey, 
                            Torus32* crs, vector<Torus32> publicKey, int id){

    auto res = new_MKLweSample_array(WIDTH, LWEparams, MKparams);
    auto zero = new_MKLweSample_array(WIDTH, LWEparams, MKparams);
    auto one = new_MKLweSample_array(WIDTH, LWEparams, MKparams);

    IntAsymEncrypt(zero, 0, crs, publicKey[id]);
    genExtCipherArr(zero, zero, id, publicKey);

    IntAsymEncrypt(one, 1, crs, publicKey[id]);
    genExtCipherArr(one, one, id, publicKey);

    CktMax(bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, res, zero, x, one);
    
    return res;
}


vector<vector<vector<MKLweSample*>>> enc_relu(vector<vector<vector<MKLweSample*>>> &inputs, const MKLweBootstrappingKeyFFT_v2 *bkFFT, 
                                const LweParams* LWEparams, const LweParams *extractedLWEparams, const TLweParams* RLWEparams, 
                                const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey, 
                                Torus32* crs, vector<Torus32> publicKey, int id){
    int depth = inputs.size();
    int height = inputs[0].size();
    int width = inputs[0][0].size();

    vector<vector<vector<MKLweSample*>>> outputs(depth, vector<vector<MKLweSample*>>(height, vector<MKLweSample*>(width)));

    for (int d = 0; d < depth; d++)
    {
        for (int i = 0; i < height; i++)
        {
            for (int j=0; j < width; j++)
            {
                outputs[d][i][j] = CktRelu(inputs[d][i][j], bkFFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, crs, publicKey, id);
            }
        }
    }

        return outputs;
}

vector<vector<int>> transpose_matrix(vector<vector<int>> original_vector){

    vector<vector<int>> transposed_vector((int)original_vector[0].size(), vector<int>((int)original_vector.size()));

    for (int i = 0; i < (int)original_vector.size(); ++i) {
        for (int j = 0; j < (int)original_vector[0].size(); ++j) {
            transposed_vector[j][i] = original_vector[i][j];
        }
    }
    
    return transposed_vector;
}

vector<vector<vector<int>>> load_test_data()
{
    ifstream file_test_data("test/mnist_data/t10k-images.idx3-ubyte", std::ios::binary);

    if (!file_test_data.is_open()) {
        std::cerr << "Failed to open image file\n";
    }

    // Skip the first 16 bytes of the file, which contain header information
    file_test_data.seekg(16);

    // Read the contents of the file into a vector
    std::vector<unsigned char> train_data((std::istreambuf_iterator<char>(file_test_data)),
                                    std::istreambuf_iterator<char>());

    // Create a 3D array to store the image data
    std::vector<std::vector<std::vector<int>>> images(NUM_TEST_IMAGES, std::vector<std::vector<int>>(INPUT_HEIGHT, std::vector<int>(INPUT_WIDTH)));

    // Fill the image data from the vector
    for (int i = 0; i < NUM_TEST_IMAGES; i++) {
        for (int j = 0; j < INPUT_HEIGHT; j++) {
            for (int k = 0; k < INPUT_WIDTH; k++) {
                images[i][j][k] = ((int)train_data[i * INPUT_HEIGHT * INPUT_WIDTH + j * INPUT_WIDTH + k] > 0) ? 1 : 0;
            }
        }
    }

    return images;
}

std::vector<int> load_test_labels()
{
    std::ifstream file_test_label("test/mnist_data/t10k-labels.idx1-ubyte", std::ios::binary);

    if (!file_test_label.is_open()) {
        std::cerr << "Failed to open label file\n";
    }

    // Skip the first 8 bytes of the file, which contain header information
    file_test_label.seekg(8);

    // Read the contents of the file into a vector
    std::vector<unsigned char> train_label((std::istreambuf_iterator<char>(file_test_label)),
                                    std::istreambuf_iterator<char>());

    // Create a 1D array to store the label data
    std::vector<int> labels(NUM_TEST_LABELS);

    // Fill the label data from the vector
    for (int i = 0; i < NUM_TEST_LABELS; i++) {
        labels[i] = (int)train_label[i];
    }

    return labels;
}

vector<vector<vector<int>>> load_cnn_weight()
{
    ifstream file("test/model_params/cnn_weight.txt");
    if (!file.is_open()) {
        cerr << "Failed to open the file" << endl;
    }

    int rows, columns, depth;
    file >> rows >> columns >> depth;

    // Create a 3D vector to store the array
    vector<vector<vector<int>>> arr(depth, vector<vector<int>>(rows, vector<int>(columns)));

    // Read the data from the file and store it in the 3D vector
    for (int d = 0; d < depth; ++d) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                file >> arr[d][i][j];
            }
        }
    }

    // Close the file
    file.close();

    return arr;
}

vector<int> load_cnn_bias()
{
    ifstream file("test/model_params/cnn_bias.txt");
    if (!file.is_open()) {
        cerr << "Failed to open the file" << endl;
    }

    int columns;
    file >> columns;

    // Create a 1D vector to store the array
    vector<int> arr(columns);

    // Read the data from the file and store it in the 1D vector
    for (int j = 0; j < columns; ++j) {
        file >> arr[j];
    }

    // Close the file
    file.close();

    return arr;
}

vector<vector<int>> load_fc_weight()
{
    ifstream file("test/model_params/fc_weight.txt");
    if (!file.is_open()) {
        cerr << "Failed to open the file" << endl;
    }

    int rows, columns;
    file >> rows >> columns;

    // Create a 2D vector to store the array
    vector<vector<int>> arr(rows, vector<int>(columns));

    // Read the data from the file and store it in the 2D vector
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            file >> arr[i][j];
        }
    }

    // Close the file
    file.close();

    return arr;
}

vector<int> load_fc_bias()
{
    ifstream file("test/model_params/fc_bias.txt");
    if (!file.is_open()) {
        cerr << "Failed to open the file" << endl;
    }

    int columns;
    file >> columns;

    // Create a 1D vector to store the array
    vector<int> arr(columns);

    // Read the data from the file and store it in the 1D vector
    for (int j = 0; j < columns; ++j) {
        file >> arr[j];
    }

    // Close the file
    file.close();

    return arr;
}



int32_t main(int32_t argc, char **argv) {

    // generate params
    static const int32_t k = 1;
    static const double ks_stdev = 2.44e-5; //3.05e-5;// 2.44e-5; //standard deviation
    static const double bk_stdev = 3.29e-10; // 3.72e-9; // 3.29e-10; //standard deviation
    static const double max_stdev = 0.012467; //max standard deviation for a 1/4 msg space 0.012467
    static const int32_t n = 560; //500;            // LWE modulus
    static const int32_t n_extract = 1024;    // LWE extract modulus (used in bootstrapping)
    static const int32_t hLWE = 0;         // HW secret key LWE --> not used
    static const double stdevLWE = 0.012467;      // LWE ciphertexts standard deviation 0.012467
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
    static const int32_t parties = 2;      // number of parties

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

    auto test_images = load_test_data();
    auto test_labels = load_test_labels();

    cout << "Data loaded." << endl;

    auto cnn_kernel = load_cnn_weight();
    auto cnn_bias = load_cnn_bias();
    auto fc_weight = load_fc_weight();
    auto fc_weightT = transpose_matrix(fc_weight);
    auto fc_bias = load_fc_bias();

    cout << "Params loaded." << endl;

    auto zero = new_MKLweSample_array(WIDTH, LWEparams, MKparams);
    IntAsymEncrypt(zero, 0, crs_a, PublicKey[0]);
    genExtCipherArr(zero, zero, 0, PublicKey);

    vector<vector<MKLweSample*>> enc_image(INPUT_HEIGHT, vector<MKLweSample*>(INPUT_WIDTH));
    
    cout << "\nEncrypting input images . . . " << endl;

    struct timespec img_start = {0, 0};
	struct timespec img_end = {0, 0};

    clock_gettime(CLOCK_MONOTONIC, &img_start);
    for (int i=0; i < INPUT_HEIGHT; i++)
    {
        for (int j=0; j < INPUT_WIDTH; j++)
        {
            int img_val = test_images[0][i][j];
            auto enc_img_val = new_MKLweSample_array(WIDTH, LWEparams, MKparams);
            IntAsymEncrypt(enc_img_val, img_val, crs_a, PublicKey[0]);
            genExtCipherArr(enc_img_val, enc_img_val, 0, PublicKey);
            enc_image[i][j] = enc_img_val;
        }
    }
    clock_gettime(CLOCK_MONOTONIC, &img_end);
    cout << "Input images encryption done." << ": Time --> " << (((double)img_end.tv_nsec + 1.0e+9 * img_end.tv_sec) - ((double)img_start.tv_nsec + 1.0e+9 * img_start.tv_sec)) * 1.0e-9 << " sec" << endl;


    cout << "\nEncrypting CNN model . . . " << endl;

    struct timespec model_start = {0, 0};
	struct timespec model_end = {0, 0};

    vector<vector<vector<MKLweSample*>>> enc_cnn_kernel(NUM_FILTERS, vector<vector<MKLweSample*>>(FILTER_SIZE, vector<MKLweSample*>(FILTER_SIZE)));

    clock_gettime(CLOCK_MONOTONIC, &model_start);
    for (int d = 0; d < NUM_FILTERS; d++)
    {
        for (int i=0; i < FILTER_SIZE; i++)
        {   
            for (int j=0; j < FILTER_SIZE; j++)
            {   
                int cnn_weight_val = cnn_kernel[d][i][j];
                auto enc_cnn_weight = new_MKLweSample_array(WIDTH, LWEparams, MKparams);
                IntAsymEncrypt(enc_cnn_weight, cnn_weight_val, crs_a, PublicKey[1]);
                genExtCipherArr(enc_cnn_weight, enc_cnn_weight, 1, PublicKey);
                enc_cnn_kernel[d][i][j] = enc_cnn_weight;
            }
        }
    }

    vector<MKLweSample*> enc_cnn_bias(NUM_FILTERS);

    for (int d = 0; d < NUM_FILTERS; d++)
    {
        int cnn_bias_val = cnn_bias[d];
        auto enc_cnn_bias_val = new_MKLweSample_array(WIDTH, LWEparams, MKparams);
        IntAsymEncrypt(enc_cnn_bias_val, cnn_bias_val, crs_a, PublicKey[1]);
        genExtCipherArr(enc_cnn_bias_val, enc_cnn_bias_val, 1, PublicKey);
        enc_cnn_bias[d] = enc_cnn_bias_val;
    }

    vector<vector<MKLweSample*>> enc_fc_weight((int)fc_weight.size(), vector<MKLweSample*>((int)fc_weight[0].size()));

    for (int i=0; i < (int)fc_weight.size(); i++)
    {   
        for (int j=0; j < (int)fc_weight[0].size(); j++)
        {   
            int fc_weight_val = fc_weightT[i][j];
            auto enc_fc_weight_val = new_MKLweSample_array(WIDTH, LWEparams, MKparams);
            IntAsymEncrypt(enc_fc_weight_val, fc_weight_val, crs_a, PublicKey[1]);
            genExtCipherArr(enc_fc_weight_val, enc_fc_weight_val, 1, PublicKey);
            enc_fc_weight[i][j] = enc_fc_weight_val;
        }
    }

    vector<MKLweSample*> enc_fc_bias(NUM_CLASSES);

    for (int d = 0; d < NUM_CLASSES; d++)
    {
        int fc_bias_val = fc_bias[d];
        auto enc_fc_bias_val = new_MKLweSample_array(WIDTH, LWEparams, MKparams);
        IntAsymEncrypt(enc_fc_bias_val, fc_bias_val, crs_a, PublicKey[1]);
        genExtCipherArr(enc_fc_bias_val, enc_fc_bias_val, 1, PublicKey);
        enc_fc_bias[d] = enc_fc_bias_val;
    }
    clock_gettime(CLOCK_MONOTONIC, &model_end);
    cout << "CNN model encryption done." << ": Time --> " << (((double)model_end.tv_nsec + 1.0e+9 * model_end.tv_sec) - ((double)model_start.tv_nsec + 1.0e+9 * model_start.tv_sec)) * 1.0e-9 << " sec" << endl;


    time_t start_time = time(NULL);
    printf("Clock starts at: %s", ctime(&start_time));

    cout << "\nPerforming Inference on Encrypted CNN . . ." << endl;

    auto outputs = enc_conv2d(enc_image, enc_cnn_kernel, enc_cnn_bias, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, crs_a, PublicKey, 1);
    cout << "\nConvolution Layer Done." << endl;

    outputs = enc_relu(outputs, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, crs_a, PublicKey, 1);
    cout << "\nRelu Activation Done." << endl;
    
    outputs = enc_max_pooling2d(outputs, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, crs_a, PublicKey, 1);
    cout << "\nMax-Pooling Done." << endl;

    auto fc_input = enc_flatten(outputs);

    auto fc_output = enc_fc_layer(fc_input, enc_fc_weight, enc_fc_bias, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey, crs_a, PublicKey, 1);
    cout << "\nFully-Connected Layer Done . . ." << endl;

    time_t end_time = time(NULL);
    printf("Clock ends at: %s", ctime(&end_time));

    for (int i=0; i < (int)fc_output.size(); i++)
    {
        cout << IntSymDecrypt(fc_output[i], MKlwekey) << " ";
    }

    return 0;

}