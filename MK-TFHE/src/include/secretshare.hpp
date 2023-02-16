#pragma once

#include <map>
#include <iostream>
#include "tfhe.h"
#include "tfhe_io.h"
#include "lwe-functions.h"
#include "numeric_functions.h"
#include "tlwe_functions.h"
#include "tfhe_garbage_collector.h"
#include <random>
#include <bits/stdc++.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/blas.hpp>
#include <omp.h>
#include <cblas.h>

#define MSIZE 2
#define NSAMPLES 20
namespace ublas = boost::numeric::ublas;

class MKKeyShare {
public:
    std::map<int, LweKey*> shared_key_repo;	/* Stores <group_id>: <key_share> */
    LweKey* GetShare(int group_id);
    // void PartialDecrypt(TLweSample* ciphertext, TLweParams* params, TorusPolynomial* partial_ciphertext, std::vector<int> parties, int t, int p, double sd);
};

class SecretSharing {
private:
    std::map<std::pair<int, int>, int> ncr_cacheT;	            /* Stores <<n, r>: C(n, r)> */
    std::map<std::pair<int, int>, LweKey*> shared_key_repo;	/* Stores <<party_id, group_id>: key_share> */

public:
    void distributeShares(ublas::matrix<int>& S, int t, int k, int p, LweParams *params);
    void shareSecret(int t, int p, LweKey *key, LweParams *params);
    TFheGateBootstrappingSecretKeySet *sk;
    TFheGateBootstrappingCloudKeySet *bk;    
    void GetShareSet(int party, MKKeyShare *share);    
};

void findParties(std::vector<int>& pt, int gid, int t, int p);
int findGroupId(std::vector<int> parties, int t, int p);
// int finalDecrypt(TLweSample* ciphertext, TorusPolynomial** partial_ciphertexts, TLweParams* params, std::vector<int> parties, int t, int p);
void TLweFromLwe(TLweSample *ring_cipher, LweSample *cipher, TLweParams *tlwe_params);
TFheGateBootstrappingParameterSet *initialize_gate_bootstrapping_params();