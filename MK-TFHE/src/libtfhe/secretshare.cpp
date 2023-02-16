#include "secretshare.hpp"
#include <iostream>

std::map<std::pair<int, int>, int> ncr_cacheT;

int ncrT(int n, int r){
    if (ncr_cacheT.find({n, r}) == ncr_cacheT.end()){
        if (r > n || n < 0 || r < 0)
            return 0;
        else{
            if (r == 0 || r == n){
                ncr_cacheT[{n, r}] = 1;
            }else if (r == 1 || r == n - 1){
                ncr_cacheT[{n, r}] = n;
            }else{
                ncr_cacheT[{n, r}] = ncrT(n - 1, r) + ncrT(n - 1, r - 1);
            }
        }
    }
    return ncr_cacheT[{n, r}];
}

// Copies the entirety of src matrix inside dst
// starting from the index (dstR, dstC) in dst
void matrixCopy(ublas::matrix<int> &dst, ublas::matrix<int> &src, int dstR, int dstC){
	ublas::matrix_range<ublas::matrix<int>> dstSlice (dst, ublas::range(dstR, dstR + src.size1()), ublas::range(dstC, dstC + src.size2()));
	dstSlice = src;
}

ublas::matrix<int> optAndCombineT(int t, int k){
	ublas::matrix<int> I;
	I = ublas::identity_matrix<int>(k);

	ublas::matrix<int> Mf;
	int kt = k*t;
	Mf = ublas::zero_matrix<int>(kt, kt);

	// Divide the matrix into t*t chunks of k*k sub-matrices
	for (int r = 0; r < t; r++){
		for (int c = 0; c < t; c++){
			if (r == 0 || c == t - r){
				matrixCopy(Mf, I, r*k, c*k);
			}
		}
	}

	return Mf;
}

ublas::matrix<int> optOrCombineT(int k, int t, int l, ublas::matrix<int> &A){
	ublas::matrix<int> F, R;
	F = ublas::zero_matrix<int>(A.size1(), k);
	R = ublas::zero_matrix<int>(A.size1(), A.size2() - k);

	for (uint r = 0; r < A.size1(); r++){
		for (uint c = 0; c < A.size2(); c++){
			if (c < (uint)k){
				F(r, c) = A(r, c);
			}else{
				R(r, c - k) = A(r, c);
			}
		}
	}

	ublas::matrix<int> M;
	M = ublas::zero_matrix<int>(l*k*t, k*(t-1) * l + k);
	for (int i = 0; i < l; i++){
		matrixCopy(M, F, i*k*t, 0);
		matrixCopy(M, R, i*k*t, k + i*k*(t-1));
	}

	return M;
}

/* Build Distribution matrix for OR-ing C(p,t) number of x1x2...xt like terms. */
void buildDistributionMatrix(int t, int k, int p, ublas::matrix<int>& M){
	ublas::matrix<int> M1;
	M1 = optAndCombineT(t, k);
	M = optOrCombineT(k, t, ncrT(p, t), M1);
}

/* rho is a random binary matrix with first k rows coming from the k rows of the secret key */
ublas::matrix<int> buildRho(int k, int p, int e, LweKey *key, ublas::matrix<int>& rho){
	int n = key->params->n;
	rho = ublas::matrix<int>(e, n);
	for(int row = 0; row < k; row++){
		for(int col = 0; col < n; col++){
			rho(row,col) = key->key[col];
		}
	}
	std::default_random_engine gen;
    std::uniform_int_distribution<int> dist(0, 1);
	for(int row = k; row < e; row++){
		for(int col = 0; col < n; col++){
			rho(row,col) = dist(gen);
		}
	} 
	return rho;
}

/* Naive Matrix Multiplication */
void multiply(ublas::matrix<int>& C, ublas::matrix<int>& A, ublas::matrix<int>& B){
	int m = A.size1();
	int k = A.size2();
	int n = B.size2();

	double *_A = (double *)malloc(m * k * sizeof(double));
	double *_B = (double *)malloc(k * n * sizeof(double));
	double *_C = (double *)malloc(m * n * sizeof(double));

	memset(_A, 0, m * k * sizeof(double));
	memset(_B, 0, k * n * sizeof(double));
	memset(_C, 0, m * n * sizeof(double));

	for (int r = 0; r < m; r++)
		for (int c = 0; c < k; c++)
			_A[r * k + c] = A(r, c);

	for (int r = 0; r < k; r++)
		for (int c = 0; c < n; c++)
			_B[r * n + c] = B(r, c);


	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, _A, k, _B, n, 0.0, _C, n);

	for (int r = 0; r < m; r++)
		for (int c = 0; c < n; c++)
			C(r, c) = _C[r * n + c];

	free(_A);
	free(_B);
	free(_C);
}

/* Given a group_id, find the party_ids present in (group_id)^th combination out of C(p,t) combinations */
void findParties(std::vector<int>& pt, int gid, int t, int p)
{
	int mem = 0, tmp;
	pt.clear();
	for(int i = 1; i < p; i++){
		tmp = ncrT(p - i, t - mem -1);
		if(gid > tmp){
			gid -= tmp;
		}
		else{
			pt.push_back(i);
			mem += 1;
		}
		if(mem + (p-i) == t){
			for(int j = i + 1; j <= p; j++){
				pt.push_back(j);
			}
			break;
		}
	}
}

/* Given a t-sized list of party-ids compute its rank among total C(p,t) combinations */
int findGroupId(std::vector<int> parties, int t, int p)
{
	int mem = 0;
	int group_count = 1;
	for(int i = 1; i <= p; i++){
		if(std::find(parties.begin(), parties.end(), i) != parties.end()){
			mem += 1;
		}
		else{
			group_count += ncrT(p - i, t - mem - 1);
		}
		if(mem == t){
			break;
		}
	}
	return group_count;
}



/* Get and store the actual shares of each party */
void SecretSharing::distributeShares(ublas::matrix<int>& S, int t, int k, int p, LweParams *params)
{
	int r = S.size1(), n = params->n;
	int row = 1, group_id, row_count;
	std::vector<int> parties;
	while(row <= r){
		group_id = ceil(row/(floor)(k*t));
		findParties(parties, group_id, t, p);
		for(int it = 1; it <= t; it++){
			row_count = row + (it - 1) * k;
			LweKey *key_share = new_LweKey(params);
			for(int i = 0; i < k; i++){
				for(int j = 0; j < n; j++){
					key_share->key[j] = S(row_count + i - 1, j);
				}
			}
			shared_key_repo[{parties[it-1], group_id}] = key_share;
		}
		row += (k*t);
	}
}

/* Preprocess to share the secret key among p parties */
void SecretSharing::shareSecret(int t, int p, LweKey *key, LweParams *params)
{
	int k = 1;
	int n = key->params->n;

	ublas::matrix<int> M;
	buildDistributionMatrix(t, k, p, M);
	int d = M.size1();
	int e = M.size2();

	ublas::matrix<int> rho;
	rho = buildRho(k, p, e, key, rho);

	ublas::matrix<int> shares(d, n);
	multiply(shares, M, rho);	/* shares = M . rho */

	this->distributeShares(shares, t, k, p, params);
}


TFheGateBootstrappingParameterSet *initialize_gate_bootstrapping_params()
{
    static const int32_t N = 1024;
    static const int32_t k = 1;
    static const int32_t n = 560;
    static const int32_t bk_l = 3;
    static const int32_t bk_Bgbit = 7;
    static const int32_t ks_basebit = 2;
    static const int32_t ks_length = 8;
    static const double ks_stdev = pow(2.,-15); //standard deviation
    static const double bk_stdev = pow(2.,-25);; //standard deviation
    static const double max_stdev = 0.012467; //max standard deviation for a 1/4 msg space

    LweParams *params_in = new_LweParams(n, ks_stdev, max_stdev);
    TLweParams *params_accum = new_TLweParams(N, k, bk_stdev, max_stdev);
    TGswParams *params_bk = new_TGswParams(bk_l, bk_Bgbit, params_accum);

    TfheGarbageCollector::register_param(params_in);
    TfheGarbageCollector::register_param(params_accum);
    TfheGarbageCollector::register_param(params_bk);

    return new TFheGateBootstrappingParameterSet(ks_length, ks_basebit, params_in, params_bk);
}

void TLweFromLwe(TLweSample *ring_cipher, LweSample *cipher, TLweParams *tlwe_params)
{
	int N = tlwe_params->N;
	ring_cipher->a[0].coefsT[0] = cipher->a[0];
	ring_cipher->b->coefsT[0] = cipher->b;
	for(int i = 1; i < N; i++){
		ring_cipher->a[0].coefsT[i] = -cipher->a[N-i];
	}
}

void TLweKeyFromLweKey(const LweKey *lwe_key, TLweKey *tlwe_key){
	int N = tlwe_key->params->N;
	tlwe_key->key[0].coefs[0] = lwe_key->key[0];
	for(int i = 0; i < N; i++){
		tlwe_key->key[0].coefs[i] = lwe_key->key[i];
	}
}


void SecretSharing::GetShareSet(int party, MKKeyShare *share)
{
    for (auto &pr : this->shared_key_repo){
        if (pr.first.first == party){
            share->shared_key_repo[pr.first.second] = pr.second;
        }
    }
}

LweKey* MKKeyShare::GetShare(int group_id)
{
    return this->shared_key_repo[group_id];
}
