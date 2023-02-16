# MKFHE
Multi-Key Homomophic Encryption from TFHE

### Installation of MKFHE

After cloning the repository, do the following: 
```
cd MK-TFHE 
git submodule init
git submodule update 
mkdir build
cd build
cmake ../src -DENABLE_TESTS=on -DENABLE_FFTW=off -DENABLE_NAYUKI_PORTABLE=off -DENABLE_NAYUKI_AVX=off -DCMAKE_BUILD_TYPE=release -DCMAKE_CXX_FLAGS=-fopenmp
make
```

To test (from build):

Threshold decryption on single bit messages,
```
./test/PublicMKTFHE-spqlios-fma
```
Threshold decryption on packed ciphertext,
```
./test/PublicMKTFHEPack-spqlios-fma
```
Threshold decryption on multiple ciphertexts without packing,
```
./test/PublicMKTFHELBit-spqlios-fma
```

CNN inference in clear,
```
./test/cnnInfer_test-spqlios-fma
```

Encrypted CNN inference,
```
./test/encCnnInfer_test-spqlios-fma
```
