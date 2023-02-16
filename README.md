# MKFHE
Multi-Key Homomophic Encryption from TFHE

### Installation of MKFHE

After cloning the repository, do the following: 
```
cd MK-TFHE 
git submodule init
git submodule update 
cd OpenBLAS
make
sudo make PREFIX=/usr/local install
cd ..
mkdir build
cd build
cmake ../src -DENABLE_TESTS=on -DENABLE_FFTW=off -DENABLE_NAYUKI_PORTABLE=off -DENABLE_NAYUKI_AVX=off -DCMAKE_BUILD_TYPE=release -DCMAKE_CXX_FLAGS=-fopenmp
make
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
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
For our CNN model training run trainCNN.py file in src/test. We use PyTorch libarary of Python to train our CNN model. Install PyTorch from: https://pytorch.org/
```
pip3 install torch torchvision torchaudio
```
trainCNN.py contains the training of our CNN model and saving the trained parameters in src/test/model_params directory. The MNIST dataset is stored in src/test/mnist_data. Move both the data in src/test/mnist_data and model parameters in src/test/model_params to build/test/mnist_data and build/test/model_params directory respectively.

CNN inference in clear,
```
./test/cnnInfer_test-spqlios-fma
```

Encrypted CNN inference,
```
./test/encCnnInfer_test-spqlios-fma
```
