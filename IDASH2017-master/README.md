# IDASH2017

IDASH2017 is a project for implementing our Logistic Regression Traning on  encrypted datasets (Privacy-Preserving Logistic Regression Training with A Faster Gradient Variant )

## Pythond codes 
The Python source codes are 'xperiment 11. AdagradWith.25XTXasG vs. Adagrad MNIST.py' and 'Experiment 11. NesterovWith.25XTXasG vs. Nesterov MNIST.py', in the '/data' folder.

## How to run this program? 

### Dependencies

On a Ubuntu cloud, our implementation requires the following libraries in order:
* `g++`:      
```sh
apt install g++ 
```

* `make`:       
```sh
apt install make
```

* `m4`: #        
```sh
apt install m4
```

* `GMP`(ver. 6.1.2):      
```sh
cd gmp-x.x.x  
./configure --enable-cxx  
make
make install
ldconfig
```

* `NTL`(ver. 11.3.0): 
```sh
cd ntl-x.x.x
cd src
./configure NTL_THREADS=on NTL_THREAD_BOOST=on NTL_EXCEPTIONS=on
make
make install
```

### Running IDASH2017

You need to configure and build the CNNinference project. 

After that, in the 'Debug' folder, you can run our project by the following command lines:

```sh
make clean
make all
./MyIDASH2017
``` 

You can change the source codes and then repeat the above lines to debug your own project.

## Running a test source code

In the 'Debug' folder, you can find the C++ running results for six datasets:   

        'MyIDASH2017ArchiveFile20240128_SetNumThreads(36)_kdeg5_numIter4_fold10_Blogp40_idash18x1579.txt_nohup.out'  
        
        'MyIDASH2017ArchiveFile20240128_SetNumThreads(36)_kdeg5_numIter4_fold05_Blogp40_edin.txt_nohup.out'  
        'MyIDASH2017ArchiveFile20240128_SetNumThreads(36)_kdeg5_numIter4_fold05_Blogp40_lbw.txt_nohup.out'  
        'MyIDASH2017ArchiveFile20240128_SetNumThreads(36)_kdeg5_numIter4_fold05_Blogp40_nhanes3.txt_nohup.out'  
        'MyIDASH2017ArchiveFile20240128_SetNumThreads(36)_kdeg5_numIter4_fold05_Blogp40_pcs.txt_nohup.out'  
        'MyIDASH2017ArchiveFile20240128_SetNumThreads(36)_kdeg5_numIter4_fold05_Blogp40_uis.txt_nohup.out'  