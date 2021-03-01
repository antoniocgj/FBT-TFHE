# 1. Installing TFHE

If you already have TFHE installed, you may skip this step.

## 1.1. Clone:

    git clone https://github.com/tfhe/tfhe.git
    cd tfhe

## 1.2. Building:

    mkdir build
    cd build
    cmake ../src -DCMAKE_BUILD_TYPE=optim
    make

# 2. Compiling our code

The following command compiles our code assuming that the TFHE repository is in the parent folder (`../`). If that is not the case, you'll need to adjust the paths `../tfhe/src/include/` and `../tfhe/build/libtfhe/` to point to your install of TFHE. 

    g++ main.cpp -O3 -std=c++11 -funroll-all-loops -march=native -I"../tfhe/src/include/" -L"../tfhe/build/libtfhe/" -ltfhe-spqlios-fma -lm -g -o main

By default, our code will use the `5_5_6_2` parameter set and run each function 10 times. You can change that by defining the constants `P6_4_6_3` and `NUM_EXE`. Examples: 

Example 1 (Run 5 times each function instead of 10):

    g++ main.cpp -O3 -std=c++11 -funroll-all-loops -march=native -I"../tfhe/src/include/" -L"../tfhe/build/libtfhe/" -ltfhe-spqlios-fma -lm -g -o main -DNUM_EXE=5

Example 2 (Use the `6_4_6_3` parameter set instead of the `5_5_6_2`):

    g++ main.cpp -O3 -std=c++11 -funroll-all-loops -march=native -I"../tfhe/src/include/" -L"../tfhe/build/libtfhe/" -ltfhe-spqlios-fma -lm -g -o main -DP6_4_6_3

Exemple 3 (Use the `6_4_6_3` parameter set and run 100 times each function):

    g++ main.cpp -O3 -std=c++11 -funroll-all-loops -march=native -I"../tfhe/src/include/" -L"../tfhe/build/libtfhe/" -ltfhe-spqlios-fma -lm -g -o main -DP6_4_6_3 -DNUM_EXE=100


# 3. Running our code

First, you need to include the `libtfhe` folder in the `LD_LIBRARY_PATH` environment variable. For example (The path to the `libtfhe` folder may need to be adjusted depending of your environment, and it might also need to be an absolute path):

    export LD_LIBRARY_PATH="../tfhe/build/libtfhe/"

Then, you can simply run the compiled program.
  
    ./main

It will execute each function described in the paper. The execution reported for each one of them is the sum of all executions (to calculate the mean, you need to divide the time by the number of executions, `NUM_EXE`).


