# vector commitment scheme with efficient updates
This is an implementation of the vector commitment scheme that is described in Section 5 of the paper [Edrax: A Cryptocurrency with Stateless Transaction Validation](https://eprint.iacr.org/2018/968) by Alexander Chepurnoy (Ergo Platform and IOHK), Charalampos Papamanthou (University of Maryland) and Yupeng Zhang (UC Berkeley). The impelmented vector commitment scheme allows one to (i) represent a vector of 2^l elements with a succinct digest; (b) compute individual proofs for all vector elements; (c) update the digest directly and in constant time whenever a vector element changes; (d) update individual proofs directly and in logarithmic time whenever a vector element changes.

## lead implementor
Yupeng Zhang (UC Berkeley). Please email comments and feedback to jasonzhangyp@gmail.com or cpap@umd.edu.

## build instructions
The implementation relies on the following libraries:
* ate-pairing for bilinear groups. https://github.com/herumi/ate-pairing.
* xbyak, required by ate-pairing https://github.com/herumi/xbyak.
* gmp for big numbers. ```$ sudo apt-get install libgmp3-dev```
* Cmake

The default cmake file assumes ate-pairing is located in the same directory. 

To compile the code, run ```$ cmake .```, then ```$ make```. 

To run the code, execute ```$ ./test [number]```, where 2^[number] specifies the size of the vector. 

test.cpp runs key generation, stores the keys in files, initializes the vector to 0, performs 100 updates and verifies the updated proofs.

