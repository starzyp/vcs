# Vector Commitment Scheme
## Authors
Alexander Chepurnoy, Charalampos Papamanthou, Yupeng Zhang
## Build Instructions
VCS relies on the following libraries
* ate-pairing for bilinear groups. https://github.com/herumi/ate-pairing.
* xbyak, required by ate-pairing https://github.com/herumi/xbyak.
* gmp for big numbers. ```$ sudo apt-get install libgmp3-dev```
* Cmake

The default cmake file assumes ate-pairing is located in the same directory. To compile the code, run ```$ cmake .```, then ```$ make```. 
To run the code, execute ```$ ./test [number]```, where 2^[number] is the size of the vector. test.cpp runs key generation, stores the keys in files, initializes the vector to 0, perform 100 updates (update value, digest, proof) and verifies the proofs.

## References

[Edrax: A Cryptocurrency with Stateless Transaction Validation](https://eprint.iacr.org/2018/968). Alexander Chepurnoy, Charalampos Papamanthou and Yupeng Zhang.
