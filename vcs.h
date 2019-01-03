#include <vector>
#include "test_point.hpp"
#include "bn.h"
#include <gmp.h>
#include <gmpxx.h>
#include <fstream>
#include <thread>

using namespace std;
using namespace bn;


class vcs{
	public:
	vcs(int, mpz_class, Ec1, Ec2);
	~vcs();
	
	
	mpz_class p; //p is the prime that defines the field, and it must be the same as the base group of the bilinear group.
	
	int L,N,P;

	Ec1 g1;
	Ec2 g2;
	
	//L is the number of variables and N=2^L is the number of elements in the vector.
	//P is the number of bits in p. It is used for fast exponentiation during keygen
	

	vector<Ec1> calc_update_key(long long int index, vector<vector<Ec1> >& prk);
	vector<vector<Ec1> > calc_update_key_batch(vector<long long int> index, vector<vector<Ec1> >& prk);

	void keygen(vector<vector<Ec1> >& prk, vector<Ec2>& vrk);
	void load_key(vector<vector<Ec1> >& prk, vector<Ec2>& vrk);
	
	Ec1 setup(vector<mpz_class>& a, vector<vector<Ec1> >& prk);

	vector<Ec1> prove(int index, vector<mpz_class>& a, vector<vector<Ec1> >& prk);
	bool verify(Ec1 digest, int index, mpz_class a_i, vector<Ec1> proof, vector<Ec2> vrk);
	bool batch_verify(Ec1 digest, vector<int> index, vector<mpz_class> a_i, vector<vector<Ec1> > proof, vector<Ec2> vrk);
	Ec1 update_digest(Ec1 digest, int updateindex, mpz_class delta, vector<Ec1> upk_u);
	vector<Ec1> update_proof(vector<Ec1> proof, int updateindex, int index, mpz_class delta, vector<Ec1> upk_u);
};

