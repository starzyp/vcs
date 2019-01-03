#include "vcs.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <time.h>
#include <map>

#include "test_point.hpp"
#include "bn.h"

#include <gmp.h>
#include <gmpxx.h>

using namespace std;
using namespace bn;

int main(int argc, char** argv){
	
	int L = atoi(argv[1]);

	bn::CurveParam cp = bn::CurveFp254BNb;
	Param::init(cp);
	const Point& pt = selectPoint(cp);
	const Ec2 g2(
		Fp2(Fp(pt.g2.aa), Fp(pt.g2.ab)),
		Fp2(Fp(pt.g2.ba), Fp(pt.g2.bb))
	);
	const Ec1 g1(pt.g1.a, pt.g1.b);
	
	mpz_class p;
	p.set_str("16798108731015832284940804142231733909759579603404752749028378864165570215949",10);
	
	srand(time(NULL));
	gmp_randstate_t r_state;
	unsigned long int seed = rand();
	gmp_randinit_default (r_state);
    gmp_randseed_ui(r_state, seed);
	
	
	vcs a(L,p,g1,g2);
	
	auto t2 = chrono::high_resolution_clock::now();
	
	vector<vector<Ec1> > prk;
	vector<Ec2> vrk;

	a.keygen(prk, vrk);
	a.load_key(prk,vrk);
	
	auto t3 = chrono::high_resolution_clock::now();
	auto t4 = t3 - t2;
	cout << "keygen time: " << chrono::duration<double, milli>(t4).count()/1000 << "s" << endl;
	
	
	map<int, mpz_class> b;
	
	Ec1 digest = g1*0;
	
	
	int num_txs = 1024;
	
	vector<int> index(num_txs);
	
	
	vector<vector<Ec1> > proof(num_txs);
	for(int i=0;i<num_txs;i++){
		index[i] = rand()%a.N;
		b[index[i]] = 0;
		proof[i]=vector<Ec1>(L,g1*0);
	}
	
	
	t2 = chrono::high_resolution_clock::now();
	
	vector<long long int> updateindex(100);
	for(int round =0;round<updateindex.size();round++)
		updateindex[round] = rand()%a.N;
	
	vector<vector<Ec1> > upk_ui = a.calc_update_key_batch(updateindex, prk);
	
	t3 = chrono::high_resolution_clock::now();
	t4 = t3 - t2;
	cout << "load update key time: " << chrono::duration<double, milli>(t4).count()/1000 << "s" << endl;
	
	double update_digest_time=0.0, update_proof_time = 0.0;
	
	for(int round=0;round<updateindex.size();round++){
		
		mpz_class delta;
		delta = rand()%10000;
		if(b.find(updateindex[round]) == b.end())
			b[updateindex[round]]=delta;
		else
			b[updateindex[round]]+=delta;
		
		
		t2 = chrono::high_resolution_clock::now();
		
		digest = a.update_digest(digest, updateindex[round], delta, upk_ui[round]);
		
		t3 = chrono::high_resolution_clock::now();
		t4 = t3 - t2;
		
		update_digest_time+=chrono::duration<double, milli>(t4).count()/1000;
		
		t2 = chrono::high_resolution_clock::now();
		
		for(int j=0;j<num_txs;j++){
			proof[j] = a.update_proof(proof[j], updateindex[round], index[j], delta, upk_ui[round]);
		}
		
		
		t3 = chrono::high_resolution_clock::now();
		t4 = t3 - t2;
		
		update_proof_time+=chrono::duration<double, milli>(t4).count()/1000/num_txs;
		

	}
	
	cout<<"average update digest time: "<<update_digest_time<<"s\n";
	cout<<"average update proof time: "<<update_proof_time<<"s\n";
	
	
	clock_t t1 = clock();
	bool tf3 = 1;
	for(int i=0;i<num_txs;i++)
		tf3 = tf3 & a.verify(digest, index[i], b[index[i]], proof[i], vrk);
	cout<<"verification: "<<tf3<<endl;
	cout<<"verify time: "<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s\n";
	
	
	vector<mpz_class> a_i(num_txs);
	for(int i=0; i<num_txs;i++){
		a_i[i]=b[index[i]];
	}
	t1 = clock();
	t2 = chrono::high_resolution_clock::now();
	
	
	
	
	bool tf4 = a.batch_verify(digest, index, a_i, proof, vrk);
	
	cout<<"batch verification: "<<tf4<<endl;
	
	t3 = chrono::high_resolution_clock::now();
	t4 = t3 - t2;
	cout << "verify time 2: " << chrono::duration<double, milli>(t4).count()/1000 << "s" << endl;
	
	
	
	return 1;
}