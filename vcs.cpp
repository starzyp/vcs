#include "vcs.h"

#include <cstring>
#include <string>
#include <cmath>
#include <sys/stat.h>


#define ncore 16
#define nfiles 8
#define lognfiles 3

string path = "pkvk/";

vector<bool> to_binary(int index, int L){ //LSB first
	vector<bool> binary(L);
	for(int i=0;i<L;i++){
		binary[i] = index%2;
		index/=2;
	}
	return binary;
}


void precompute_g1(Ec1 g1, vector<Ec1>& g1_pre, int P){
	g1_pre.resize(P);
	g1_pre[0] = g1;
	for(int i=1;i<P;i++){
		g1_pre[i]=g1_pre[i-1]+g1_pre[i-1];
	}
	return;
}


template< class T >
T pre_exp(vector<T>& pre, mpz_class n){
	T temp = pre[0]*0;
	int length = mpz_sizeinbase(n.get_mpz_t(), 2);

	for(int i=0;i<length;i++){
		if(mpz_tstbit(n.get_mpz_t(),i)==1)
			temp = temp + pre[i];		
	}
	
	
	return temp;

}


vcs::vcs(int d, mpz_class p, Ec1 g1, Ec2 g2){
	
	L = d;
	N = (int)pow(2,L);
	
	
	this->p = p;
	//p.set_str("16798108731015832284940804142231733909759579603404752749028378864165570215949",10);
	P=mpz_sizeinbase(p.get_mpz_t(),2);
	
	this->g1 = g1;
	this->g2 = g2;
}

vcs::~vcs(){}

void vcs::keygen(vector<vector<Ec1> >& prk, vector<Ec2>& vrk){
	
	unsigned long int seed;
	gmp_randstate_t r_state;
	short size = sizeof(seed);
	ifstream urandom("/dev/urandom", ios::in|ios::binary);
	urandom.read((char*)&seed,size);
	urandom.close();
	
	
    gmp_randinit_default (r_state);
    gmp_randseed_ui(r_state, seed);
	
	
	
	vector<Ec1> g1_pre;
	precompute_g1(g1, g1_pre,P);
	
	//secret keys
	vector<mpz_class> s(L);
	
	for(int i=0;i<L;i++)
		mpz_urandomm(s[i].get_mpz_t(),r_state,p.get_mpz_t());
	
	
	
	//compute public keys
	
	prk.resize(L+1);
	vrk.resize(L);

	
	for(int i=0;i<L+1;i++){
		if(i>lognfiles)
			prk[i].resize((int)pow(2,i-lognfiles));
		else
			prk[i].resize((int)pow(2,i));
	}
	
	prk[0][0] = g1;
	
	
	
	vector<vector<mpz_class> > vars(L+1);
	for(int i=0;i<L+1;i++){
		if(i>lognfiles)
			vars[i].resize((int)pow(2,i-lognfiles));
		else
			vars[i].resize((int)pow(2,i));
	}
	

	
	vars[0][0]=1;

	
	
	for(int i=1;i<lognfiles+1;i++){
		for(int j=0;j<(int)pow(2,i-1);j++){
			vars[i][2*j+1] = (vars[i-1][j]*s[i-1])%p;			
			vars[i][2*j] = (vars[i-1][j]*(1-s[i-1]))%p;
		}
	}
	
	
	for(int i=1;i<lognfiles+1;i++){
		for(int j=0;j<(int)pow(2,i-1);j++){
				
			if(vars[i][2*j+1]<0)
				vars[i][2*j+1]+=p;

				
			prk[i][2*j+1] = pre_exp(g1_pre,vars[i][2*j+1]);
			prk[i][2*j] = prk[i-1][j]-prk[i][2*j+1];
			
			
		}
		
		
	}
	
	auto f = [](int i, int x, int y, mpz_class p, vector<Ec1>* g1_pre, vector<vector<mpz_class> >* vars, vector<vector<Ec1> >* prk) {
        for (int j = x; j < y; j++){
			if((*vars)[i][2*j+1]<0)
				(*vars)[i][2*j+1]+=p;

				
			(*prk)[i][2*j+1] = pre_exp(*g1_pre,(*vars)[i][2*j+1]);
			(*prk)[i][2*j] = (*prk)[i-1][j]-(*prk)[i][2*j+1];
		
		}
    };
	

	for(int batch = 0; batch < nfiles; batch++){
		
		
		vars[lognfiles+1][1] = (vars[lognfiles][batch]*s[lognfiles])%p;			
		vars[lognfiles+1][0] = (vars[lognfiles][batch]*(1-s[lognfiles]))%p;

		
		for(int i=lognfiles+2;i<L+1;i++){
			
			for(int j=0;j<(int)pow(2,i-1-lognfiles);j++){
				vars[i][2*j+1] = (vars[i-1][j]*s[i-1])%p;			
				vars[i][2*j] = (vars[i-1][j]*(1-s[i-1]))%p;
			}
		}
		
		if(vars[lognfiles+1][1]<0)
			vars[lognfiles+1][1]+=p;
					
		prk[lognfiles+1][1] = pre_exp(g1_pre,vars[lognfiles+1][1]);
		prk[lognfiles+1][0] = prk[lognfiles][batch]-prk[lognfiles+1][1];
				
	
		
		for(int i=lognfiles+2;i<lognfiles+(int)log2(ncore)+1;i++){
			for(int j=0;j<(int)pow(2,i-1-lognfiles);j++){
					
				if(vars[i][2*j+1]<0)
					vars[i][2*j+1]+=p;

					
				prk[i][2*j+1] = pre_exp(g1_pre,vars[i][2*j+1]);
				prk[i][2*j] = prk[i-1][j]-prk[i][2*j+1];
				
				
			}
		
		
		}
		
		thread th[ncore];
	
		for(int i=lognfiles+(int)log2(ncore)+1;i<L+1;i++){

			int total_size = (int)pow(2,i-1-lognfiles);
			for(int k=0;k<ncore;k++)
				th[k]=thread(f,i,total_size/ncore*k, total_size/ncore*(k+1),p ,&g1_pre, &vars, &prk);
				
			for(int k=0;k<ncore;k++)
				th[k].join();	
				
		
				
		}
		
		mkdir(path.c_str(),S_IRWXU);
		
		ofstream OutFile;
		string filename = path+"pk"+to_string(batch)+".txt";
		
		
		
		OutFile.open(filename, ios::out | ios::binary);
		
		
		
		for(int i= lognfiles+1; i<L+1;i++){
			for(int j=0;j<prk[i].size();j++){
				OutFile.write( (char*)&prk[i][j], sizeof(Ec1));
			
			}
		}
		
		OutFile.close();
	
	}
	
	for(int i=lognfiles+1;i<L+1;i++){
		prk[i].resize(0);
	}
	
	ofstream OutFile;
	string filename = path+"pk.txt";
	OutFile.open(filename, ios::out | ios::binary);
	
	
	
	for(int i=0 ; i<lognfiles+1;i++){
		for(int j=0;j<prk[i].size();j++){
			OutFile.write( (char*)&prk[i][j], sizeof(Ec1));
		
		}
	}
	
	OutFile.close();
	
	
	
	for(int i=0;i<L;i++){
		const mie::Vuint temp((s[i].get_str()).c_str());
		vrk[i]=g2*temp;
	}
	
	filename = path+"vrk.txt";
	OutFile.open(filename, ios::out | ios::binary);
	
	
	
	for(int i=0 ; i<vrk.size();i++){
		OutFile.write( (char*)&vrk[i], sizeof(Ec2));
	}
	
	OutFile.close();
	
	
	/*
	upk.resize(N);
	
	for(int i=0;i<N;i++){
		upk[i].resize(L);
		for(int j=L-1;j>=0;j--){
			upk[i][j] = prk[j+1][i>>(L-j-1)];
		}
	} */
	
	return;
}

void vcs::load_key(vector<vector<Ec1> >& prk, vector<Ec2>& vrk){
	prk.resize(lognfiles+1);
	vrk.resize(L);
	
	ifstream InFile;
	string filename = path+"pk.txt";
	InFile.open(filename, ios::in | ios::binary);
	
	
	for(int k=0 ; k<lognfiles+1;k++){
		prk[k].resize((int)pow(2,k));
		for(int t=0;t<prk[k].size();t++){
			InFile.read( (char*)&prk[k][t], sizeof(Ec1));
		
		}
	}
	
	InFile.close();
	
	filename = path+"vrk.txt";
	InFile.open(filename, ios::in | ios::binary);
	
	
	for(int k=0 ; k<L;k++){
		InFile.read( (char*)&vrk[k], sizeof(Ec2));
	}
	
	InFile.close();
	return;

}

vector<Ec1> vcs::calc_update_key(long long int index, vector<vector<Ec1> >& prk){
    vector<Ec1> upk;
    upk.resize(L);
	
	
	for(int j=lognfiles-1;j>=0;j--){
    	upk[j] = prk[j+1][index >> (L-j-1)];
    }
	
	
    for(int j=L-1;j>=lognfiles;j--){
		
		
		
		int filenum = (index >> (L-j-1))*nfiles/(int)pow(2,j+1);
		int keynum = (index >> (L-j-1))% ((int)pow(2,j+1)/nfiles);
		
		//cout<<index<<" "<<filenum<<" "<<keynum<<endl;
		
		
		vector<vector<Ec1> > load_prk(L+1);
		for(int k=lognfiles+1;k<L+1;k++)
			load_prk[k].resize((int)pow(2,k-lognfiles));
		
		
		ifstream InFile;
		string filename = path+"pk"+to_string(filenum)+".txt";
		InFile.open(filename, ios::in | ios::binary);
		
		
		for(int k= lognfiles+1; k<L+1;k++){
			for(int t=0;t<load_prk[k].size();t++){
				InFile.read( (char*)&load_prk[k][t], sizeof(Ec1));
			
			}
		}
		
		InFile.close();
	
    	upk[j] = load_prk[j+1][keynum];
    }
    return upk;
}

vector<vector<Ec1> > vcs::calc_update_key_batch(vector<long long int> index, vector<vector<Ec1> >& prk){
    vector<vector<Ec1> > upk;
    upk.resize(index.size());
	for(int i=0;i<index.size();i++)
		upk[i].resize(L);
	
	for(int i=0;i<upk.size();i++){
		for(int j=lognfiles-1;j>=0;j--){
			upk[i][j] = prk[j+1][index[i] >> (L-j-1)];
		}
	}
	
	
	for(int filenum = 0; filenum<nfiles; filenum++){
	
		
		
		vector<vector<Ec1> > load_prk(L+1);
		for(int k=lognfiles+1;k<L+1;k++)
			load_prk[k].resize((int)pow(2,k-lognfiles));
		
		ifstream InFile;
		string filename = path+"pk"+to_string(filenum)+".txt";
		InFile.open(filename, ios::in | ios::binary);
		
		
		for(int k= lognfiles+1; k<L+1;k++){
			for(int t=0;t<load_prk[k].size();t++){
				InFile.read( (char*)&load_prk[k][t], sizeof(Ec1));
			
			}
		}
		
		
		InFile.close();
		

	
		for(int i=0;i<index.size();i++){
			for(int j=L-1;j>=lognfiles;j--){

			
				int filenum_temp = (index[i] >> (L-j-1))*nfiles/(int)pow(2,j+1);
				int keynum = (index[i] >> (L-j-1))% ((int)pow(2,j+1)/nfiles);
				
				if(filenum_temp == filenum){
					upk[i][j] = load_prk[j+1][keynum];
				
				}
				
			}
		}
	
	}

    return upk;
}


Ec1 vcs::setup(vector<mpz_class>& a, vector<vector<Ec1> >& prk){

	Ec1 digest = g1*0;
	for(int i=0;i<N;i++){
		if(a[i]!=0){
			if(a[i]==1){
				digest = digest+ prk[L][i];
			}
			else{
				if(a[i]>=0){
					const mie::Vuint temp((a[i].get_str()).c_str());
					digest = digest+prk[L][i]*temp;
				}
				else{
					mpz_class temp2= -a[i];
					const mie::Vuint temp((temp2.get_str()).c_str());
					digest = digest-prk[L][i]*temp;
				}
			}
		}
	}
	
	return digest;

}


vector<Ec1> vcs::prove(int index, vector<mpz_class>& a, vector<vector<Ec1> >& prk){

	vector<bool> index_binary = to_binary(index,L);
	
	vector<mpz_class> witness_coeffs(N), temp_coeffs = a;
	
	int start_index = 0;
	
	for(int i=0;i<L;i++){
		for(int j=0;j<pow(2,L-i-1);j++){
			witness_coeffs[start_index+j] = (-temp_coeffs[2*j]+temp_coeffs[2*j+1])%p;
			temp_coeffs[j] = (-temp_coeffs[2*j]*(index_binary[i]-1)+temp_coeffs[2*j+1]*index_binary[i])%p;
		}
		temp_coeffs.resize(pow(2,L-i-1));
		start_index+=pow(2,L-i-1);
	}
	
	
	vector<Ec1> witness(L);
	
	start_index = 0;
	
	for(int i=0;i<L;i++){
		witness[i] = g1*0;
		
		for(int j=0;j<pow(2,L-i-1);j++){
			if(witness_coeffs[start_index+j]>0){
				const mie::Vuint temp((witness_coeffs[start_index+j].get_str()).c_str());
				witness[i] += prk[L-i-1][j]*temp;	
			}
			else{
				witness_coeffs[start_index+j]=-witness_coeffs[start_index+j];
				const mie::Vuint temp((witness_coeffs[start_index+j].get_str()).c_str());
				witness[i] -= prk[L-i-1][j]*temp;	
			}
		}
		
		start_index+=pow(2,L-i-1);
		
	}
	
	return witness;

}

bool vcs::verify(Ec1 digest, int index, mpz_class a_i, vector<Ec1> proof, vector<Ec2> vrk){
	
	
	
	Fp12 e1,e3=1;
	vector<Fp12> e2(L);
	vector<bool> index_binary=to_binary(index,L);
	
	if(a_i>=0){
		const mie::Vuint temp2((a_i.get_str()).c_str());
		Ec1 temp = g1*temp2;
		opt_atePairing(e1,g2, digest-temp);
	}
	else{
		a_i=-a_i;
		const mie::Vuint temp2((a_i.get_str()).c_str());
		Ec1 temp = g1*temp2;
		opt_atePairing(e1,g2, digest+temp);
	}
	
	for(int i=0;i<L;i++){
		Ec2 temp1 = vrk[L-i-1]-g2*(int)index_binary[i];
		
		
		//clock_t t1=clock();
		opt_atePairing(e2[i],temp1, proof[i]);
		//if(i==0){
			//cout<<"pairing time: "<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s\n";
		//}
		
		e3*=e2[i];
		
	}
	
	/*
	
	unsigned long int seed;
	gmp_randstate_t r_state;
	short size = sizeof(seed);
	ifstream urandom("/dev/urandom", ios::in|ios::binary);
	urandom.read((char*)&seed,size);
	urandom.close();
	mpz_class rdomain;
	rdomain.set_str("1208925819614629174706175",10); //2^80-1
	//rdomain.set_str("340282366920938463463374607431768211455",10); //2^128-1
    gmp_randinit_default (r_state);
    gmp_randseed_ui(r_state, seed);
	
	mpz_class test;
	mpz_urandomm(test.get_mpz_t(),r_state,rdomain.get_mpz_t());
	const mie::Vuint temp3((test.get_str()).c_str());
	clock_t t1=clock();
	proof[0] = proof[0]*temp3;
	cout<<"exp time: "<<(double)(clock()-t1)/CLOCKS_PER_SEC<<"s\n";
	*/
	
	return (e1==e3);
	
}

bool vcs::batch_verify(Ec1 digest, vector<int> index, vector<mpz_class> a_i, vector<vector<Ec1> > proof, vector<Ec2> vrk){
	
	clock_t t1=clock();
	
	// to binary
	vector<vector<bool> > index_binary(index.size());
	for(int i=0;i<index.size();i++)
		index_binary[i]=to_binary(index[i],L);
	
	//random seed
	unsigned long int seed;
	gmp_randstate_t r_state;
	short size = sizeof(seed);
	ifstream urandom("/dev/urandom", ios::in|ios::binary);
	urandom.read((char*)&seed,size);
	urandom.close();
	mpz_class rdomain;
	//rdomain.set_str("1208925819614629174706175",10); //2^80-1
	rdomain.set_str("340282366920938463463374607431768211455",10); //2^128-1
    gmp_randinit_default (r_state);
    gmp_randseed_ui(r_state, seed);

	
	//randomness
	
	vector<mpz_class> r(index.size());
	for(int i=0;i<index.size();i++)
		mpz_urandomm(r[i].get_mpz_t(),r_state,rdomain.get_mpz_t());
	

	
	//proof^randomness
	
	auto t2 = chrono::high_resolution_clock::now();
	
	auto f = [](int x, int y, vector<vector<Ec1> >* proof, vector<mpz_class>* r, int L) {
        for(int i=x;i<y;i++){
			const mie::Vuint temp(((*r)[i].get_str()).c_str());
			for(int j=0;j<L;j++){
				(*proof)[i][j] = (*proof)[i][j]*temp;
			}
		}
    };
	
	thread th[ncore];
	
	for(int k=0;k<ncore;k++)
		th[k]=thread(f,index.size()/ncore*k, index.size()/ncore*(k+1),&proof, &r, L);
		
	/*
	
	vector<vector<Ec1> > proof(L);
	for(int i=0;i<L;i++){
		proof[i].resize(index.size());
		for(int j=0;j<index.size();j++)
			proof[i][j]=proof2[j][i];
	
	}
	
	auto f = [](vector<Ec1>* proof, vector<mpz_class> *r) {
		
		for(int j=0;j<(*proof).size();j++){
			const mie::Vuint temp(((*r)[j].get_str()).c_str());
			(*proof)[j] = (*proof)[j]*temp;
		}
    };
	
	thread th[L];
	
	
	
	for(int k=0;k<L;k++)
		th[k]=thread(f,&proof[k], &r);
	*/
		
	for(int k=0;k<ncore;k++)
		th[k].join();	
	
	auto t3 = chrono::high_resolution_clock::now();
	auto t4 = t3 - t2;
	//cout << "verify time: " << chrono::duration<double, milli>(t4).count()/1000 << "s" << endl;
	
	/*
	for(int i=0;i<index.size();i++){
		const mie::Vuint temp((r[i].get_str()).c_str());
		for(int j=0;j<L;j++){
			proof[i][j] = proof[i][j]*temp;
		
		}
	}
	*/
	
	
	//left side
	
	Fp12 e1,e3=1,e2;
	//vector<Fp12> e2(L);
	
	mpz_class a=0, r_sum=0;
	for(int i=0;i<index.size();i++){
		a+=a_i[i]*r[i];
		r_sum+=r[i];
	}
	
	if(a<0)
		a+=p;
		
	if(r_sum<0)
		r_sum+=p;
	
	
	const mie::Vuint temp2((a.get_str()).c_str());
	const mie::Vuint temp3((r_sum.get_str()).c_str());
	opt_atePairing(e1,g2, digest*temp3-g1*temp2);
	
	
	//right side
	vector<Ec1> proof_combined(2*L, g1*0);
	for(int i=0;i<index.size();i++){
		for(int j=0;j<L;j++){
			if(index_binary[i][j]==0){
				proof_combined[2*j]=proof_combined[2*j]+proof[i][j];
			}
			else
				proof_combined[2*j+1]=proof_combined[2*j+1]+proof[i][j];
		}
	}
	
	for(int i=0;i<L;i++){
		Ec2 temp1 = vrk[L-i-1]-g2*0;
		opt_atePairing(e2,temp1, proof_combined[2*i]);
		
		e3*=e2;
		
		temp1 = vrk[L-i-1]-g2*1;
		opt_atePairing(e2,temp1, proof_combined[2*i+1]);
		
		e3*=e2;
		
	}
	
	
	return (e1==e3);
	
}

Ec1 vcs::update_digest(Ec1 digest, int updateindex, mpz_class delta, vector<Ec1> upk_u){
	if(delta>=0){
		const mie::Vuint temp((delta.get_str()).c_str());
		return digest+upk_u[L-1]*temp;
	}
	else{
		delta=-delta;
		const mie::Vuint temp((delta.get_str()).c_str());
		return digest-upk_u[L-1]*temp;
	}

}

vector<Ec1> vcs::update_proof(vector<Ec1> proof, int updateindex, int index, mpz_class delta, vector<Ec1> upk_u){
	vector<Ec1> new_proof=proof;
	vector<bool> index_binary=to_binary(index,L), updateindex_binary=to_binary(updateindex,L);
	
	if(delta>=0){
		const mie::Vuint temp((delta.get_str()).c_str());
		for(int i=0;i<L;i++){
		
			if(i<L-1){
				if(updateindex_binary[i] == 0 && index_binary[i]==1){
					new_proof[i]=proof[i]-upk_u[L-i-2]*temp;
					break;
				}
				else if(updateindex_binary[i] == 1 && index_binary[i]==0){
					new_proof[i]=proof[i]+upk_u[L-i-2]*temp;
					break;
				}
				else if(updateindex_binary[i] == 0 && index_binary[i]==0){
					new_proof[i]=proof[i]-upk_u[L-i-2]*temp;
				}
				else{
					new_proof[i]=proof[i]+upk_u[L-i-2]*temp;
				}
			}
			else{
				if(updateindex_binary[i] == 0 && index_binary[i]==1)
					new_proof[i]=proof[i]-g1*temp;
				else if(updateindex_binary[i] == 1 && index_binary[i]==0)
					new_proof[i]=proof[i]+g1*temp;
				else if(updateindex_binary[i] == 0 && index_binary[i]==0)
					new_proof[i]=proof[i]-g1*temp;
				else
					new_proof[i]=proof[i]+g1*temp;
			}
		
		}
	}
	else{
		delta=-delta;
		const mie::Vuint temp((delta.get_str()).c_str());
		for(int i=0;i<L;i++){
		
			if(i<L-1){
				if(updateindex_binary[i] == 0 && index_binary[i]==1){
					new_proof[i]=proof[i]+upk_u[L-i-2]*temp;
					break;
				}
				else if(updateindex_binary[i] == 1 && index_binary[i]==0){
					new_proof[i]=proof[i]-upk_u[L-i-2]*temp;
					break;
				}
				else if(updateindex_binary[i] == 0 && index_binary[i]==0){
					new_proof[i]=proof[i]+upk_u[L-i-2]*temp;
				}
				else{
					new_proof[i]=proof[i]-upk_u[L-i-2]*temp;
				}
			}
			else{
				if(updateindex_binary[i] == 0 && index_binary[i]==1)
					new_proof[i]=proof[i]+g1*temp;
				else if(updateindex_binary[i] == 1 && index_binary[i]==0)
					new_proof[i]=proof[i]-g1*temp;
				else if(updateindex_binary[i] == 0 && index_binary[i]==0)
					new_proof[i]=proof[i]+g1*temp;
				else
					new_proof[i]=proof[i]-g1*temp;
			}
		
		}
	}
	
	return new_proof;
}
