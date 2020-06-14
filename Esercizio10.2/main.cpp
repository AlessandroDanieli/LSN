#include "random.h"
#include "genetic.h"
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include "/home/alessandro/Scrivania/LSN/openmpi-4.0.3/ompi/include/mpi.h"
using namespace std;

void Input(const string filename, Genetic& G, char& c, unsigned int& Ngen, unsigned int& Nmigr) {
	ifstream Input;
	Input.open(filename);
	unsigned int Npop, Nx, d;
	double r;
	Input >> d; 			Input.ignore(120, '\n');	 	
	Input >> Nx; 			Input.ignore(120, '\n');
	Input >> Npop; 			Input.ignore(120, '\n');
	G.SetParameters(Npop, Nx, d);
	Input >> Ngen;			Input.ignore(120, '\n');
	Input >> Nmigr;			Input.ignore(120, '\n');	
	Input >> c; 			Input.ignore(120, '\n'); 		c = tolower(c);
	Input >> d;				Input.ignore(120, '\n'); 		G.SetNorm(d);
	Input >> r;				Input.ignore(120, '\n');
	
	if (Ngen % Nmigr != 0) {
		cerr << "Error [main :: Input(const string, Genetic&, char&, unsigned int&, unsigned int&)]: Ngen % Nmigr != 0" << endl;
		exit(-1);
	}
	
	if (c == 'c') {
		G.RandomSphSurface(r);
	}
	else if (c == 's') {
		G.RandomCube(r);
	}
	G.GenPopulation();
	
	c = toupper(c);
	Input.close();
	return;
}



void OutputL(const string filename, vector<double> L, vector<double> MeanL) {
	//cout << "> Saving <L> to file <"+filename+">..."<< endl;
	ofstream out(filename);
	if (out.is_open()) {
		out << "#Generations: " << MeanL.size()<< endl;
		out << "#Individuals: " << L.size()<< endl;
		for (unsigned int i = 0; i < MeanL.size(); i++) {
			out << MeanL[i] << endl; 
		}
		out << "#" << endl;
		for (unsigned int i = 0; i < L.size()-1; i++) {
			out << L[i] << endl; 
		}
		out << L[L.size()-1];
		out.close();
	}
	else {
		cerr << "Error [OutputMeanL(const string, vector<double>, vector<double>)]: cannot open input file" << endl;
		exit(-1);
	}
	return;
}

void SetRandGen(Random* Rnd, int rank) {
	int p1, p2;
	int seed[4];

	//******* Input: Primes *******//
	ifstream Primes("Primes");
	if (Primes.is_open()) {
		for (int i = 0; i <= rank; i++) Primes >> p1 >> p2;
	}
	else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();
	//******* Input: seed.in *******//
	ifstream Seed("seed.in");
	string property;
	Seed >> property;
	if (Seed.is_open()) {
		Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
		//cout << seed[0] << " " << seed[3] << endl;
		Rnd->SetRandom(seed, p1, p2);
	}
	else cerr << "PROBLEM: Unable to open seed.in" << endl;
	Seed.close();
	return;
}

void Swap(int* V, int p, int q) {
	int temp = V[p];
	V[p] = V[q];
	V[q] = temp;
	return;
}

void OutputPath(const string filename, Genetic& G, int* Path, double L) { // Output the best path
	int Nx = G.GetNx();
	int d = G.GetD();
	//cout << "> Saving to file <"+filename+">..."<< endl;
	ofstream out(filename);
	if (out.is_open()) {
		out << "#SpaceDimension: " << d << endl << "#NumberOfPoints: " << Nx << endl;
		out << "#L: " << L << endl;
		for (unsigned int i = 0; i < Nx-1; i++) {
			out << i+1 << "   ";
			for (unsigned int k = 0; k < d-1; k++) out << G.GetX()[Path[i]][k] << " ";
			out << G.GetX()[Path[i]][d-1] << endl;
		}	
		out << Nx << "   ";
		for (unsigned int k = 0; k < d-1; k++) out << G.GetX()[Path[Nx-1]][k] << " ";
		out << G.GetX()[Path[Nx-1]][d-1];
	}
	else {
		cerr << "Error [Genetic :: OutputPath(const string)]: cannot open input file" << endl;
		exit(-1);
	}
	return;
}


/*
	int n = 8;
			
	int* recv = new int[size*n];
	for(int i = 0; i < size*n; i++) {
		recv[i]=i;
	}
	
	int* send = new int[n];
	for (int k=0; k<n; k++){
		send[k]=10*rank+k;
	}
	
	MPI_Gather(send, n, MPI_INTEGER, recv, n, MPI_INTEGER, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		for(int i = 0; i < size*n; i++) {
			cout<< recv[i] <<endl;
		}
	}// <<" "<<recv[1] <<" " <<recv[2] <<endl;
	
	MPI_Finalize();
	return 0;
	*/
	
/*
for (int i = 0; i < size*Nx; i++) {
	cout << PathArray[i] << " ";
} cout << endl << endl;
*/
//MPI_Bcast(PathArray, Nx*size, MPI_INTEGER, 0, MPI_COMM_WORLD);


int main(int argc, char* argv[]) {
	Random Rnd, Exc;
	Exc.RandGenSetup();
	Genetic G;
	
	int res = system("rm -rf seed.out");
	bool exec = true;
	
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status* stat = new MPI_Status[size];
	MPI_Request* req = new MPI_Request[size];
	double t1 = MPI_Wtime();
	
	Rnd.RandGenSetup();  
	
	G.SetRandGen(&Rnd); // Same seed -> Same points
	unsigned int Ngen;
	unsigned int Nmigr;
	char c;
	Input("input.dat", G, c, Ngen, Nmigr); // Initialization and point generation
	
	vector<double> L, MeanL;
	
	SetRandGen(&Rnd, rank);  // Different seeds -> Same points, but different path exploration
	
	int tag = 1;
	int Nx = G.GetNx();
	int* PathArray = new int[Nx*size];
	int *RecvSinglePath = new int[Nx];
	int *SinglePath = new int[Nx];
	
	for (int cycle = 0; cycle < Ngen; cycle++) {
		if (cycle > 0 && cycle % Nmigr == 0) { // Exchange between nodes
			
			MPI_Gather(G.GetBest(), Nx, MPI_INTEGER, PathArray, Nx, MPI_INTEGER, 0, MPI_COMM_WORLD);
			
			// Random exchanges
			if (rank == 0) {
				for (int i = 0; i < size; i++) {
					if (Exc.Uniform(0, 1) < 0.5) {
						for (int k = 0; k < Nx; k++) {
							Swap(PathArray, i*Nx+k, ((i+1)%size)*Nx+k);
							//cout << i*Nx+k << " " <<  ((i+1)%size)*Nx+k << endl;
						}
					}
				}
				for (int i = 0; i < size; i+=2) {
					if (Exc.Uniform(0, 1) < 0.5) {
						for (int k = 0; k < Nx; k++) {
							Swap(PathArray, i*Nx+k, ((i+1)%size)*Nx+k);
							//cout << i*Nx+k << " " <<  ((i+1)%size)*Nx+k << endl;
						}
					}
				}
				/*
				for (int i = 0; i < size*Nx; i++) {
					cout << PathArray[i] << " ";
				} cout << endl << endl;
				*/
				// SETTING (rank = 0)'S OWN PATH
				for (int i = 0; i < Nx; i++) {
					SinglePath[i] = PathArray[i];
				}
				G.SetBest(SinglePath);
				
				// SENDING THE OTHER PATHS TO THE rank!=0 NODES
				for (int r = 1; r < size; r++) {
					for (int i = 0; i < Nx; i++) {
						SinglePath[i] = PathArray[i+r*Nx];
					}
					
					MPI_Isend(SinglePath, Nx, MPI_INTEGER, r, tag, MPI_COMM_WORLD, &req[r]);
				}
			}
			
			if (rank != 0) {
				MPI_Recv(RecvSinglePath, Nx, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, &stat[rank]);
				/*
				cout << "Rank = " << rank << endl;
				for (int i = 0; i < Nx; i++) {
					cout << RecvSinglePath[i] << " ";
				}cout << endl;
				*/
				G.SetBest(RecvSinglePath);
			}	
		}
			
		G.Measure();
		G.Sort();
		G.Select();
		//G.PrintBestL();
	}
	G.SetBestPath();

	double* BestL = new double[size];
	double ActualBestL = G.GetBestL();
	MPI_Gather(&ActualBestL, 1, MPI_DOUBLE_PRECISION, BestL, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
	MPI_Gather(G.GetBest(), Nx, MPI_INTEGER, PathArray, Nx, MPI_INTEGER, 0, MPI_COMM_WORLD);
	
	if (rank == 0) {
		for (int r = 0; r < size; r++) {
			for (int i = 0; i < Nx; i++) {
				SinglePath[i] = PathArray[i+r*Nx];
			}
			OutputPath("Path"+to_string(G.GetNorm())+string(1, c)+"_rank"+to_string(r)+".out", G, SinglePath, BestL[r]);
			//OutputL("L"+to_string(G.GetNorm())+string(1, c)+"_rank"+to_string(rank)+".out", L, MeanL);
		}
	}	
	
	double t2 = MPI_Wtime();
	cout << "[Rank " << rank << "] Elapsed time = " << t2-t1 << " s " << endl;
	MPI_Finalize();	
	
	
	
	Rnd.SaveSeed();
	return 0;
	
}		

	//Rnd.SaveSeed();
	
	
	// Input: seed.out
	/*
	ifstream Seed;
	Seed.open("seed.out");
	if (Seed.is_open()) {
		if (Seed.is_open()) {
			int seed[4];
			Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
			Rnd.SetRandom(seed, p1, p2);
		}
		Seed.close();
	}
	*/
	
	
	
	//const int n = 300000;
	
	
	/*
	int size, rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status* stat = new MPI_Status[4];
	MPI_Request* req = new MPI_Request[4];
	
	int* imesg1 = new int[n]; 
	int* imesg2 = new int[n];
	int* imesg3 = new int[n]; 
	int* imesg4 = new int[n];
	int itag=1; 
	int itag2=2;
	for(int i=0; i<n; i++){
		imesg1[i] = rank; 
		imesg2[i] = rank+1;
		imesg3[i] = rank+2; 
		imesg4[i] = rank+3;
	}
	if (rank==0) {
		MPI_Isend(&imesg1[0], n, MPI_INTEGER, 1, 00, MPI_COMM_WORLD, &req[0]);
		MPI_Recv(&imesg2[0], n, MPI_INTEGER, 1, 01, MPI_COMM_WORLD, &stat[0]);
		//MPI_Wait(&req, &stat1);
		cout << rank << endl;
		cout << "messaggio = " << imesg2[0] << endl;
	}
	else if (rank==1) {
		MPI_Isend(&imesg2[0], n, MPI_INTEGER, 0, 01, MPI_COMM_WORLD, &req[1]);
		MPI_Recv(&imesg1[0], n, MPI_INTEGER, 0, 00, MPI_COMM_WORLD, &stat[1]);
		cout << rank << endl;
		cout << "messaggio = " << imesg1[0] << endl;
		//MPI_Send(&imesg2[0], n, MPI_INTEGER, 1, itag2, MPI_COMM_WORLD);
	}
	else if (rank==2) {
		MPI_Isend(&imesg3[0], n, MPI_INTEGER, 3, 02, MPI_COMM_WORLD, &req[2]);
		
		MPI_Recv(&imesg4[0], n, MPI_INTEGER, 3, 03, MPI_COMM_WORLD, &stat[2]);
		cout << rank << endl;
		cout << "messaggio = " << imesg4[0] << endl;
		//MPI_Send(&imesg2[0], n, MPI_INTEGER, 3, itag2, MPI_COMM_WORLD);
	}
	else if (rank==3) {
		MPI_Isend(&imesg4[0], n, MPI_INTEGER, 2, 03, MPI_COMM_WORLD, &req[3]);
		//MPI_Send(&imesg[0], n, MPI_INTEGER, 2, itag2, MPI_COMM_WORLD);
		MPI_Recv(&imesg3[0], n, MPI_INTEGER, 2, 02, MPI_COMM_WORLD, &stat[3]);
		//MPI_Wait(&req, &stat1);
		cout << rank << endl;
		cout << "messaggio = " << imesg3[0] << endl;
	}
	MPI_Finalize();
	*/


/*
	int size, rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status stat1, stat2;
	int* imesg = new int[n]; 
	int* imesg2 = new int[n];
	int itag=1; 
	int itag2=2;
	for(int i=0; i<n; i++) {
		imesg[i]=rank; 
		imesg2[i]=rank+1;
	}
	if(rank==1){
		MPI_Send(&imesg[0], n, MPI_INTEGER, 0, itag, MPI_COMM_WORLD);
		MPI_Recv(&imesg2[0], n, MPI_INTEGER, 0, itag2, MPI_COMM_WORLD, &stat2);
		cout<<"messaggio = "<<imesg2[0]<<endl;
	}
	else if(rank==0){
		MPI_Send(&imesg2[0], n, MPI_INTEGER, 1, itag2, MPI_COMM_WORLD);
		MPI_Recv(&imesg[0], n, MPI_INTEGER, 1, itag, MPI_COMM_WORLD, &stat1);
		cout<<"messaggio = "<<imesg[0]<<endl;
	}
	MPI_Finalize();
	
	*/
