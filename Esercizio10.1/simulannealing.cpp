#include "simulannealing.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <algorithm>
using namespace std;

SimulAnnealing :: SimulAnnealing(): Rnd(nullptr), Nx(0), d(0), L(0), T(1) {
	FinalPath = vector<unsigned int>();
}

SimulAnnealing :: ~SimulAnnealing() {}

SimulAnnealing :: SimulAnnealing(Random* r) {
	Rnd = r;
	Nx = 0;
	d = 0;
	L = 0;
	T = 1;
}

SimulAnnealing :: SimulAnnealing(int _Nx, int _d) {
	Rnd = nullptr;
	if (_Nx < 0) {
		cerr << "Error [SimulAnnealing :: SimulAnnealing(int, int, int)]: input number of points < 0" << endl;
		exit(-1);
	}
	if (_d < 0) {
		cerr << "Error [SimulAnnealing :: SimulAnnealing(int, int, int)]: dimension < 0" << endl;
		exit(-1);
	}
	Nx = (unsigned int)_Nx;
	d = (unsigned int)_d;
	
	for (unsigned int i = 0; i < Nx; i++) {
		X.emplace_back(Vector(d));
	}
	Path = vector<unsigned int>(Nx, 0);
	
	L = 0;
	//FinalPath = vector<unsigned int>();
}

SimulAnnealing :: SimulAnnealing(Random* r, int _Nx, int _d): SimulAnnealing(_Nx, _d) {
	Rnd = r;
}

void SimulAnnealing :: SetParameters(int _Nx, int _d) {
	if (_Nx < 0) {
		cerr << "Error [SimulAnnealing :: SetParameters(int, int)]: input number of points < 0" << endl;
		exit(-1);
	}
	if (_d < 0) {
		cerr << "Error [SimulAnnealing :: SetParameters(int, int)]: dimension < 0" << endl;
		exit(-1);
	}
	
	if (Nx != abs(_Nx)) {
		for (unsigned int i = 0; i < Nx; i++) {
			X[i] = Vector((unsigned int)_d); 
		}
		X.resize(_Nx, Vector((unsigned int)_d));
	}
	
	Nx = _Nx;
	d = _d;
	if (Path.size() != Nx) {
		Path.resize(Nx, 0);
	}
	
	FinalPath = vector<unsigned int>();
}

void SimulAnnealing :: SetNorm(unsigned int p) {
	if (p == 1) Norm = 1;
	else if (p == 2) Norm = 2;
	else {
		cerr << "Error [SimulAnnealing :: SetNorm(unsigned int)]: invalid norm" << endl;
		exit(-1);
	}
	return;
}

void SimulAnnealing :: SetT(double _T) {
	if (_T <= 0) {
		cerr << "Error [SimulAnnealing :: SetT(double)]: invalid input temperature" << endl;
		exit(-1);
	}
	else T = _T;
	return;
}

void SimulAnnealing :: SetCycles(int N) {
	if (N <= 0) {
		cerr << "Error [SimulAnnealing :: SetCycles(int)]: invalid input number of cycles" << endl;
		exit(-1);
	}
	else Cycles = N;
	return;
}

// ***** // ***** // ***** // ***** // *****   POINT GENERATION   ***** // ***** // ***** // ***** // ***** //

void SimulAnnealing :: RandomSphSurface(double R) { // if d = 2 reduces to the random generation on a circumference with radius R
	vector<double> angles(d-1, 0);
	for (unsigned int n = 0; n < Nx; n++) {
		for (unsigned int i = 0; i < d-2; i++) {
			angles[i] = Rnd->Uniform(0, M_PI);
		}
		angles[d-2] = Rnd->Uniform(-M_PI, M_PI); // Last angle (periodicity 2pi)
		X[n].Set(R, angles);
	}
	return;	
}

void SimulAnnealing :: RandomCube(double l) { // if d = 2 reduces to the random generation on a circumference with radius R	
	double half_l = l/2.;
	for (unsigned int n = 0; n < Nx; n++) {
		for (unsigned int i = 0; i < d; i++) {
			X[n][i] = Rnd->Uniform(-half_l, half_l);
		}
	}
	return;	
}

void SimulAnnealing :: GenPath() {
	unsigned int r1, r2, t;
	for (unsigned int i = 0; i < Nx; i++) Path[i] = i;
	for (unsigned int i = 0; i < (unsigned int)Nx/2; i++) {
		r1 = (unsigned int) Rnd->Uniform(1, Nx);
		r2 = (unsigned int) Rnd->Uniform(1, Nx);
		t = Path[r1];
		Path[r1] = Path[r2];
		Path[r2] = t;
	}
	Measure(); // First value of L (ensures that L is initialized)
	return;
}


// ***** // ***** // ***** // ***** // *****   SET METHODS   ***** // ***** // ***** // ***** // ***** //

void SimulAnnealing :: Measure() {
	if (Norm == 1) L1();
	else if (Norm == 2) L2();
	else {
		cerr << "Error [SimulAnnealing :: Measure()]: invalid selection of a L^p norm" << endl;
		exit(-1);
	}
	return; 
}

void SimulAnnealing :: L1() {
	L = 0;
	for (unsigned int i = 0; i < Nx-1; i++) {
		L += X[Path[i+1]].Distance(X[Path[i]]);
	}
	L += X[Path[Nx-1]].Distance(X[Path[0]]);
	return;
}

void SimulAnnealing :: L2() {
	L = 0;
	for (unsigned int i = 0; i < Nx-1; i++) {
		L += X[Path[i+1]].Distance2(X[Path[i]]);
	}
	L += X[Path[Nx-1]].Distance2(X[Path[0]]);
	return;
}

double SimulAnnealing :: Measure(vector<unsigned int>& NewPath) {
	if (Norm == 1) return L1(NewPath);
	else if (Norm == 2) return L2(NewPath);
	else {
		cerr << "Error [SimulAnnealing :: Measure(vector<unsigned int>&)]: invalid selection of a L^p norm" << endl;
		exit(-1);
	}
}

double SimulAnnealing :: L1(vector<unsigned int>& NewPath) {
	double NewL = 0;
	for (unsigned int i = 0; i < Nx-1; i++) {
		NewL += X[NewPath[i+1]].Distance(X[NewPath[i]]);
	}
	NewL += X[NewPath[Nx-1]].Distance(X[NewPath[0]]);
	return NewL;
}

double SimulAnnealing :: L2(vector<unsigned int>& NewPath) {
	double NewL = 0;
	for (unsigned int i = 0; i < Nx-1; i++) {
		NewL += X[NewPath[i+1]].Distance2(X[NewPath[i]]);
	}
	NewL += X[NewPath[Nx-1]].Distance2(X[NewPath[0]]);
	return NewL;
}




// ***** // ***** // ***** // ***** // *****   MUTATION ALGORITHMS   ***** // ***** // ***** // ***** // ***** //
vector<double> SimulAnnealing :: Cycle(double threshold) {
	vector<double> L_cycle;
	for (unsigned int i = 0; i < Cycles; i++) {
		Mutate(threshold);
		L_cycle.push_back(L);
	}
	return L_cycle;
}


void SimulAnnealing :: Mutate(double threshold) { // Probability of each simple mutation (Swap, Permute, Shift, Invert)
	double r;
	int pos, blocks, returned = 0;
	
	vector<unsigned int> NewPath = Path;
	for (unsigned int Nmut = 0; Nmut < 1; Nmut++) {
		r = Rnd->Uniform(0, 1);
		if (r < threshold) { // Simple Swap
			returned = Swap(NewPath, (int)Rnd->Uniform(1, Nx), (int)Rnd->Uniform(1, Nx));
		}
		else if (r < 2.*threshold) { // Permutation
			pos = (int)Rnd->Uniform(1, (int)Nx-1);
			returned = Permute(NewPath, (int)Rnd->Uniform(1, (int)Nx-1-pos)/2, pos);
		}
		else if (r < 3.*threshold) { // Shift
			if (Nx >= 3) {
				pos = (int)Rnd->Uniform(1, (int)Nx);
				blocks = (int)Rnd->Uniform(1,(int) Nx-pos);
				returned = Shift(NewPath, pos, blocks, (int)Rnd->Uniform(-(pos-1), (int)Nx+1-(pos+blocks)));
			}
		}
		else if (r < 4.*threshold) { // Inversion
			pos = (int)Rnd->Uniform(1, (int)Nx-1);
			returned = Invert(NewPath, pos, (int)Rnd->Uniform(1, (int)Nx-pos));
		}
		/*
		else if (r < 10.*threshold) { // Crossing Over (probability ~60%)
			CrossingOver(NewPath, (unsigned int)Rnd->Uniform(first, last), (unsigned int)Rnd->Uniform(1, Nx-1));
		}
		*/
		if (returned == -1) {
			cerr << "Warning [SimulAnnealing :: Mutate]: invalid mutation attempt" << endl;
			exit(-1);
		}
	}
	
	double NewL = Measure(NewPath);
	double p = exp((L-NewL)/T); 
	
	if (p >= Rnd->Uniform()) {
		Path = NewPath;
		L = NewL;
	}
	return;
}

int SimulAnnealing :: Swap(vector<unsigned int>& NewPath, unsigned int i, unsigned int j) {
	if (i >= Nx || j >= Nx) {
		cerr << "Error [SimulAnnealing :: Swap("<<i<<", "<<j<<")]: indices out of range" << endl << endl;
		exit(-1);
	}
	if (i == 0 || j == 0) {
		cerr << "Warning [SimulAnnealing :: Swap("<<i<<", "<<j<<")]: invalid swap." << endl << endl;
		return -1;
	}
	unsigned int t = NewPath[i];
	NewPath[i] = NewPath[j];
	NewPath[j] = t;
	return 0;
}

int SimulAnnealing :: Permute(vector<unsigned int>& NewPath, int n_block, int i, int j) {
	if (i <= 0 || (j <= 0 && j != -1) || (j == -1 && abs(i+2*n_block) > Nx) || (j > 0 && (abs((i>j?i:j)+n_block) > Nx || abs(i-j)<n_block))) {
		cerr << "Warning [SimulAnnealing :: Permute("<<n_block<<", "<<i;
		if (j!=-1) cerr << ", " << j;
		cerr <<")]: invalid permutation" << endl << endl;
		return -1;
	}
	else if (j == -1) j = i+n_block; // Contiguous block permutation (by default, since j = -1 is a default parameter)
	unsigned int t;
	for (int n = 0; n < n_block; n++) {
		t = NewPath[i+n];
		NewPath[i+n] = NewPath[j+n];
		NewPath[j+n] = t;
	}
	return 0;
}

// Shift n_shift indices by n_places places (>< 0), starting from the i-th index
int SimulAnnealing :: Shift(vector<unsigned int>& NewPath, int i, int n_block, int n_shift) { 
	if (i <= 0 || abs(i) >= Nx || abs(i+n_shift+n_block) > Nx || i+n_shift <= 0) {
		cerr << "Warning [SimulAnnealing :: Shift("<<i<<", "<<n_block<<", "<<n_shift<<")]: invalid shift" << endl << endl;
		return -1;
	}
	else {
		unsigned int t;
		if (n_shift >= 0) {
			for (int k = n_block-1; k >= 0; k--) {
				for (int j = i; j < i+n_shift; j++) {
					t = NewPath[j+k];
					NewPath[j+k] = NewPath[j+k+1];
					NewPath[j+k+1] = t;
				}
			}
		}
		else {
			for (int k = 0; k <= n_block-1; k++) {
				for (int j = i; j > i+n_shift; j--) {
					t = NewPath[j+k];
					NewPath[j+k] = NewPath[j+k-1];
					NewPath[j+k-1] = t;
				}
			}
		}
	}
	return 0;
}

int SimulAnnealing :: Invert(vector<unsigned int>& NewPath, unsigned int i, unsigned int n_inv) {
	if (i == 0 || i+n_inv > Nx) {
		cerr << "Warning [SimulAnnealing :: Invert("<<i<<", "<<n_inv<<")]: invalid inversion" << endl << endl;
		return -1;
	}
	else {
		unsigned int t;
		for (unsigned int k = 0; k < (unsigned int)(n_inv/2); k++) {
			t = NewPath[i+k];
			NewPath[i+k] = NewPath[i+n_inv-1-k];
			NewPath[i+n_inv-1-k] = t;
		}
	}
	return 0;
}


/*
void SimulAnnealing :: CrossingOver(unsigned int p, unsigned int q, unsigned int cut) { // cut = first index of the mutable section
	if (p >= Npop || q >= Npop) {
		cerr << "Error [SimulAnnealing :: CrossingOver(unsigned int, unsigned int, unsigned int)]: population index out of range" << endl;
		exit(-1);
	}
	if (cut == 0 || cut >= Nx) {
		cerr << "Warning [SimulAnnealing :: CrossingOver("<<p<<", "<<q<<", "<<cut<<")]: invalid crossing over" << endl << endl;
		return;
	}
	if (cut == Nx-1) return; // No valid crossing over with only one mutable index can be performed
	else {
		
		vector<unsigned int> exchangep, exchangeq;
		for (unsigned int i = cut; i < Nx; i++) {
			exchangep.push_back(Population[p][i]);
			exchangeq.push_back(Population[q][i]);
		}
		
		sort(exchangep.begin(), exchangep.end());
		sort(exchangeq.begin(), exchangeq.end());
		
		// Search the p-th individual's indices in the q-th individual
		int found, c = 0; // Found index, #found (count)
		for (unsigned int find = 1; find < Nx; find++) {
			found = BinarySearch(exchangep, c, exchangep.size()-1, Population[q][find]);
			if (found != -1) {
				Shift(exchangep, found, c);
				c++;
			}
			if (abs(c) >= exchangep.size()) break;
		}
		//PrintV(exchangep);

		// Search the q-th individual's indices in the p-th individual
		c = 0; 
		for (unsigned int find = 1; find < Nx; find++) {
			found = BinarySearch(exchangeq, c, exchangeq.size()-1, Population[p][find]);
			if (found != -1) {
				Shift(exchangeq, (unsigned int)found, c);
				c++;
			}
			if (abs(c) >= exchangeq.size()) break;
		} 
		
		for (unsigned int i = cut; i < Nx; i++) {
			Population[p][i] = exchangep[(int)(i-cut)];
			Population[q][i] = exchangeq[(int)(i-cut)];
		}
	}
	return;
}

*/

// ***** // ***** // ***** // ***** // *****   AUXILIARY ALGORITHMS   ***** // ***** // ***** // ***** // ***** //
/*
// Shift n_shift indices by n_places places, starting from the i-th index
void SimulAnnealing :: Shift(vector<unsigned int>& V, int i, int i_new) const { 
	if (abs(i>i_new?i:i_new) > V.size() || (i<i_new?i:i_new) < 0) {
		cerr << "Warning [SimulAnnealing :: Shift(vector<unsigned int>, "<<i<<", "<<i_new<<")]: invalid shift" << endl << endl;
		return;
	}
	else {
		unsigned int t;
		if (i_new-i >= 0) {
			for (int j = i; j < i_new; j++) {
				t = V[j];
				V[j] = V[j+1];
				V[j+1] = t;
			}
		}
		else {
			for (int j = i; j > i_new; j--) {
				t = V[j];
				V[j] = V[j-1];
				V[j-1] = t;
			}
		}
	}
	return;
}

int SimulAnnealing :: BinarySearch(vector<unsigned int> V, unsigned int first, unsigned int last, unsigned int find) const {
	if (first == last) {
		if (V[first] == find) return first;
		else return -1; // Not found;
	}
	else {
		if (V[first] > find || V[last] < find) return -1; // Not found
		else if (V[first] == find) return first;
		else if (V[last] == find) return last;
		else {
			unsigned int m;
			do {
				m = (unsigned int)(first+last)/2;
				if (V[m] == find) return m;
				else if (V[m] < find) {
					first = m;
				}
				else if (V[m] > find) {
					last = m;
				}
			} while (first < last-1);
			// first, last are contiguous indices
			if (V[first] == find) return first;
			else if (V[last] == find) return last;
			else return -1; // Not found
		}
	}
}
*/

// ***** // ***** // ***** // ***** // *****   SCREEN OUTPUT METHODS   ***** // ***** // ***** // ***** // ***** //

void SimulAnnealing :: PrintPath() const {
	cout << "Path = " << "(";
	for (unsigned int n = 0; n < Nx-1; n++) cout << Path[n] << ", ";
	cout << Path[Nx-1] << ")" << endl;
	return;
}
void SimulAnnealing :: PrintPath(vector<unsigned int> NewPath) const {
	cout << "Path = " << "(";
	for (unsigned int n = 0; n < Nx-1; n++) cout << NewPath[n] << ", ";
	cout << NewPath[Nx-1] << ")" << endl;
	return;
}

void SimulAnnealing :: PrintX() const {
	cout << "*** X array ***" << endl;
	for (unsigned int n = 0; n < Nx; n++) {
		cout << "[" << n << "] -> ";
		X[n].Print();
	}
	cout << endl;
}


/*

void SimulAnnealing :: PrintV(vector<unsigned int> V) const {
	cout << "Vector:  ";
	for (unsigned int i = 0; i < V.size(); i++) {
		cout << V[i] << " ";
	}
	cout << endl;
}

void SimulAnnealing :: PrintBest() const {
	cout << "> Best Individual: "; PrintIndividual(Order[0]); // Order[0] -> the individuals are not ordered by Sort()
	cout << "  L = " << L[0] << endl << endl; // [0] -> the array L |>IS<| ordered by Sort() 
}

void SimulAnnealing :: PrintBestL() const {
	cout << "> Best L = " << L[0] << endl << endl; // [0] -> the array L |>IS<| ordered by Sort() 
}
*/

void SimulAnnealing :: LoadPoints(const string filename) {
	ifstream in(filename);
	if (in.is_open()) {
		string buff;
		unsigned int _d, _Nx;
		in >> buff >> _d >> buff >> _Nx;
		SetParameters(_Nx, _d);
		Vector v(d);
		double t;
		for (unsigned int i = 0; i < Nx; i++) {
			for (unsigned int j = 0; j < d; j++) {
				in >> t;
				v[j] = t;
			}
			X[i] = v;
		}	
	}
	else {
		cerr << "Error [SimulAnnealing :: LoadPoints(const string)]: cannot open input file" << endl;
		exit(-1);
	}
	return;
}

void SimulAnnealing :: OutputPoints(const string filename) const {
	ofstream out(filename);
	if (out.is_open()) {
		out << "#SpaceDimension: " << d << endl << "#NumberOfPoints: " << Nx;
		for (unsigned int i = 0; i < Nx; i++) {
			out << endl;
			for (unsigned int j = 0; j < d-1; j++) {
				out << X[i][j] << " ";
			}
			out << X[i][d-1];
		}
		out.close();
	}
	else {
		cerr << "Error [SimulAnnealing :: OutputPoints(const string)]: cannot open input file" << endl;
		exit(-1);
	}
	return;
}

void SimulAnnealing :: LoadPathIndices(const string filename) {
	cout << "> Loading path indices from file <"+filename+">..."<< endl;
	ifstream in(filename);
	if (in.is_open()) {
		string buff;
		unsigned int _Nx, _d;
		in >> buff >> _d;
		in >> buff >> _Nx;
		in >> buff >> buff;
		SetParameters(_Nx, _d);
		
		for (unsigned int i = 0; i < Nx; i++) in >> Path[i];
		in.close();
	}
	else {
		cerr << "Error [SimulAnnealing :: LoadPathIndices(const string)]: cannot open input file" << endl;
		exit(-1);
	}
	return;
}
void SimulAnnealing :: OutputPathIndices(const string filename) const {
	cout << "> Saving path indices to file <"+filename+">..."<< endl;
	ofstream out(filename);
	if (out.is_open()) {
		out << "#SpaceDimension: " << d << endl;
		out << "#NumberOfPoints: " << Nx << endl;
		out << "#L: " << L << endl;
		
		for (unsigned int i = 0; i < Nx-1; i++) out << Path[i] << endl;
		out << Path[Nx-1];
		out.close();
	}
	else {
		cerr << "Error [SimulAnnealing :: OutputPathIndices(const string)]: cannot open input file" << endl;
		exit(-1);
	}
	return;
}

void SimulAnnealing :: OutputPathCoord(const string filename) const {
	cout << "> Saving path coordinates to file <"+filename+">..."<< endl;
	ofstream out(filename);
	if (out.is_open()) {
		out << "#SpaceDimension: " << d << endl;
		out << "#NumberOfPoints: " << Nx << endl;
		out << "#L: " << L << endl;
		for (unsigned int i = 0; i < Nx-1; i++) {
			out << i+1 << "   ";
			for (unsigned int k = 0; k < d-1; k++) out << X[Path[i]][k] << " ";
			out << X[Path[i]][d-1] << endl;
		}	
		out << Nx << "   ";
		for (unsigned int k = 0; k < d-1; k++) out << X[Path[Nx-1]][k] << " ";
		out << X[Path[Nx-1]][d-1];
		out.close();
	}
	else {
		cerr << "Error [SimulAnnealing :: OutputPathCoord(const string)]: cannot open input file" << endl;
		exit(-1);
	}
	return;
}
