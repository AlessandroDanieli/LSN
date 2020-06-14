#include "genetic.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <algorithm>
using namespace std;

Genetic :: Genetic(): Rnd(nullptr), Npop(0), Nx(0), d(0) {
	FinalPath = vector<unsigned int>();
}

Genetic :: ~Genetic() {}

Genetic :: Genetic(int _Npop, int _Nx, int _d) {
	Rnd = nullptr;
	if (_Npop < 0) {
		cerr << "Error [Genetic :: Genetic(int, int, int)]: input size of population < 0" << endl;
		exit(-1);
	}
	if (_Nx < 0) {
		cerr << "Error [Genetic :: Genetic(int, int, int)]: input number of points < 0" << endl;
		exit(-1);
	}
	if (_d < 0) {
		cerr << "Error [Genetic :: Genetic(int, int, int)]: dimension < 0" << endl;
		exit(-1);
	}
	
	for (unsigned int i = 0; i < abs(_Npop); i++) {
		Population.push_back(vector<unsigned int>(_Nx, 0));
		L.push_back(0);
		Order.push_back(i);
	}
	for (unsigned int i = 0; i < abs(_Nx); i++) {
		X.emplace_back(Vector((unsigned int)_d)); 
	}
	FinalPath = vector<unsigned int>();
	Npop = _Npop;
	Nx = _Nx;
	d = _d;
}

Genetic :: Genetic(Random* r, int _Npop, int _Nx, int _d): Genetic(_Npop, _Nx, _d) {
	Rnd = r;
}

void Genetic :: GenPopulation() {
	unsigned int r1, r2;
	for (unsigned int p = 0; p < Npop; p++) {
		for (unsigned int i = 0; i < Nx; i++) Population[p][i] = i;
		for (unsigned int i = 0; i < (unsigned int)Nx/2; i++) {
			r1 = (unsigned int) Rnd->Uniform(1, Nx);
			r2 = (unsigned int) Rnd->Uniform(1, Nx);
			Swap(p, r1, r2);
		}
	}
	return;
}

// ***** // ***** // ***** // ***** // *****   POINT GENERATION   ***** // ***** // ***** // ***** // ***** //

void Genetic :: RandomSphSurface(double R) { // if d = 2 reduces to the random generation on a circumference with radius R
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

void Genetic :: RandomCube(double l) { // if d = 2 reduces to the random generation on a circumference with radius R	
	double half_l = l/2.;
	for (unsigned int n = 0; n < Nx; n++) {
		for (unsigned int i = 0; i < d; i++) {
			X[n][i] = Rnd->Uniform(-half_l, half_l);
		}
	}
	return;	
}

// ***** // ***** // ***** // ***** // *****   SET METHODS   ***** // ***** // ***** // ***** // ***** //

void Genetic :: SetParameters(int _Npop, int _Nx, int _d) {
	if (_Npop < 0) {
		cerr << "Error [Genetic :: SetParameters(int, int, int)]: input size of population < 0" << endl;
		exit(-1);
	}
	if (_Nx < 0) {
		cerr << "Error [Genetic :: SetParameters(int, int, int)]: input number of points < 0" << endl;
		exit(-1);
	}
	if (_d < 0) {
		cerr << "Error [Genetic :: SetParameters(int, int, int)]: dimension < 0" << endl;
		exit(-1);
	}
	FinalPath = vector<unsigned int>();
	if (Npop != abs(_Npop)) {
		Population.resize(_Npop);
		L.resize(_Npop);
		Order.resize(_Npop);
		for (unsigned int i = 0; i < abs(_Npop); i++) {
			Population[i] = vector<unsigned int>(_Nx, 0);
			L[i] = 0;
			Order[i] = i;
		}
		Npop = _Npop;
	}
	if (Nx != abs(_Nx)) {
		for (unsigned int i = 0; i < Nx; i++) {
			X[i] = Vector((unsigned int)_d); 
		}
		X.resize(_Nx, Vector((unsigned int)_d));
	}
	
	Nx = _Nx;
	d = _d;
}

void Genetic :: SetNorm(unsigned int p) {
	if (p == 1) Norm = 1;
	else if (p == 2) Norm = 2;
	else {
		cerr << "Error [Genetic :: SetNorm(unsigned int)]: invalid norm" << endl;
		exit(-1);
	}
	return;
}

void Genetic :: Measure() {
	if (Norm == 1) L1();
	if (Norm == 2) L2();
	return; 
}

void Genetic :: L1() {
	for (unsigned int p = 0; p < Npop; p++) {
		L[p] = 0;
		for (unsigned int i = 0; i < Nx-1; i++) {
			L[p] += X[Population[p][i+1]].Distance(X[Population[p][i]]);
		}
		L[p] += X[Population[p][Nx-1]].Distance(X[Population[p][0]]);
	}
	return;
}

void Genetic :: L2() {
	for (unsigned int p = 0; p < Npop; p++) {
		L[p] = 0;
		for (unsigned int i = 0; i < Nx-1; i++) {
			L[p] += X[Population[p][i+1]].Distance2(X[Population[p][i]]);
		}
		L[p] += X[Population[p][Nx-1]].Distance2(X[Population[p][0]]);
	}
	return;
}

void Genetic :: Sort(int sx, int dx) {
	if (sx == -1 && dx == -1) {
		sx = 0; dx = (int)Npop-1;
	}
	else if (dx <= sx || sx < 0 || dx < 0) return;
	int i = sx, j = dx;
	double m;
	m = (double)(L[sx] + L[dx])/2.;
	
	unsigned int TempIndex;
	double TempL;
	do {
		while(L[i] < m) i++;
		while(L[j] > m) j--;
		
		if (i <= j) {
			TempL = L[i];
			L[i] = L[j];
			L[j] = TempL;
			
		    TempIndex = Order[i];
		    Order[i] = Order[j];
		    Order[j] = TempIndex;
		    
		    i++;
		    j--;
		}
	} while (j >= i);
	
	if (sx < j) Sort(sx, j);
	if (i < dx) Sort(i, dx);
	return;
}


double Genetic :: MeanL() const {
	double ML = 0;
	for (unsigned int i = 0; i < Npop/2; i++) {
		ML += L[i];
	}
	return 2.*ML/(double)(Npop);
}

// ***** // ***** // ***** // ***** // *****   SELECTION ALGORITHMS   ***** // ***** // ***** // ***** // ***** //

void Genetic :: Select() { // Select the best individuals and generate the remaining
	unsigned int F = (unsigned int)(Npop/4);
	unsigned int R = (unsigned int)((int)Npop-(int)F);

	vector<vector<unsigned int>> New_Generation;
	for (unsigned int i = 0; i < F; i++) {
		New_Generation.push_back(Population[Order[i]]);
	}
	for (unsigned int i = 0; i < R; i++) {
		New_Generation.push_back(Population[Order[i]]);
	}
	for (unsigned int i = 0; i < Npop; i++) {
		Order[i] = i;
	}
	//cout << "Select " << endl;
	//PrintOrder();
	Population.clear();
	Population = New_Generation;
	
	Mutate(F, Npop);
}

void Genetic :: SetBestPath() {
	FinalPath = Population[Order[0]];
	return;
}


// ***** // ***** // ***** // ***** // *****   MUTATION ALGORITHMS   ***** // ***** // ***** // ***** // ***** //

void Genetic :: Mutate(int first, int last) {
	double threshold = 0.08; // Probability of each simple mutation (Swap, Permute, Shift, Invert). P(Crossing Over) ~50%
	double r;
	int pos, blocks, returned;
	
	if (first == -1 && last == -1) { // Default: every individual is subject to mutation
		first = 0; last = Npop;
	}
	
	for (int p = first; p < last; p++) {
		r = Rnd->Uniform(0, 1);
		
		if (r < threshold) { // Simple Swap
			returned = Swap(p, (int)Rnd->Uniform(1, Nx), (int)Rnd->Uniform(1, Nx));
		}
		else if (r < 2.*threshold) { // Permutation
			pos = (int)Rnd->Uniform(1, (int)Nx-1);
			returned = Permute(p, (int)Rnd->Uniform(1, (int)Nx-1-pos)/2, pos);
		}
		else if (r < 3.*threshold) { // Shift
			if (Nx >= 3) {
				pos = (int)Rnd->Uniform(1, (int)Nx);
				blocks = (int)Rnd->Uniform(1,(int) Nx-pos);
				returned = Shift(p, pos, blocks, (int)Rnd->Uniform(-(pos-1), (int)Nx+1-(pos+blocks)));
			}
		}
		else if (r < 4.*threshold) { // Inversion
			pos = (int)Rnd->Uniform(1, (int)Nx-1);
			returned = Invert(p, pos, (int)Rnd->Uniform(1, (int)Nx-pos));
		}
		else if (r < 10.*threshold) { // Crossing Over (probability ~60%)
			CrossingOver(p, (unsigned int)Rnd->Uniform(first, last), (unsigned int)Rnd->Uniform(1, Nx-1));
		}
		if (returned == -1) {
			cerr << "Warning [Genetic :: Mutate]: invalid operation at Population Index = " << p << endl;
			exit(-1);
		}
	}
}

int Genetic :: Swap(unsigned int p, unsigned int i, unsigned int j) {
	if (p >= Npop) {
		cerr << "Error [Genetic :: Swap(unsigned int, unsigned int, unsigned int)]: population index out of range" << endl;
		exit(-1);
	}
	if (i >= Nx || j >= Nx) {
		cerr << "Error [Genetic :: Swap("<<p<<", "<<i<<", "<<j<<")]: indices out of range" << endl << endl;
		exit(-1);
	}
	if (i == 0 || j == 0) {
		cerr << "Warning [Genetic :: Swap("<<p<<", "<<i<<", "<<j<<")]: invalid swap." << endl << endl;
		return -1;
	}
	unsigned int t = Population[p][i];
	Population[p][i] = Population[p][j];
	Population[p][j] = t;
	return 0;
}

int Genetic :: Permute(unsigned int p, int n_block, int i, int j) {
	if (p >= Npop) {
		cerr << "Error [Genetic :: Permute(unsigned int, unsigned int, unsigned int)]: population index out of range" << endl;
		exit(-1);
	}
	if (i <= 0 || (j <= 0 && j != -1) || (j == -1 && abs(i+2*n_block) > Nx) || (j > 0 && (abs((i>j?i:j)+n_block) > Nx || abs(i-j)<n_block))) {
		cerr << "Warning [Genetic :: Permute("<<p<<", "<<n_block<<", "<<i;
		if (j!=-1) cerr << ", " << j;
		cerr <<")]: invalid permutation" << endl << endl;
		return -1;
	}
	else if (j == -1) j = i+n_block; // Contiguous block permutation (by default, since j = -1 is a default parameter)
	unsigned int t;
	for (int n = 0; n < n_block; n++) {
		t = Population[p][i+n];
		Population[p][i+n] = Population[p][j+n];
		Population[p][j+n] = t;
	}
	return 0;
}

// Shift n_shift indices by n_places places (>< 0), starting from the i-th index
int Genetic :: Shift(unsigned int p, int i, int n_block, int n_shift) { 
	if (p >= Npop || p < 0 || i <= 0 || abs(i) >= Nx || abs(i+n_shift+n_block) > Nx || i+n_shift <= 0) {
		cerr << "Warning [Genetic :: Shift("<<p<<", "<<i<<", "<<n_block<<", "<<n_shift<<")]: invalid shift" << endl << endl;
		return -1;
	}
	else {
		unsigned int t;
		if (n_shift >= 0) {
			for (int k = n_block-1; k >= 0; k--) {
				for (int j = i; j < i+n_shift; j++) {
					t = Population[p][j+k];
					Population[p][j+k] = Population[p][j+k+1];
					Population[p][j+k+1] = t;
				}
			}
		}
		else {
			for (int k = 0; k <= n_block-1; k++) {
				for (int j = i; j > i+n_shift; j--) {
					t = Population[p][j+k];
					Population[p][j+k] = Population[p][j+k-1];
					Population[p][j+k-1] = t;
				}
			}
		}
	}
	return 0;
}

int Genetic :: Invert(unsigned int p, unsigned int i, unsigned int n_inv) {
	if (p >= Npop) {
		cerr << "Error [Genetic :: Invert(unsigned int, unsigned int, unsigned int)]: population index out of range" << endl;
		exit(-1);
	}
	if (i == 0 || i+n_inv > Nx) {
		cerr << "Warning [Genetic :: Invert("<<p<<", "<<i<<", "<<n_inv<<")]: invalid inversion" << endl << endl;
		return -1;
	}
	else {
		for (unsigned int k = 0; k < (unsigned int)(n_inv/2); k++) Swap(p, i+k, i+n_inv-1-k);	
	}
	return 0;
}

void Genetic :: CrossingOver(unsigned int p, unsigned int q, unsigned int cut) { // cut = first index of the mutable section
	if (p >= Npop || q >= Npop) {
		cerr << "Error [Genetic :: CrossingOver(unsigned int, unsigned int, unsigned int)]: population index out of range" << endl;
		exit(-1);
	}
	if (cut == 0 || cut >= Nx) {
		cerr << "Warning [Genetic :: CrossingOver("<<p<<", "<<q<<", "<<cut<<")]: invalid crossing over" << endl << endl;
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



// ***** // ***** // ***** // ***** // *****   AUXILIARY ALGORITHMS   ***** // ***** // ***** // ***** // ***** //

// Shift n_shift indices by n_places places, starting from the i-th index
void Genetic :: Shift(vector<unsigned int>& V, int i, int i_new) const { 
	if (abs(i>i_new?i:i_new) > V.size() || (i<i_new?i:i_new) < 0) {
		cerr << "Warning [Genetic :: Shift(vector<unsigned int>, "<<i<<", "<<i_new<<")]: invalid shift" << endl << endl;
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

int Genetic :: BinarySearch(vector<unsigned int> V, unsigned int first, unsigned int last, unsigned int find) const {
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

// ***** // ***** // ***** // ***** // *****   SCREEN OUTPUT METHODS   ***** // ***** // ***** // ***** // ***** //

void Genetic :: PrintX() const {
	cout << "*** X array ***" << endl;
	for (unsigned int n = 0; n < Nx; n++) {
		cout << "[" << n << "] -> ";
		X[n].Print();
	}
	cout << endl;
}

void Genetic :: PrintPopulation() const {
	cout << "*** Population array ***" << endl;
	for (unsigned int n = 0; n < Npop; n++) PrintIndividual(n);
	cout << endl;
}

void Genetic :: PrintIndividual(unsigned int i) const {
	cout << "[" << i << "] = " << "(";
	for (unsigned int n = 0; n < Nx-1; n++) cout << Population[i][n] << ", ";
	cout << Population[i][Nx-1] << ")" << endl;
}

void Genetic :: PrintL() const {
	cout << "*** L array ***" << endl;
	cout << "L = (";
	for (unsigned int n = 0; n < Npop-1; n++) cout << L[n] << ", ";
	cout << L[Npop-1] << ")" << endl;
	cout << endl;
}

void Genetic :: PrintOrder() const {
	cout << "*** Order array ***" << endl;
	cout << "Order = (";
	for (unsigned int n = 0; n < Npop-1; n++) cout << Order[n] << ", ";
	cout << Order[Npop-1] << ")" << endl;
	cout << endl;
}

void Genetic :: PrintV(vector<unsigned int> V) const {
	cout << "Vector:  ";
	for (unsigned int i = 0; i < V.size(); i++) {
		cout << V[i] << " ";
	}
	cout << endl;
}

void Genetic :: PrintBest() const {
	cout << "> Best Individual: "; PrintIndividual(Order[0]); // Order[0] -> the individuals are not ordered by Sort()
	cout << "  L = " << L[0] << endl << endl; // [0] -> the array L |>IS<| ordered by Sort() 
}

void Genetic :: PrintBestL() const {
	cout << "> Best L = " << L[0] << endl << endl; // [0] -> the array L |>IS<| ordered by Sort() 
}


void Genetic :: LoadPoints(const string filename) {
	ifstream in(filename);
	if (in.is_open()) {
		string buff;
		unsigned int _d, _Nx;
		in >> buff >> _d >> buff >> _Nx;
		SetParameters(Npop, _Nx, _d);
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
		cerr << "Error [Genetic :: LoadPoints(const string)]: cannot open input file" << endl;
		exit(-1);
	}
	return;
}

void Genetic :: OutputPoints(const string filename) const {
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
	}
	else {
		cerr << "Error [Genetic :: OutputPoints(const string)]: cannot open input file" << endl;
		exit(-1);
	}
	return;
}

void Genetic :: OutputBestPath(const string filename) const { // Output the best path
	if (FinalPath.size() == 0) {
		cerr << "Error [Genetic :: OutputFinalPath(const string)]: final path not set" << endl;
		exit(-1);
	}
	cout << "> Saving to file <"+filename+">..."<< endl;
	ofstream out(filename);
	if (out.is_open()) {
		out << "#SpaceDimension: " << d << endl << "#NumberOfPoints: " << Nx << endl;
		out << "#L: " << L[0] << endl;
		for (unsigned int i = 0; i < Nx-1; i++) {
			out << i+1 << "   ";
			for (unsigned int k = 0; k < d-1; k++) out << X[FinalPath[i]][k] << " ";
			out << X[FinalPath[i]][d-1] << endl;
		}	
		out << Nx << "   ";
		for (unsigned int k = 0; k < d-1; k++) out << X[FinalPath[Nx-1]][k] << " ";
		out << X[FinalPath[Nx-1]][d-1];
	}
	else {
		cerr << "Error [Genetic :: OutputPath(const string)]: cannot open input file" << endl;
		exit(-1);
	}
	return;
}


