#include <iostream>
#include <complex>
#include <cmath>
#include <ctime>
#include <random>
#include <fstream>

using namespace std;


double randn()
{
	// 1. Initializing the device to receive random seed
		static random_device rd;

	// 2. Initializing the random number generator Mersenne Twister
		static mt19937 gen(rd());

	// 3. Creating a standard normal distribution (N(0,1))
		static normal_distribution<> dist(0, 1);

	// 4. Generating and returning a random number
		return dist(gen);
}


void de2bi(int* Output, int* M)
{
	int k = log2(*M);									// Number of bits in a symbol
	int* vsp = new int[*M];								// Array  [0 : M-1]
	for (int i = 0; i < *M; i++) { *(vsp + i) = i; }	// Array  [0 : M-1]
	
	for (int i = 0; i < *M; i++) {						// |
		for (int j = k-1; j >= 0; j--) {				// |
			*(Output + k*i+j) = *(vsp+i) % 2;			// | The main loop of function
			*(vsp + i) = *(vsp + i) / 2;				// |
		}												// |
	}													// |

	delete[] vsp;										// Clear memory
}


class Mapper
{
	public:
		int M;								// The order of modulation

		complex<double>* constellation;		// The signal constellation


		void CalcAuxParams()
		{
			// Some checks
				// Maybe in another version
				
			// Calculation of the signal constellation
				constellation = new complex<double>[M];
				int k = log2(M);							// Number of symbols
				int sizeCB = M*k;							// Number of bits required to getting the constellation
				int* ConstBits = new int[sizeCB];			// Array of bits to create a constellation
			
				de2bi(ConstBits, &M);						// Filling ConstBits
				StepTx(constellation, ConstBits, &M);		// Getting the constellation

				delete[] ConstBits;							// Clear memory
		}

		void StepTx(complex<double>* Output, int* Input, int* Ns)
		{
			if (M == 4) {
				for (int i = 0; i < *Ns; i++) {
					*(Output + i) = complex<double>(
						(2*(*(Input+2*i+0) == 1)-1) ,
						(2*(*(Input+2*i+1) == 0)-1) );
				}
			}
			else if (M == 16) {
				for (int i = 0; i < *Ns; i++) {
					*(Output + i) = complex<double>(
						(2*(*(Input+4*i+0) == 1)-1)  *  (2 + 2*(*(Input+4*i+1) == 0)-1) ,
						(2*(*(Input+4*i+2) == 0)-1)  *  (2 + 2*(*(Input+4*i+3) == 0)-1) );
				}
			}
			else if (M == 64) {
				for (int i = 0; i < *Ns; i++) {
					*(Output + i) = complex<double>(
						 ( 2*(*(Input+6*i+0)==1)-1 )  *  (4 + (2*(*(Input+6*i+1)==0)-1)*(2+(2*(*(Input+6*i+2)==0)-1))), 
                         ( 2*(*(Input+6*i+3)==0)-1 )  *  (4 + (2*(*(Input+6*i+4)==0)-1)*(2+(2*(*(Input+6*i+5)==0)-1))) );
				}
			}
		}

		void StepRx(int* Output, complex<double>* Input, int* Ns)
		{
			if (M == 4) {
				for (int i = 0; i < *Ns; i++) {
					*(Output + 2*i+0) = real(*(Input+i)) > 0;
					*(Output + 2*i+1) = imag(*(Input+i)) < 0;
				}
			}
			else if (M == 16) {
				for (int i = 0; i < *Ns; i++) {
					*(Output + 4*i+0) = real(*(Input+i)) > 0;
					*(Output + 4*i+1) = abs(real(*(Input+i))) < 2;
					*(Output + 4*i+2) = imag(*(Input+i)) < 0;
					*(Output + 4*i+3) = abs(imag(*(Input+i))) < 2;
				}
			}
			else if (M == 64) {
				for (int i = 0; i < *Ns; i++) {
					*(Output + 6*i+0) = real(*(Input+i)) > 0;
					*(Output + 6*i+1) = abs(real(*(Input+i))) < 4;
					*(Output + 6*i+2) = abs(abs(real(*(Input+i)))-4) < 2;
					*(Output + 6*i+3) = imag(*(Input+i)) < 0;
					*(Output + 6*i+4) = abs(imag(*(Input+i))) < 4;
					*(Output + 6*i+5) = abs(abs(imag(*(Input+i)))-4) < 2;
				}
			}
		}

		void ClearMemory()
		{
			delete[] constellation;
		}
	
};


class Channel
{
	public:
		void Step(complex<double>* Output, complex<double>* Input, double* sigma2, int* Ns)
		{
			// randn() generates a normally distributed RV with zero mean and variance of 1, MSE 
			// is also equal to 1. In our case, complex noise with variance 2*sigma^2, that is, 
			// for each quadrature noise with variance sigma^2. Remembering from probability theory
			// the following identity: D[a*ksi] = a^2 * D[ksi], it becomes clear that randn needs
			// multiply by sigma, to get the noise with the required variance sigma^2

			double sigma = sqrt(*sigma2);

			for (int i = 0; i < *Ns; i++) {
				*(Output+i) = *(Input+i) + complex<double>(sigma*randn(), sigma*randn());
			}
		}
};


int main()
{
	srand(time(NULL));

	// Initialization of modeling parameters
		Mapper QAM;												// Initialization of the modulator object
		Channel AWGN;											// Initialization of the channel object
		QAM.M = 64;												// The order of modulation
		int Ns = 1e4;											// Number of modulation symbols
		int N = Ns*log2(QAM.M);									// Number of information bits

	// Initialization of noise immunity curve parameters
		int Nav = 1e2;											// Number of averages
		double VARmin = 0.05;									// Initial variance
		double VARstep = 0.025;									// The variance step
		double VARmax = 0.5;									// The final variance
		int VARsize = floor((VARmax - VARmin) / VARstep) + 1;	// The length of the variance array

	// Starting the calculation of auxiliary parameters and checking the parameters
		QAM.CalcAuxParams();

	// Initializing modeling arrays
		int* Tx_InfBit = new int[N];							// Array of information bits
		int* Rx_InfBit = new int[N];							// Array of received inf. bits
		complex<double>* Tx_ModSym = new complex<double>[Ns];	// Array of modulation symbols
		complex<double>* Rx_ModSym = new complex<double>[Ns];	// Array of received mod. symbols

	// Initializing noise immunity curve arrays
		double* VAR = new double[VARsize];						// Array of variances
		double* BER = new double[VARsize];						// Array of bit error probability

		for (int i = 0; i < VARsize; i++) { *(VAR + i) = VARmin + VARstep*i; }	// Initial filling VAR
		for (int i = 0; i < VARsize; i++) { *(BER + i) = 0; }					// Initial filling BER

	// Main loop
		for (int i = 0; i < VARsize; i++) {
			double Sigma2 = *(VAR+i);
			cout << Sigma2 << "\t";
			for (int j = 0; j < Nav; j++) {
				// Source
					for (int h = 0; h < N; h++) { 
						*(Tx_InfBit + h) = rand() % 2; 
					}

				// Mapper
					QAM.StepTx(Tx_ModSym, Tx_InfBit, &Ns);

				// Channel
					AWGN.Step(Rx_ModSym, Tx_ModSym, &Sigma2, &Ns);

				// Demapper
					QAM.StepRx(Rx_InfBit, Rx_ModSym, &Ns);

				// BER
					for (int h = 0; h < N; h++) {
						*(BER+i) = *(BER+i) + (*(Rx_InfBit+h) != *(Tx_InfBit+h));
					}
			}
		}

	// Normalize and display BER
		cout << endl;
		for (int h = 0; h < VARsize; h++) {
			*(BER + h) = *(BER + h) / N / Nav;
			cout << *(BER + h) << "\t";
		} cout << endl;
		
	// Writing in file
		ofstream fout;
		fout.open("BER_from_Cpp.txt");

		if (!fout.is_open()) {
			cout << "File opening error" << endl;
		}
		else {
			for (int h = 0; h < VARsize; h++) {
				fout << *(VAR + h) << "\t";
			} fout << endl;
			for (int h = 0; h < VARsize; h++) {
				fout << *(BER + h) << "\t";
			} fout << endl;
			for (int h = 0; h < QAM.M; h++) {
				fout << *(QAM.constellation + h) << "\t";
			} fout << endl;
		}

		fout.close();

	// Clear memory
		delete[] Tx_InfBit;
		delete[] Rx_InfBit;
		delete[] Tx_ModSym;
		delete[] Rx_ModSym;
		delete[] VAR;
		delete[] BER;
		QAM.ClearMemory();

	system("pause");
	return 0;
}