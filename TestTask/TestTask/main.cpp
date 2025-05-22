#include <iostream>
#include <complex>
#include <cmath>
#include <ctime>
#include <random>
#include <fstream>

using namespace std;


double randn()
{
	// 1. ������������� ���������� ��� ��������� ���������� ����� (seed)
		static random_device rd;	// �������� �������� (������� �� �������)

	// 2. ������������� ���������� ��������� ����� Mersenne Twister
		static mt19937 gen(rd());	// ��������� � seed �� random_device

	// 3. �������� ������������ ����������� ������������� (N(0,1))
		static normal_distribution<> dist(0, 1);

	// 4. ��������� � ������� ���������� �����
		return dist(gen);

	//// 1. ������������� ���������� ��������� �����
	//	default_random_engine generator;

	//// 2. �������� ������������ ����������� ������������� (N(0,1))
	//	normal_distribution<double> distribution(0, 1);

	//// 3. ��������� � ������� ���������� �����
	//	return distribution(generator);
}


void de2bi(int* Output, int* M)
{
	int k = log2(*M);									// ���-�� ��� � �������
	int* vsp = new int[*M];								// ������ 0 : M-1
	for (int i = 0; i < *M; i++) { *(vsp + i) = i; }	// ������ 0 : M-1
	
	for (int i = 0; i < *M; i++) {						// |
		for (int j = k-1; j >= 0; j--) {				// |
			*(Output + k*i+j) = *(vsp+i) % 2;			// | �������� ���� �������
			*(vsp + i) = *(vsp + i) / 2;				// |
		}												// |
	}													// |

	delete[] vsp;										// ������������ ������
}


class Mapper
{
	public:
		int M;

		complex<double>* constellation;


		void CalcAuxParams()
		{
			// ������ ��������
				
			// ������ ����������� ���������
				constellation = new complex<double>[M];
				int k = log2(M);							// ���-�� � �������
				int sizeCB = M*k;							// ���-�� ���, ����������� ��� ��������� ����������� ���������
				int* ConstBits = new int[sizeCB];			// ������ ���, ��� �������� ���������
			
				de2bi(ConstBits, &M);						// ���������� ConstBits
				StepTx(constellation, ConstBits, &M);		// ��������� ����������� ���������

				delete[] ConstBits;							// ������������ ������
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
			// randn() ���������� ���������-������������� �� � ������� ������� � ���������� 1, ��� 
			// ������������� ���� �������. � ����� ������ ����������� ��� � ���������� 2*sigma^2,
			// �.�. ��� ������ ���������� ��� � ���������� sigma^2. �������� �� ������ ����������� 
			// ��������� ���������: D[a*ksi] = a^2 * D[ksi], ���������� �������, ��� randn �����
			// ��������� �� sigma, ��� ��������� ���� � ��������� ���������� sigma^2

			double sigma = sqrt(*sigma2);

			for (int i = 0; i < *Ns; i++) {
				double a = randn();
				double b = randn();
				*(Output+i) = *(Input+i) + complex<double>(sigma*a, sigma*b);
			}
		}
};


int main()
{
	setlocale(LC_ALL, "rus");
	srand(time(NULL));

	// ������������� ���������� �������������
		Mapper QAM;												// ������������� ������� ����������
		Channel AWGN;											// ������������� ������� ������
		QAM.M = 64;												// ������� ���������
		int Ns = 1e5;											// ���-�� ������������� ��������
		int N = Ns*log2(QAM.M);									// ���-�� �������������� ���

	// ������������� ���������� ������ ������������������
		int Nav = 1e2;											// ���-�� ����������
		double VARmin = 0.05;									// ��������� ���������
		double VARstep = 0.025;									// ��� ���������
		double VARmax = 0.5;									// �������� ���������
		int VARsize = floor((VARmax - VARmin) / VARstep) + 1;	// ����� ������� ���������

	// ��������� ������ ��������������� ���������� � �������� ����������
		QAM.CalcAuxParams();

	// ������������� �������� �������������
		int* Tx_InfBit = new int[N];							// ������ �������������� ���
		int* Rx_InfBit = new int[N];							// ������ �������� ���. ���
		complex<double>* Tx_ModSym = new complex<double>[Ns];	// ������ ������������� ��������
		complex<double>* Rx_ModSym = new complex<double>[Ns];	// ������ �������� ���. ��������

	// ������������� �������� ������ ������������������
		double* VAR = new double[VARsize];						// ������ ���������
		double* BER = new double[VARsize];						// ������ ������� ����������� ������

		for (int i = 0; i < VARsize; i++) { *(VAR + i) = VARmin + VARstep*i; }	// ����������
		for (int i = 0; i < VARsize; i++) { *(BER + i) = 0; }					// ����������

	// �������� ���� �������������
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

	// ���������� � ����� BER
		cout << endl;
		for (int h = 0; h < VARsize; h++) {
			*(BER + h) = *(BER + h) / N / Nav;
			cout << *(BER + h) << "\t";
		} cout << endl;
		
	// ������ � ����
		ofstream fout;
		fout.open("BER_from_Cpp.txt");

		if (!fout.is_open()) {
			cout << "������ �������� �����" << endl;
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

	// ������������ ������
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