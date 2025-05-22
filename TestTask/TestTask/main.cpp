#include <iostream>
#include <complex>
#include <cmath>
#include <ctime>
#include <random>
#include <fstream>

using namespace std;


double randn()
{
	// 1. Инициализация устройства для получения случайного зерна (seed)
		static random_device rd;	// Источник энтропии (зависит от системы)

	// 2. Инициализация генератора случайных чисел Mersenne Twister
		static mt19937 gen(rd());	// Генератор с seed от random_device

	// 3. Создание стандартного нормального распределения (N(0,1))
		static normal_distribution<> dist(0, 1);

	// 4. Генерация и возврат случайного числа
		return dist(gen);

	//// 1. Инициализация генератора случайных чисел
	//	default_random_engine generator;

	//// 2. Создание стандартного нормального распределения (N(0,1))
	//	normal_distribution<double> distribution(0, 1);

	//// 3. Генерация и возврат случайного числа
	//	return distribution(generator);
}


void de2bi(int* Output, int* M)
{
	int k = log2(*M);									// Кол-во бит в символе
	int* vsp = new int[*M];								// Массив 0 : M-1
	for (int i = 0; i < *M; i++) { *(vsp + i) = i; }	// Массив 0 : M-1
	
	for (int i = 0; i < *M; i++) {						// |
		for (int j = k-1; j >= 0; j--) {				// |
			*(Output + k*i+j) = *(vsp+i) % 2;			// | Основной цикл функции
			*(vsp + i) = *(vsp + i) / 2;				// |
		}												// |
	}													// |

	delete[] vsp;										// Освобождение памяти
}


class Mapper
{
	public:
		int M;

		complex<double>* constellation;


		void CalcAuxParams()
		{
			// Всякие проверки
				
			// Расчёт сигнального созвездия
				constellation = new complex<double>[M];
				int k = log2(M);							// Кол-во в символе
				int sizeCB = M*k;							// Кол-во бит, необходимое для получения сигнального созвездия
				int* ConstBits = new int[sizeCB];			// Массив бит, для создания созвездия
			
				de2bi(ConstBits, &M);						// Заполнение ConstBits
				StepTx(constellation, ConstBits, &M);		// Получение сигнального созвездия

				delete[] ConstBits;							// Освобождение памяти
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
			// randn() генерирует нормально-распределённую СВ с нулевым средним и дисперсией 1, СКО 
			// следовательно тоже единица. В нашем случае комплексный шум с дисперсией 2*sigma^2,
			// т.е. для каждой квадратуры шум с дисперсией sigma^2. Вспомнив из теории вероятности 
			// следующее тождество: D[a*ksi] = a^2 * D[ksi], становится понятно, что randn нужно
			// домножить на sigma, для получения шума с требуемой дисперсией sigma^2

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

	// Инициализация параметров моделирования
		Mapper QAM;												// Инициализация объекта модулятора
		Channel AWGN;											// Инициализация объекта канала
		QAM.M = 64;												// Порядок модуляции
		int Ns = 1e5;											// Кол-во модуляционных символов
		int N = Ns*log2(QAM.M);									// Кол-во информационных бит

	// Инициализация параметров кривой помехоустойчивости
		int Nav = 1e2;											// Кол-во усреднений
		double VARmin = 0.05;									// Начальная дисперсия
		double VARstep = 0.025;									// Шаг дисперсии
		double VARmax = 0.5;									// Конечная дисперсия
		int VARsize = floor((VARmax - VARmin) / VARstep) + 1;	// Длина массива дисперсий

	// Запускаем расчёт вспомогательных параметров и проверку параметров
		QAM.CalcAuxParams();

	// Инициализация массивов моделирования
		int* Tx_InfBit = new int[N];							// Массив информационных бит
		int* Rx_InfBit = new int[N];							// Массив принятых инф. бит
		complex<double>* Tx_ModSym = new complex<double>[Ns];	// Массив модуляционных символов
		complex<double>* Rx_ModSym = new complex<double>[Ns];	// Массив принятых мод. символов

	// Инициализация массивов кривой помехоустойчивости
		double* VAR = new double[VARsize];						// Массив дисперсий
		double* BER = new double[VARsize];						// Массив битовой вероятности ошибки

		for (int i = 0; i < VARsize; i++) { *(VAR + i) = VARmin + VARstep*i; }	// Заполнение
		for (int i = 0; i < VARsize; i++) { *(BER + i) = 0; }					// Заполнение

	// Основной цикл моделирования
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

	// Нормировка и вывод BER
		cout << endl;
		for (int h = 0; h < VARsize; h++) {
			*(BER + h) = *(BER + h) / N / Nav;
			cout << *(BER + h) << "\t";
		} cout << endl;
		
	// Запись в файл
		ofstream fout;
		fout.open("BER_from_Cpp.txt");

		if (!fout.is_open()) {
			cout << "Ошибка открытия файла" << endl;
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

	// Освобождение памяти
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