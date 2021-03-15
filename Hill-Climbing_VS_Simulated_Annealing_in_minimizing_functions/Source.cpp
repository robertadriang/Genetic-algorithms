#include <math.h>
#include <iostream>
#include <vector>
#include <time.h>
#include <thread>
#include <sstream>
#include <fstream>
using namespace std;
#define PI 3.14159
double a, b;
int prec = 5,l;
unsigned n;

int getId(thread::id id)
{
	stringstream buffer;
	buffer << id;
	return stoull(buffer.str());
}

string generateFileName(int f, int best = 0)
{
	int id = getId(this_thread::get_id());
	srand(id * clock());
	string name;
	if (best == 0)
		name = "HC_BEST_";
	else if (best == 1)
		name = "HC_FIRST_";
	else
		name = "SA_";
	switch (f)
	{
	case 1:
		name.append("D_");
		break;
	case 2:
		name.append("S_");
		break;
	case 3:
		name.append("R_");
		break;
	case 4:
		name.append("M_");
		break;
	}
	name.append(to_string(n));
	name.append("_");
	name.append(to_string(id).c_str());
	name.append(".txt");
	return name;
}

double deJong(vector<double>& values)
{
	double result = 0;
	int size = values.size();
	for (int i = 0; i < size; ++i)
		result += values[i] * values[i];
	return result;
}

double Schwefel(vector<double>& values)
{
	double result = 0;
	int size = values.size();
	for (int i = 0; i < size; ++i)
	{
		result += (-values[i] * sin(sqrt(abs(values[i]))));
	}
	return result;
}

double Rastrigin(vector<double>& values)
{
	double result = 10 * values.size();
	int size = values.size();
	for (int i = 0; i < size; ++i)
	{
		result += (values[i] * values[i] - 10 * cos(2 * PI * values[i]));
	}
	return result;
}

double Michalewicz(vector<double>& values)
{
	double result = 0;
	int size = values.size();
	for (int i = 1; i <= size; ++i)
	{
		result += (sin(values[i-1]) * pow(sin(i * values[i-1] * values[i-1] / PI), 20));
	}
	return -result;
}

int computeComponentLength(double a, double b, int prec)
{
	return ceil(log((b - a) * pow(10, prec)) / log(2));
}

double decodeDimension(vector<char>::iterator itStart, vector<char>::iterator itEnd, int l, int a, int b)
{
	unsigned long xInt = 0;
	double xScaled;
	for (auto it = itStart; it != itEnd; ++it)
	{
		xInt *= 2;
		xInt += *it;
	}
	xScaled = xInt / (pow(2, l) - 1);
	return (xScaled * (b - a) + a);
}

vector<double> decode(vector<char>& bits, int l, unsigned n, double a, double b)
{
	vector<double> values;
	vector<char>::iterator itStart, itEnd;
	for (int i = 0; i < n; ++i)
	{
		itStart = bits.begin() + i * l;
		itEnd = itStart + l;
		double x = decodeDimension(itStart, itEnd, l, a, b);
		values.push_back(x);
	}
	return values;
}

double minFinderBest(vector<char>& bits, unsigned n, int l, double a, double b, double (*func)(vector<double>&))
{
	auto aux = decode(bits, l, n, a, b);
	double min = func(aux);
	int size = bits.size();
	bool local = false;
	int minInd = 0;
	do
	{
		local = false;
		double initial = min;
		for (int i = 0; i < size; ++i)
		{
			bits[i] = !bits[i];
			aux = decode(bits, l, n, a, b);
			double val = func(aux);
			if (val < min)
			{
				min = val;
				minInd = i;
			}
			bits[i] = !bits[i];
		}
		if (min == initial)
			local = true;
		else
			bits[minInd] = !bits[minInd];
	} while (!local);
	return min;
}

void hillClimbingBest(unsigned n, double a,double b,double(*func)(vector < double>&), int functionNumber)
{
	clock_t start, end;
	start = clock();
	int id = getId(this_thread::get_id());
	string name = generateFileName(functionNumber);
	ofstream myfile(name);
	srand(id * clock());

	double globalMin = DBL_MAX,localMin;
	int itNrMin = -1;
	int L = l * n;
	vector<char> bits;
	vector<double> minBits;
	for (int t = 0; t < 10000; ++t)
	{
		bits.clear();
		for (int i = 0; i < L; ++i)
			bits.push_back(rand() % 2);
		localMin = minFinderBest(bits,n,l,a,b,func);
		if (localMin < globalMin)
		{
			globalMin = localMin;
			itNrMin = t;
			minBits = decode(bits, l, n, a, b);
		}
	}
	myfile << "GlobalMinimum=" << globalMin << " found at " << itNrMin << " iteration at coordinates: \n";
	int size = minBits.size();
	for (int i = 0; i < size; ++i)
		myfile << minBits[i] << ' ';
	myfile << '\n';
	end = clock();
	myfile << "Time taken: " << (double)((double)(end - start) / (double)CLOCKS_PER_SEC) << " seconds.";
	myfile.close();
}

double minFinderFirst(vector<char>& bits, unsigned n, int l, double a, double b, double (*func)(vector<double>&))
{

	auto aux= decode(bits, l, n, a, b);
	double min = func(aux);
	int size = bits.size();
	bool local = false;
	do
	{
		local = false;
		double initial = min;
		for (int i = 0; i < size; ++i)
		{
			bits[i] = !bits[i];
			aux = decode(bits, l, n, a, b);
			double val = func(aux);
			if (val < min)
			{
				min = val;
				break;
			}
			bits[i] = !bits[i];
		}
		if (min == initial)
		{
			local = true;
		}
	} while (!local);
	return min;
}

void hillClimbingFirst(int n, double a, double b,double(*func)(vector<double>&), int functionNumber)
{
	clock_t start, end;
	start = clock();
	int id = getId(this_thread::get_id());
	string name = generateFileName(functionNumber, 1);
	ofstream myfile(name);
	srand(id * clock());

	double globalMin = DBL_MAX, localMin;
	int itNrMin = -1;
	int L = l * n;
	vector<char> bits;
	vector<double> minBits;
	for (int t = 0; t < 10000; ++t)
	{
		bits.clear();
		for (int i = 0; i < L; ++i)
			bits.push_back(rand() % 2);
		localMin = minFinderFirst(bits, n, l, a, b, func);
		if (localMin < globalMin)
		{
			globalMin = localMin;
			itNrMin = t;
			minBits = decode(bits, l, n, a, b);
		}
	}

	myfile << "GlobalMinimum=" << globalMin << " found at " << itNrMin << " iteration at coordinates: \n";
	int size = minBits.size();
	for (int i = 0; i < size; ++i)
		myfile << minBits[i] << ' ';
	myfile << '\n';
	end = clock();
	myfile << "Time taken: " << (double)((double)(end - start) / (double)CLOCKS_PER_SEC) << " seconds.";
	myfile.close();
}

double random01()
{
	return (double)(rand() % 100000) / (double)100000;
}

void simulatedAnnealing(unsigned n, double a,double b,double(*func)(vector<double>&), int functionNumber)
{
	clock_t start, end;
	start = clock();
	int id = getId(this_thread::get_id());
	string name = generateFileName(functionNumber, 2);
	ofstream myfile(name);
	srand(id * clock());

	int L = l * n;
	vector<char> bits;
	for (int i = 0; i < L; ++i)
		bits.push_back(rand() % 2);
	auto aux = decode(bits, l, n, a, b);
	double min = func(aux);
	vector<double>minBits;
	double temperature = 100;
	
	while (temperature > 0.0000001)
	{
		for (int i = 0; i < 10000; ++i)
		{
			int rPosition = rand() % L;
			bits[rPosition] = !bits[rPosition];
			aux = decode(bits, l, n, a, b);
			double local = func(aux);
			if (local < min)
			{
				min = local;
				minBits = decode(bits, l, n, a, b);
			}
			else if (random01() < exp(-1 * abs(local - min) / temperature))
			{
				min = local;
				minBits = decode(bits, l, n, a, b);
			}
			else
				bits[rPosition] = !bits[rPosition];
		}
		temperature *= .995;
	}
	
	myfile << "GlobalMinimum=" << min << " found at coordinates: \n";
	int size = minBits.size();
	for (int i = 0; i < size; ++i)
		myfile << minBits[i] << ' ';
	myfile << '\n';
	end = clock();
	myfile << "Time taken: " << (double)((double)(end - start) / (double)CLOCKS_PER_SEC) << " seconds.";
	myfile.close();
}

int main()
{
	/*cout << "Starting Hill Climbing Best FUNCTION\n";
	l = computeComponentLength(a, b, prec);
	cout << "Preparing DeJong function \n";
	{
		cout << "Preparing threads for DeJong size 5 SA\n";
		n = 5;
		a = -5.12;
		b = 5.12;
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingBest, n, a, b, deJong, 1);
			thread t2(hillClimbingBest, n, a, b, deJong, 1);
			thread t3(hillClimbingBest, n, a, b, deJong, 1);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[DeJong] Size 5 completed! \n";

		n = 10;
		cout << "Preparing threads for DeJong size 10 best improvement\n";
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingBest, n, a, b, deJong, 1);
			thread t2(hillClimbingBest, n, a, b, deJong, 1);
			thread t3(hillClimbingBest, n, a, b, deJong, 1);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[DeJong] Size 10 completed! \n";

		n = 30;
		cout << "Preparing threads for DeJong size 30 best improvement\n";
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingBest, n, a, b, deJong, 1);
			thread t2(hillClimbingBest, n, a, b, deJong, 1);
			thread t3(hillClimbingBest, n, a, b, deJong, 1);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[DeJong] Size 30 completed! \n";
		cout << "[Best improvement] DeJong function finished! \n";
	}
	//Here end DeJong SA

	cout << "Preparing Schwefel function! \n";
	{
		cout << "Preparing threads for Schwefel size 5 best\n";
		n = 5;
		a = -500;
		b = 500;
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingBest, n, a, b, Schwefel, 2);
			thread t2(hillClimbingBest, n, a, b, Schwefel, 2);
			thread t3(hillClimbingBest, n, a, b, Schwefel, 2);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Schwefel] Size 5 completed! \n";

		n = 10;
		cout << "Preparing threads for Schwefel size 10 best improvement\n";
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingBest, n, a, b, Schwefel, 2);
			thread t2(hillClimbingBest, n, a, b, Schwefel, 2);
			thread t3(hillClimbingBest, n, a, b, Schwefel, 2);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Schwefel] Size 10 completed! \n";

		n = 30;
		cout << "Preparing threads for Schwefel size 30 best improvement\n";
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingBest, n, a, b, Schwefel, 2);
			thread t2(hillClimbingBest, n, a, b, Schwefel, 2);
			thread t3(hillClimbingBest, n, a, b, Schwefel, 2);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Schwefel] Size 30 completed! \n";
		cout << "[Best improvement] Schwefel function finished! \n";
	}
	//Here ends Schwefel

	cout << "Preparing Rastrigin function!\n";
	{
		cout << "Preparing threads for Rastrigin size 5 SA\n";
		n = 5;
		a = -5.12;
		b = 5.12;
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingBest, n, a, b, Rastrigin, 3);
			thread t2(hillClimbingBest, n, a, b, Rastrigin, 3);
			thread t3(hillClimbingBest, n, a, b, Rastrigin, 3);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Rastrigin] Size 5 completed! \n";

		n = 10;
		cout << "Preparing threads for Rastrigin size 10 best improvement\n";
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingBest, n, a, b, Rastrigin, 3);
			thread t2(hillClimbingBest, n, a, b, Rastrigin, 3);
			thread t3(hillClimbingBest, n, a, b, Rastrigin, 3);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Rastrigin] Size 10 completed! \n";

		n = 30;
		cout << "Preparing threads for Rastrigin size 30 best improvement\n";
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingBest, n, a, b, Rastrigin, 3);
			thread t2(hillClimbingBest, n, a, b, Rastrigin, 3);
			thread t3(hillClimbingBest, n, a, b, Rastrigin, 3);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Rastrigin] Size 30 completed! \n";
		cout << "[SIMULATED ANNEALING] Rastrigin function finished! \n";
	}

	cout << "Preparing Michalewicz's function!\n";
	{
		cout << "Preparing threads for Michalewicz size 5 best\n";
		n = 5;
		a = 0;
		b = PI;
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingBest, n, a, b, Michalewicz, 4);
			thread t2(hillClimbingBest, n, a, b, Michalewicz, 4);
			thread t3(hillClimbingBest, n, a, b, Michalewicz, 4);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Michalewicz] Size 5 completed! \n";

		n = 10;
		cout << "Preparing threads for Michalewicz size 10 best improvement\n";
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingBest, n, a, b, Michalewicz, 4);
			thread t2(hillClimbingBest, n, a, b, Michalewicz, 4);
			thread t3(hillClimbingBest, n, a, b, Michalewicz, 4);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Michalewicz] Size 10 completed! \n";

		n = 30;
		cout << "Preparing threads for Michalewicz size 30 Best improvement\n";
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingBest, n, a, b, Michalewicz, 4);
			thread t2(hillClimbingBest, n, a, b, Michalewicz, 4);
			thread t3(hillClimbingBest, n, a, b, Michalewicz, 4);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Michalewicz] Size 30 completed! \n";
		cout << "[Best Improvement] Michalewicz function finished! \n";
	}
	cout << "Best Improvement FINISHED!";
	*/

	/*	cout << "Starting Hill Climbing First FUNCTION\n";
	l = computeComponentLength(a, b, prec);
	cout << "Preparing DeJong function \n";
	{
		cout << "Preparing threads for DeJong size 5 SA\n";
		n = 5;
		a = -5.12;
		b = 5.12;
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingFirst, n, a, b, deJong, 1);
			thread t2(hillClimbingFirst, n, a, b, deJong, 1);
			thread t3(hillClimbingFirst, n, a, b, deJong, 1);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[DeJong] Size 5 completed! \n";

		n = 10;
		cout << "Preparing threads for DeJong size 10 first improvement\n";
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingFirst, n, a, b, deJong, 1);
			thread t2(hillClimbingFirst, n, a, b, deJong, 1);
			thread t3(hillClimbingFirst, n, a, b, deJong, 1);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[DeJong] Size 10 completed! \n";

		n = 30;
		cout << "Preparing threads for DeJong size 30 first improvement\n";
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingFirst, n, a, b, deJong, 1);
			thread t2(hillClimbingFirst, n, a, b, deJong, 1);
			thread t3(hillClimbingFirst, n, a, b, deJong, 1);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[DeJong] Size 30 completed! \n";
		cout << "[First improvement] DeJong function finished! \n";
	}
	//Here end DeJong SA

	cout << "Preparing Schwefel function! \n";
	{
		cout << "Preparing threads for Schwefel size 5 SA\n";
		n = 5;
		a = -500;
		b = 500;
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingFirst, n, a, b, Schwefel, 2);
			thread t2(hillClimbingFirst, n, a, b, Schwefel, 2);
			thread t3(hillClimbingFirst, n, a, b, Schwefel, 2);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Schwefel] Size 5 completed! \n";

		n = 10;
		cout << "Preparing threads for Schwefel size 10 first improvement\n";
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingFirst, n, a, b, Schwefel, 2);
			thread t2(hillClimbingFirst, n, a, b, Schwefel, 2);
			thread t3(hillClimbingFirst, n, a, b, Schwefel, 2);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Schwefel] Size 10 completed! \n";

		n = 30;
		cout << "Preparing threads for Schwefel size 30 first improvement\n";
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingFirst, n, a, b, Schwefel, 2);
			thread t2(hillClimbingFirst, n, a, b, Schwefel, 2);
			thread t3(hillClimbingFirst, n, a, b, Schwefel, 2);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Schwefel] Size 30 completed! \n";
		cout << "[First improvement] Schwefel function finished! \n";
	}
	//Here ends Schwefel

	cout << "Preparing Rastrigin function!\n";
	{
		cout << "Preparing threads for Rastrigin size 5 SA\n";
		n = 5;
		a = -5.12;
		b = 5.12;
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingFirst, n, a, b, Rastrigin, 3);
			thread t2(hillClimbingFirst, n, a, b, Rastrigin, 3);
			thread t3(hillClimbingFirst, n, a, b, Rastrigin, 3);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Rastrigin] Size 5 completed! \n";

		n = 10;
		cout << "Preparing threads for Rastrigin size 10 first improvement\n";
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingFirst, n, a, b, Rastrigin, 3);
			thread t2(hillClimbingFirst, n, a, b, Rastrigin, 3);
			thread t3(hillClimbingFirst, n, a, b, Rastrigin, 3);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Rastrigin] Size 10 completed! \n";

		n = 30;
		cout << "Preparing threads for Rastrigin size 30 first improvement\n";
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingFirst, n, a, b, Rastrigin, 3);
			thread t2(hillClimbingFirst, n, a, b, Rastrigin, 3);
			thread t3(hillClimbingFirst, n, a, b, Rastrigin, 3);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Rastrigin] Size 30 completed! \n";
		cout << "[SIMULATED ANNEALING] Rastrigin function finished! \n";
	}

	cout << "Preparing Michalewicz's function!\n";
	{
		cout << "Preparing threads for Rastrigin size 5 SA\n";
		n = 5;
		a = 0;
		b = PI;
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingFirst, n, a, b, Michalewicz, 4);
			thread t2(hillClimbingFirst, n, a, b, Michalewicz, 4);
			thread t3(hillClimbingFirst, n, a, b, Michalewicz, 4);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Michalewicz] Size 5 completed! \n";

		n = 10;
		cout << "Preparing threads for Michalewicz size 10 first improvement\n";
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingFirst, n, a, b, Michalewicz, 4);
			thread t2(hillClimbingFirst, n, a, b, Michalewicz, 4);
			thread t3(hillClimbingFirst, n, a, b, Michalewicz, 4);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Michalewicz] Size 10 completed! \n";

		n = 30;
		cout << "Preparing threads for Michalewicz size 30 first improvement\n";
		l = computeComponentLength(a, b, prec);
		for (int i = 0; i < 10; ++i)
		{
			thread t1(hillClimbingFirst, n, a, b, Michalewicz, 4);
			thread t2(hillClimbingFirst, n, a, b, Michalewicz, 4);
			thread t3(hillClimbingFirst, n, a, b, Michalewicz, 4);
			t1.join();
			t2.join();
			t3.join();
			cout << "Iteration " << i << " ended\n";
		}
		cout << "[Michalewicz] Size 30 completed! \n";
		cout << "[First Improvement] Michalewicz function finished! \n";
	}
	cout << "First Improvement FINISHED!";*/

	/*cout << "Starting SIMULATED ANNEALING FUNCTION\n";
	l = computeComponentLength(a, b, prec);
	cout << "Preparing DeJong function \n";
	{
	cout << "Preparing threads for DeJong size 5 SA\n";
	n = 5;
	a = -5.12;
	b = 5.12;
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 10; ++i)
	{
		thread t1(simulatedAnnealing, n, a, b, deJong, 1);
		thread t2(simulatedAnnealing, n, a, b, deJong, 1);
		thread t3(simulatedAnnealing, n, a, b, deJong, 1);
		t1.join();
		t2.join();
		t3.join();
		cout << "Iteration " << i << " ended\n";
	}
	cout << "[DeJong] Size 5 completed! \n";

	n = 10;
	cout << "Preparing threads for DeJong size 10 first improvement\n";
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 10; ++i)
	{
		thread t1(simulatedAnnealing, n, a, b, deJong, 1);
		thread t2(simulatedAnnealing, n, a, b, deJong, 1);
		thread t3(simulatedAnnealing, n, a, b, deJong, 1);
		t1.join();
		t2.join();
		t3.join();
		cout << "Iteration " << i << " ended\n";
	}
	cout << "[DeJong] Size 10 completed! \n";

	n = 30;
	cout << "Preparing threads for DeJong size 30 first improvement\n";
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 10; ++i)
	{
		thread t1(simulatedAnnealing, n, a, b, deJong, 1);
		thread t2(simulatedAnnealing, n, a, b, deJong, 1);
		thread t3(simulatedAnnealing, n, a, b, deJong, 1);
		t1.join();
		t2.join();
		t3.join();
		cout << "Iteration " << i << " ended\n";
	}
	cout << "[DeJong] Size 30 completed! \n";
	cout << "[SIMULATED ANNEALING] Michalewicz function finished! \n";
}
//Here end DeJong SA

cout << "Preparing Schwefel function! \n";
{
	cout << "Preparing threads for Schwefel size 5 SA\n";
	n = 5;
	a = -500;
	b = 500;
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 10; ++i)
	{
		thread t1(simulatedAnnealing, n, a, b, Schwefel, 2);
		thread t2(simulatedAnnealing, n, a, b, Schwefel, 2);
		thread t3(simulatedAnnealing, n, a, b, Schwefel, 2);
		t1.join();
		t2.join();
		t3.join();
		cout << "Iteration " << i << " ended\n";
	}
	cout << "[Schwefel] Size 5 completed! \n";

	n = 10;
	cout << "Preparing threads for Schwefel size 10 first improvement\n";
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 10; ++i)
	{
		thread t1(simulatedAnnealing, n, a, b, Schwefel, 2);
		thread t2(simulatedAnnealing, n, a, b, Schwefel, 2);
		thread t3(simulatedAnnealing, n, a, b, Schwefel, 2);
		t1.join();
		t2.join();
		t3.join();
		cout << "Iteration " << i << " ended\n";
	}
	cout << "[Schwefel] Size 10 completed! \n";

	n = 30;
	cout << "Preparing threads for Schwefel size 30 first improvement\n";
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 10; ++i)
	{
		thread t1(simulatedAnnealing, n, a, b, Schwefel, 2);
		thread t2(simulatedAnnealing, n, a, b, Schwefel, 2);
		thread t3(simulatedAnnealing, n, a, b, Schwefel, 2);
		t1.join();
		t2.join();
		t3.join();
		cout << "Iteration " << i << " ended\n";
	}
	cout << "[Schwefel] Size 30 completed! \n";
	cout << "[SIMULATED ANNEALING] Schwefel function finished! \n";
}
//Here ends Schwefel

cout << "Preparing Rastrigin function!\n";
{
	cout << "Preparing threads for Rastrigin size 5 SA\n";
	n = 5;
	a = -5.12;
	b = 5.12;
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 10; ++i)
	{
		thread t1(simulatedAnnealing, n, a, b, Rastrigin, 3);
		thread t2(simulatedAnnealing, n, a, b, Rastrigin, 3);
		thread t3(simulatedAnnealing, n, a, b, Rastrigin, 3);
		t1.join();
		t2.join();
		t3.join();
		cout << "Iteration " << i << " ended\n";
	}
	cout << "[Rastrigin] Size 5 completed! \n";

	n = 10;
	cout << "Preparing threads for Rastrigin size 10 first improvement\n";
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 10; ++i)
	{
		thread t1(simulatedAnnealing, n, a, b, Rastrigin, 3);
		thread t2(simulatedAnnealing, n, a, b, Rastrigin, 3);
		thread t3(simulatedAnnealing, n, a, b, Rastrigin, 3);
		t1.join();
		t2.join();
		t3.join();
		cout << "Iteration " << i << " ended\n";
	}
	cout << "[Rastrigin] Size 10 completed! \n";

	n = 30;
	cout << "Preparing threads for Rastrigin size 30 first improvement\n";
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 10; ++i)
	{
		thread t1(simulatedAnnealing, n, a, b, Rastrigin, 3);
		thread t2(simulatedAnnealing, n, a, b, Rastrigin, 3);
		thread t3(simulatedAnnealing, n, a, b, Rastrigin, 3);
		t1.join();
		t2.join();
		t3.join();
		cout << "Iteration " << i << " ended\n";
	}
	cout << "[Rastrigin] Size 30 completed! \n";
	cout << "[SIMULATED ANNEALING] Rastrigin function finished! \n";
}

cout << "Preparing Michalewicz's function!\n";
{
	cout << "Preparing threads for Rastrigin size 5 SA\n";
	n = 5;
	a = 0;
	b = PI;
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 10; ++i)
	{
		thread t1(simulatedAnnealing, n, a, b, Michalewicz, 4);
		thread t2(simulatedAnnealing, n, a, b, Michalewicz, 4);
		thread t3(simulatedAnnealing, n, a, b, Michalewicz, 4);
		t1.join();
		t2.join();
		t3.join();
		cout << "Iteration " << i << " ended\n";
	}
	cout << "[Michalewicz] Size 5 completed! \n";

	n = 10;
	cout << "Preparing threads for Michalewicz size 10 first improvement\n";
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 10; ++i)
	{
		thread t1(simulatedAnnealing, n, a, b, Michalewicz, 4);
		thread t2(simulatedAnnealing, n, a, b, Michalewicz, 4);
		thread t3(simulatedAnnealing, n, a, b, Michalewicz, 4);
		t1.join();
		t2.join();
		t3.join();
		cout << "Iteration " << i << " ended\n";
	}
	cout << "[Michalewicz] Size 10 completed! \n";

	n = 30;
	cout << "Preparing threads for Michalewicz size 30 first improvement\n";
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 10; ++i)
	{
		thread t1(simulatedAnnealing, n, a, b, Michalewicz, 4);
		thread t2(simulatedAnnealing, n, a, b, Michalewicz, 4);
		thread t3(simulatedAnnealing, n, a, b, Michalewicz, 4);
		t1.join();
		t2.join();
		t3.join();
		cout << "Iteration " << i << " ended\n";
	}
	cout << "[Michalewicz] Size 30 completed! \n";
	cout << "[SIMULATED ANNEALING] Michalewicz function finished! \n";
}
cout << "SIMULATED ANNEALING FINISHED!";*/
return 0;
}