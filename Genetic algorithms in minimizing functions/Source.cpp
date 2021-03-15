#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <chrono>
#include <fstream>
//#include <thread>
using namespace std;
#define POP_SIZE 200
#define PI 3.14159265
double MUTATION_PROBABILITY=0.01; 
double CROSSOVER_PROBABILITY=0.2;
double a, b;
int prec;
int L;
int func;
int n;
ofstream myfile("Michalewicz_10_200pop.txt");

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
		result += (sin(values[i - 1]) * pow(sin(i * values[i - 1] * values[i - 1] / PI), 20));
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

void printCurrentPop(vector<vector<char>>& population)
{
	cout << "POPULATION GENERATED:\n";
	for (int i = 0; i < POP_SIZE; ++i)
	{
		for (int j = 0; j < L; ++j)
		{
			cout << (bool)population[i][j];
		}
		cout << '\n';
	}
}

void computeSolution(vector<vector<char>>& population, vector<double>& solutions, int l,int n, double(*func)(vector<double>&))
{
	solutions.clear();
	for (int i = 0; i < POP_SIZE; ++i)
	{
		auto aux = decode(population[i], l, n, a, b);
		solutions.push_back(func(aux));
	}
}

double computeFitness(double& a)
{
	if (func == 1)
		return pow(a, -3);
	else if (func == 2)
		return pow((double)((double)(500 * n + a) / (double)(500 * n * 2)), -10);
	else if (func == 3)
		//return pow((double)((double)(a) / (double)(46.2144*n)), -30); /worse results
		//return pow(a,-0.1-0.2*a) //
		return pow(a, -10);//-0.1-0.2*
	else if (func == 4)
		return pow(exp(a),-10);
}

void computeMin(vector<double>& solutions)
{
	double aux = DBL_MAX;
	for (int i = 0; i < solutions.size(); ++i)
		if (aux > solutions[i])
			aux = solutions[i];
	myfile << aux << '\n';
	//cout << aux << '\n';
}

vector<vector<char>> selectNext(vector<vector<char>>& population,vector<double>& solutions)
{
	vector<double> fitness;
	double T = 0;
	for (int i = 0; i < POP_SIZE; ++i)
	{
		fitness.push_back(computeFitness(solutions[i]));
		//cout <<"fitness="<<fitness[i] <<" value="<<solutions[i]<< '\n';
		T += fitness[i];
	}
	//cout <<"Sum="<< T<<'\n';
	vector<double> probabilities;
	for (int i = 0; i < POP_SIZE; ++i)
	{
		probabilities.push_back(fitness[i] / T);
		//cout << probabilities[i]<<'\n';
	}
	vector<double> q;
	q.push_back(0);
	for (int i = 0; i < POP_SIZE; ++i)
	{
		q.push_back(q[i] + probabilities[i]);
	}
	auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::mt19937_64 mt_rand(seed);
	std::uniform_real_distribution<double> subUnitValueGen(0.0, 1.0);
	vector < vector <char>> newSolution;
	for (int i = 0; i < POP_SIZE; ++i)
	{
		double rand01 = subUnitValueGen(mt_rand);
		//cout << "R01=" << rand01<<'\n';
		for (int ii = 0; ii < population.size();++ii)
			if (q[ii] < rand01 && rand01 <= q[ii + 1])
			{
				newSolution.push_back(population[ii]);
				break;
			}
	}
	return newSolution;
}

void mutate(vector<vector<char>>& population)
{
	auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::mt19937_64 mt_rand(seed);
	std::uniform_real_distribution<double> subUnitValueGen(0.0, 1.0);
	for (int i = 0; i < POP_SIZE; ++i)
		for (int ii = 0; ii < L; ++ii)
		{
			double pMutation = subUnitValueGen(mt_rand);
			if (pMutation < MUTATION_PROBABILITY)
				population[i][ii] = !population[i][ii];
		}
}

template <typename A, typename B>
void zip(vector<A>&a,vector<B>&b,vector<pair<A,B>>&zipped)
{
	for (auto i = 0; i < a.size(); ++i)
		zipped.emplace_back(a[i], b[i]);
}

template <typename A, typename B>
void unzip(vector<pair<A, B>>& zipped, vector<A>& a, vector<B>& b)
{
	for (auto i = 0; i < a.size(); ++i)
	{
		a[i] = zipped[i].first;
		b[i] = zipped[i].second;
	}
}

void crossover(vector<vector<char>>& population)
{
	auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::mt19937_64 mt_rand(seed);
	std::uniform_real_distribution<double> subUnitValueGen(0.0, 1.0);
	vector<double> probabilities;
	for (int i = 0; i < POP_SIZE; ++i)
		probabilities.push_back(subUnitValueGen(mt_rand));
	vector<pair< double,vector<char>>> zipped;
	zip(probabilities,population, zipped);
	sort(begin(zipped), end(zipped),
		[&](auto & a, auto & b)
		{
			return a.first < b.first;
		}
	);
	unzip(zipped, probabilities, population);
	std::uniform_int_distribution<int> ValueGen(1, L - 1);
	for (int i = 0; i < POP_SIZE; i=i+2)
	{
		if (probabilities[i + 1] < CROSSOVER_PROBABILITY)
		{
			int cutpoint = ValueGen(mt_rand);
			for (int ii = cutpoint; ii < L; ++ii)
				swap(population[i][ii], population[i + 1][ii]);
		}
		else
			break;
	}
}

void Genetics(int l,int n, double(*func)(vector<double>&))
{
	clock_t start, end;
	start = clock();
	int t = 0;
	L = l * n;
	srand(clock());
	//Generate the starting population
	vector<char> bits;
	vector<vector<char>> population;
	for (int i = 0; i < POP_SIZE; ++i)
	{
		bits.clear();
		for (int ii = 0; ii < L; ++ii)
		{
			bits.push_back(rand() % 2);
		}
		population.push_back(bits);
	}
	//printCurrentPop(population);
	vector<double> solutions;
	computeSolution(population, solutions, l, n,func);
	//computeMin(solutions);
	while (t < 1000)
	{
		++t; 
		auto aux=selectNext(population,solutions);
		population.clear();
		for (int i = 0; i < aux.size(); ++i)
			population.push_back(aux[i]);
		mutate(population);
		crossover(population);
		computeSolution(population, solutions, l, n,func);
		if (t % 100 == 0)
		{
			//computeMin(solutions);
			cout << t << '\n';
		}
	}
	computeMin(solutions);
	end = clock();
	myfile << ' ' << (double)((double)(end - start) / (double)CLOCKS_PER_SEC) << " seconds.\n";
}

int main()
{
	cout << "TEST" << '\n';
	int l;
	func = 4; ///1=DJ 2=Sch 3=Rastrigin 4=Mich
	//n = 30;
	//n = 10;
	a = 0;
	b = PI;
	prec = 3;
	//ofstream myfile("Schwefel_5.txt");
	n = 10;
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 30; ++i)
	{
		Genetics(l, n, Michalewicz);
		cout << "Iteration " << i << " finished.\n";
	}
	myfile.close();

	/* myfile=ofstream("Schwefel_30X2.txt");
	n = 30;
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 30; ++i)
	{
		Genetics(l, n, Schwefel);
		cout << "Iteration " << i << " finished.\n";
	}
	myfile.close();

	 myfile = ofstream("Schwefel_30X3.txt");
	n = 30;
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 30; ++i)
	{
		Genetics(l, n, Schwefel);
		cout << "Iteration " << i << " finished.\n";
	}
	myfile.close();

	myfile = ofstream("Schwefel_30X4.txt");
	n = 30;
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 30; ++i)
	{
		Genetics(l, n, Schwefel);
		cout << "Iteration " << i << " finished.\n";
	}
	myfile.close();


	myfile = ofstream("Schwefel_30X5.txt");
	n = 30;
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 30; ++i)
	{
		Genetics(l, n, Schwefel);
		cout << "Iteration " << i << " finished.\n";
	}
	myfile.close();


	myfile = ofstream("Schwefel_30X6.txt");
	n = 30;
	l = computeComponentLength(a, b, prec);
	for (int i = 0; i < 30; ++i)
	{
		Genetics(l, n, Schwefel);
		cout << "Iteration " << i << " finished.\n";
	}
	myfile.close();*/
}