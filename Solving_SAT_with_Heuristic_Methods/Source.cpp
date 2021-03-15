#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <chrono>
#include <random>
using namespace std;

#define POP_SIZE 100
#define MUTATION_PROBABILITY 0.01
#define CROSSOVER_PROBABILITY 0.2

bool readCNF(int& variableNumber, int& clauseNumber, vector<vector<int>>& clause)
{
	ifstream fin("test.cnf");
	char comment[256];
	fin.getline(comment,256);
	cout << comment << '\n';
	string problem;
	char p;
	fin >> p >> problem;
	if (p != 'p' || problem != "cnf")
	{
		cout << "Bad file" << '\n';
		return 0;
	}
	fin >> variableNumber >> clauseNumber;
	int aux,count=0;
	clause.resize(clauseNumber);
	while(count<clauseNumber)
	{
		fin >> aux;
		if (aux == 0)
			count++;
		else
			clause[count].push_back(aux);
	}
	
	return 1;
}

int solutionScore(vector<char>& bits, int& variableNumber, int& clauseNumber, vector<vector<int>>& clause)
{
	int satisfiedClauses = 0;
	for (int i = 0; i < clauseNumber; ++i)
	{
		bool ok = 1;
		for(int j=0;j<clause[i].size();++j)
			if (clause[i][j] > 0) 
			{
				if (bits[clause[i][j]-1] == 1)
				{
					satisfiedClauses++;
					ok = 0;
					break;
				}
			}
			else if (bits[abs(clause[i][j])-1] == 0)
			{
				satisfiedClauses++;
				ok = 0;
				break;
			}
	}
	return satisfiedClauses;
}

double random01()
{
	return (double)((double)(rand() % 100000) / (double)100000);
}

void simulatedAnnealing(int& variableNumber, int& clauseNumber, vector<vector<int>>& clause)
{
	vector<char>bits;
	vector<int>selectCandidates;
	srand(clock());
	for (int i = 0; i < variableNumber; ++i)
	{
		bits.push_back(rand() % 2);
	}
	double temperature = 100;
	int max = solutionScore(bits, variableNumber, clauseNumber, clause);
	int minValue = max;
	cout << max<<'\n';
	int plateaucounter = 0;
	int plateauvalue = clauseNumber-max;
	while ( plateaucounter<100 && temperature > 0.0000001)
	{
		for (int i = 0; i < 10000; ++i)
		{
			int position = rand() % variableNumber;
			bits[position] = !bits[position];
			int aux= solutionScore(bits, variableNumber, clauseNumber, clause);
			if (aux > max)
				max = aux;
			else if (random01() < exp(-1 * abs(aux - max) / temperature))
			{
				max = aux;
			}
			else
				bits[position] = !bits[position];
		}
		temperature *= .95;
		int unsat = clauseNumber - max;
		if (unsat == 0)
		{
			cout << "Solution found:\n";
			break;
		}
		if (unsat < minValue)
			minValue = unsat;
		if (plateauvalue == unsat)
			plateaucounter++;
		else
		{
			plateaucounter = 0;
			plateauvalue =unsat;
		}
		cout <<plateaucounter<<' '<<unsat<<' '<<minValue<<' '<<temperature<<'\n';
	}
}

void printCurrentPop(vector<vector<char>>& population,int & variableNumber)
{
	cout << "POPULATION GENERATED:\n";
	for (int i = 0; i < POP_SIZE; ++i)
	{
		for (int j = 0; j < variableNumber; ++j)
		{
			cout << (bool)population[i][j];
		}
		//cout << '\n';
	}
}

void computeSolution(vector<vector<char>>& population, vector<int>& solutions, int& variableNumber, int& clauseNumber, vector<vector<int>>&clause)
{
	solutions.clear();
	for (int i = 0; i < POP_SIZE; ++i)
	{
		solutions.push_back(solutionScore(population[i],variableNumber,clauseNumber,clause));
		//cout << solutions[i] << ' ';
	}
}

double computeFitness(int& solution,int& clauseNumber)
{
	return solution;	
}

vector<vector<char>> selectNext(vector<vector<char>>& population, vector<int>& solutions,int & clauseNumber)
{
	vector<double> fitness;
	int T = 0;
	for (int i = 0; i < POP_SIZE; ++i)
	{
		fitness.push_back(computeFitness(solutions[i],clauseNumber));
		T += fitness[i];
	}
	//cout <<"Sum="<< T<<'\n';
	vector<double> probabilities;
	for (int i = 0; i < POP_SIZE; ++i)
	{
		probabilities.push_back(fitness[i] / T);
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
		for (int ii = 0; ii < population.size(); ++ii)
			if (q[ii] < rand01 && rand01 <= q[ii + 1])
			{
				newSolution.push_back(population[ii]);
				break;
			}
	}
	return newSolution;
}

void mutate(vector<vector<char>>& population,int& variableNumber)
{
	auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::mt19937_64 mt_rand(seed);
	std::uniform_real_distribution<double> subUnitValueGen(0.0, 1.0);
	for (int i = 0; i < POP_SIZE; ++i)
		for (int ii = 0; ii < variableNumber; ++ii)
		{
			double pMutation = subUnitValueGen(mt_rand);
			if (pMutation < MUTATION_PROBABILITY)
				population[i][ii] = !population[i][ii];
		}
}

void greedyMutate(vector<vector<char>>& population, int& variableNumber, vector<int>& solutions, int& clauseNumber, vector<vector<int>>& clause)
{
	for (int i = 0; i < POP_SIZE; ++i)
	{
		for (int ii = 0; ii < variableNumber; ++ii)
		{
			population[i][ii] = !population[i][ii];
			if (solutions[i] > solutionScore(population[i], variableNumber, clauseNumber, clause))
				population[i][ii] = !population[i][ii];
		}
	}
}

void multiGreedyMutate(vector<vector<char>>& population, int& variableNumber, vector<int>& solutions, int& clauseNumber, vector<vector<int>>& clause)
{
	
	for (int i = 0; i < POP_SIZE; ++i)
	{
		int counter = 0;
		bool imprv = 1;
		while (imprv != 0)
		{
			for (int ii = 0; ii < variableNumber; ++ii)
			{
				population[i][ii] = !population[i][ii];
				auto aux = solutionScore(population[i], variableNumber, clauseNumber, clause);
				if (solutions[i] >= aux)
					population[i][ii] = !population[i][ii];
				else {
					solutions[i] = aux;
					imprv = 0;
				}
					
			}
			if (imprv == 0)
			{
				imprv = 1;
				counter++;
			}
			else imprv = 0;
		}
	}
}

template <typename A, typename B>
void zip(vector<A>& a, vector<B>& b, vector<pair<A, B>>& zipped)
{
	for (auto i = 0; i < a.size(); ++i)
		zipped.emplace_back(a[i], b[i]);
}

template <typename A, typename B>
void unzip(vector<pair<A, B>> & zipped, vector<A> & a, vector<B> & b)
{
	for (auto i = 0; i < a.size(); ++i)
	{
		a[i] = zipped[i].first;
		b[i] = zipped[i].second;
	}
}

void crossover(vector<vector<char>> & population,int&variableNumber)
{
	auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::mt19937_64 mt_rand(seed);
	std::uniform_real_distribution<double> subUnitValueGen(0.0, 1.0);
	vector<double> probabilities;
	for (int i = 0; i < POP_SIZE; ++i)
		probabilities.push_back(subUnitValueGen(mt_rand));
	vector<pair< double, vector<char>>> zipped;
	zip(probabilities, population, zipped);
	sort(begin(zipped), end(zipped),
		[&](auto & a, auto & b)
		{
			return a.first < b.first;
		}
	);
	unzip(zipped, probabilities, population);
	std::uniform_int_distribution<int> ValueGen(1, variableNumber - 1);
	for (int i = 0; i < POP_SIZE; i = i + 2)
	{
		if (probabilities[i + 1] < CROSSOVER_PROBABILITY)
		{
			int cutpoint = ValueGen(mt_rand);
			for (int ii = cutpoint; ii < variableNumber; ++ii)
				swap(population[i][ii], population[i + 1][ii]);
		}
		else
			break;
	}
}

void uniformCrossover(vector<vector<char>>& population, int& variableNumber)
{
	auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::mt19937_64 mt_rand(seed);
	std::uniform_real_distribution<double> subUnitValueGen(0.0, 1.0);
	vector<double> probabilities;
	for (int i = 0; i < POP_SIZE; ++i)
		probabilities.push_back(subUnitValueGen(mt_rand));
	vector<pair< double, vector<char>>> zipped;
	zip(probabilities, population, zipped);
	sort(begin(zipped), end(zipped),
		[&](auto & a, auto & b)
		{
			return a.first < b.first;
		}
	);
	unzip(zipped, probabilities, population);
	std::uniform_int_distribution<int> ValueGen(1, variableNumber - 1);
	for (int i = 0; i < POP_SIZE; i = i + 2)
	{
		for (int ii = 0; ii < variableNumber; ++ii)
		{
			if(random01()<0.5)
			swap(population[i][ii], population[i + 1][ii]);
		}
	}
}

int computeMin(vector<int>& solutions)
{
	int aux = INT_MIN;
	for (int i = 0; i < solutions.size(); ++i)
		if (aux < solutions[i])
			aux = solutions[i];
	return aux;
}

void GA(int& variableNumber, int& clauseNumber, vector<vector<int>>& clause)
{
	int t = 0;
	int SAT;
	vector<char>bits;
	vector<vector<char>>population;
	srand(clock());
	for (int i = 0; i < POP_SIZE; ++i)
	{
		bits.clear();
		for (int ii = 0; ii < variableNumber; ++ii)
		{
			bits.push_back(rand() % 2);
		}
		population.push_back(bits);
	} /*Random population initialization*/
	//printCurrentPop(population,variableNumber);
	vector<int> solutions;
	int min = 0;
	computeSolution(population, solutions, variableNumber, clauseNumber,clause);
	while (t < 1000)
	{
		++t;
		auto aux=selectNext(population,solutions,clauseNumber);
		population.clear();
		for (int i = 0; i < aux.size(); ++i)
			population.push_back(aux[i]);
		//crossover(population, variableNumber);
		uniformCrossover(population, variableNumber);
		//mutate(population,variableNumber);
		//greedyMutate(population, variableNumber);
		multiGreedyMutate(population, variableNumber, solutions,clauseNumber,clause);
		computeSolution(population, solutions, variableNumber, clauseNumber, clause);
		if ((SAT = computeMin(solutions)) == clauseNumber)
		{
			cout << "Solution found\n";
			break;
		}
		else 
		{
			cout <<t<<" "<<SAT <<" "<<min<<"\n";
			if (SAT > min)
				min = SAT;
		}
	}
}

int main()
{
	int variableNumber, clauseNumber;
	vector<vector<int>> clause;
	readCNF(variableNumber, clauseNumber, clause);
	cout << "Number of variables: " << variableNumber << ", Number of clauses: " << clauseNumber <<  '\n';
	//simulatedAnnealing(variableNumber, clauseNumber, clause);
	GA(variableNumber, clauseNumber, clause);
}