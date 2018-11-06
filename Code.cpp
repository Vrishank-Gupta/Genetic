#include <math.h>
#include <time.h>
#include <iostream>
#include <vector>
#include<unordered_map>
#include <fstream>
#include <limits.h>
#define ll long long
#define POP 100
#define WEIGHT 3030
#define FEA 23
#define TERMINATION 85
#define TOURNAMENT 20
#define CROSSRATIO 10
#define GENLIMIT 100
#define STEADY 1
#define GEN 0
#define MUTRATIO 10
using namespace std;


// Structures

// item structure
struct item {
	int w, v, q, s;
} listPop[23];

int genCount = 0;
// chromosome
struct chromosomes {
	vector<int> chromosome;
	int fitness;
	int weight;
	int gen;
	chromosomes(): chromosome(23, 0), fitness(0), weight(0), gen(0) {}
};

// roulette structure
class roulette {
	unordered_map <int, int> m;
	unordered_map <int, vector<chromosomes> > mapping;
public:
	void wrapperChromosome(chromosomes &c) {
		if (m.find(c.fitness) == m.end()) {
			addChromosome(c);
			m[c.fitness] = 1;
			return;
		}
		else {
			addChromosome(c);
			incr(c.fitness);
		}
	}

	void addChromosome(chromosomes &c) {
		mapping[c.fitness].push_back(c);
	}

	void incr(int key) {
		m[key]++;
	}


	pair<chromosomes, chromosomes> parents() {
		int size = m.size(), i = 0, flag = 0, f1, f2;
		int random1 = rand() % POP , random2 = rand() % POP;

		// roulette wheel pointer creation
		unordered_map <int, int> :: iterator prev = m.begin();
		for (unordered_map <int, int> :: iterator it = m.begin(); it != m.end(); ++it) {
			if (flag == 2)
				break;
			i += it->second;
			if (i >= random1 && !flag) {
				f1 = prev->first;
				flag++;
			}
			if (i >= random2) {
				f2 = prev->first;
				flag++;
			}
			prev = it;
		}

		// roulette wheel selection
		pair<chromosomes, chromosomes> p;
		flag = 0;
		// select parents
		for (unordered_map <int, vector<chromosomes> > :: iterator it = mapping.begin(); it != mapping.end(); ++it) {
			if (flag == 2)
				break;
			if (it->first == f1) {
				int r = rand() % it->second.size();
				p.first = it->second[r];
				flag++;
			}
			if (it->first == f2) {
				int r = rand() % it->second.size();
				p.second = it->second[r];
				flag++;
			}
		}
		return p;
	}

	// to check state of parent and offspring
	void display() {
		for (unordered_map <int, vector<chromosomes> > :: iterator it = mapping.begin(); it != mapping.end(); ++it) {
			cout << " fitness: " << it->first << " quantity: " << m[it->first] << endl;
		}
	}

};


// DEVELOPMENT FUNCTIONS

// to display chromosomes and to their state
void display(vector<chromosomes> &c) {
	for (int i = 0; i < c.size(); ++i) {
		for (int j = 0; j < FEA; ++j) {
			cout << c[i].chromosome[j] << " ";
		}
		cout << endl << "f: " << c[i].fitness << " w: " << c[i].weight << " g: " << c[i].gen << endl;
	}
}

// helper function to check if we have duplicate chromosomes
int checkDuplicate(vector<chromosomes> &c) {
	int ctr = 0;
	for (int i = 1; i < POP; ++i) {
		int j = 0;
		for (j = 0; j < FEA; ++j) {
			if (c[i - 1].chromosome[j] != c[i].chromosome[j])
				break;
		}
		if (j == FEA)
			ctr++;
	}
	return ctr;
}

// testing random generator
void randomGenerator() {
	for (int i = 0; i < 1000; ++i)
		cout << rand() % 10 + 1 << endl;
}

// wrapper function which can test any function using function pointer
void functionTester(void (*testFunction)()) {
	testFunction();
}


// HELPER FUNCTIONS

// comparator for sorting on basis of fitness
bool comparatorFitness(chromosomes a, chromosomes b) {
	return a.fitness > b.fitness;
}

// comparator for sorting on basis of age
bool comparatorAge(chromosomes a, chromosomes b) {
	return a.gen > b.gen;
}

// for randomly flipping a bit
bool coinToss() {
	return rand() % 2 == 1;
}

// provides avg fitness of a generation for plotting avg graph
float plotterAvg(vector<chromosomes> &c) {
	float  avg = 0;
	for (int i = 0; i < POP; ++i) {
		avg += c[i].fitness;
	}
	avg /= POP;
	return avg;
}

// provides best fitness of a generation for plotting best graph
int plotterBest(vector<chromosomes> &c) {
	int max = -1;
	for (int i = 0; i < POP; ++i) {
		if (max < c[i].fitness)
			max = c[i].fitness;
	}
	return max;
}

// MAIN FUNCTIONS

// set item value and weight
void setItem() {

	listPop[0].v = 500;
	listPop[0].w = 300;
	listPop[0].q = 5;
	listPop[0].s = 2;

	listPop[1].v = 300;
	listPop[1].w = 280;
	listPop[1].q = 3;
	listPop[1].s = 1;

	listPop[2].v = 700;
	listPop[2].w = 600;
	listPop[2].q = 5;
	listPop[2].s = 3;

	listPop[3].v = 650;
	listPop[3].w = 400;
	listPop[3].q = 5;
	listPop[3].s = 3;

	listPop[4].v = 400;
	listPop[4].w = 300;
	listPop[4].q = 5;
	listPop[4].s = 3;

	listPop[5].v = 900;
	listPop[5].w = 500;
	listPop[5].q = 7;
	listPop[5].s = 8;

	listPop[6].v = 700;
	listPop[6].w = 400;
	listPop[6].q = 5;
	listPop[6].s = 8;


	listPop[7].v = 650;
	listPop[7].w = 350;
	listPop[7].q = 3;
	listPop[7].s = 8;


	listPop[8].v = 300;
	listPop[8].w = 300;
	listPop[8].q = 3;
	listPop[8].s = 3;


	listPop[9].v = 390;
	listPop[9].w = 350;
	listPop[9].q = 4;
	listPop[9].s = 3;

	listPop[10].v = 450;
	listPop[10].w = 380;
	listPop[10].q = 6;
	listPop[10].s = 3;

	listPop[11].v = 100;
	listPop[11].w = 310;
	listPop[11].q = 2;
	listPop[11].s = 6;

	listPop[12].v = 180;
	listPop[12].w = 160;
	listPop[12].q = 4;
	listPop[12].s = 6;

	listPop[13].v = 500;
	listPop[13].w = 250;
	listPop[13].q = 5;
	listPop[13].s = 6;

	listPop[14].v = 100;
	listPop[14].w = 300;
	listPop[14].q = 3;
	listPop[14].s = 1;

	listPop[15].v = 120;
	listPop[15].w = 200;
	listPop[15].q = 3;
	listPop[15].s = 1;

	listPop[16].v = 50;
	listPop[16].w = 20;
	listPop[16].q = 2;
	listPop[16].s = 4;

	listPop[17].v = 200;
	listPop[17].w = 100;
	listPop[17].q = 4;
	listPop[17].s = 4;

	listPop[18].v = 400;
	listPop[18].w = 250;
	listPop[18].q = 5;
	listPop[18].s = 4;

	listPop[19].v = 800;
	listPop[19].w = 400;
	listPop[19].q = 7;
	listPop[19].s = 7;


	listPop[20].v = 500;
	listPop[20].w = 200;
	listPop[20].q = 5;
	listPop[20].s = 7;

	listPop[21].v = 300;
	listPop[21].w = 200;
	listPop[21].q = 4;
	listPop[21].s = 4;

	listPop[22].v = 600;
	listPop[22].w = 350;
	listPop[22].q = 6;
	listPop[22].s = 4;
	cout << "hey";

}

// make intial population
void make_population(vector<chromosomes> &c) {
	for (int i = 0; i < POP; ++i) {
		int w =  0, v = 0, s = 0, q = 0;
		while (w <= WEIGHT) {
			int random = rand() % FEA;
			if (c[i].chromosome[random] == 0) {
				w += listPop[random].w;
				q += listPop[random].q;
				s += listPop[random].s;
				c[i].chromosome[random] = 1;
				v += listPop[random].v;
			}
		}
		c[i].weight = w;
		c[i].fitness = v;
	}
}

// helper function to calculate fitness
void fitnessFunction(vector<chromosomes> &c) {
	int ctr = 0;
	for (int i = 0; i < POP; ++i) {
		int v = 0, w = 0, j = 0, in = 0;
		while (j < FEA) {
			for (; j < FEA && w <= WEIGHT; ++j) {
				if (c[i].chromosome[j]) {
					v += listPop[j].v;
					w += listPop[j].w;
					in = j;
				}
			}

			if (w > WEIGHT) {
				v -= listPop[in].v;
				w -= listPop[in].w;
				c[i].chromosome[in] = 0;
			}
		}
		c[i].fitness = v;
		c[i].weight = w;
	}
}

// selection operator - roulette
pair<chromosomes, chromosomes> selectionRoulette(vector<chromosomes> &c) {

	// create roulette
	roulette r;

	// select parents
	int key;
	for (int i = 0; i < POP; ++i) {
		r.wrapperChromosome(c[i]);
	}
	r.display();
	return r.parents();
}


// selection operator - tournament
pair<chromosomes, chromosomes> selectionTournament(vector<chromosomes> &c) {
	int k = TOURNAMENT, random;
	vector<chromosomes> selected(k);
	for (int i = 0; i < k; ++i) {
		random = rand() % POP;
		selected[i] = c[random];
	}
	sort(selected.begin(), selected.end(), comparatorFitness);
	pair<chromosomes, chromosomes> p;
	p.first = selected[0];
	p.second = selected[1];
	return p;
}


// crossover operator - one point
pair<chromosomes, chromosomes> crossoverOnePoint(pair<chromosomes, chromosomes> &p, int &gen) {
	pair<chromosomes, chromosomes> children;
	int random = rand() % FEA;
	for (int i = 0; i < FEA; ++i) {
		if (i >= random) {
			children.first.chromosome[i] = p.second.chromosome[i];
			children.second.chromosome[i] = p.first.chromosome[i];
		} else {
			children.first.chromosome[i] = p.first.chromosome[i];
			children.second.chromosome[i] = p.second.chromosome[i];
		}
	}
	for (int i = 0; i < FEA; ++i) {
		if (children.first.weight > WEIGHT) {
			children.first.fitness = -1;
			break;
		}
		if (children.first.chromosome[i]) {
			children.first.fitness += listPop[i].v;
			children.first.weight += listPop[i].w;
		}
	}
	if (children.first.weight > WEIGHT) {
		children.first.fitness = -1;
	}
	for (int i = 0; i < FEA; ++i) {
		if (children.second.weight > WEIGHT) {
			children.second.fitness = -1;
			break;
		}
		if (children.second.chromosome[i]) {
			children.second.fitness += listPop[i].v;
			children.second.weight += listPop[i].w;
		}
	}
	if (children.second.weight > WEIGHT) {
		children.second.fitness = -1;
	}
	children.first.gen = children.second.gen = gen;
	return children;
}

// crossover operator - two point
pair<chromosomes, chromosomes> crossoverTwoPoint(pair<chromosomes, chromosomes> &p, int &gen) {
	pair<chromosomes, chromosomes> children;
	int random1 = rand() % FEA;
	int random2 = rand() % FEA;
	if (random1 > random2) {
		int t = random1;
		random1 = random2;
		random2 = t;
	}
	for (int i = 0; i < FEA; ++i) {
		if (i >= random1 && i < random2) {
			children.first.chromosome[i] = p.second.chromosome[i];
			children.second.chromosome[i] = p.first.chromosome[i];
		} else {
			children.first.chromosome[i] = p.first.chromosome[i];
			children.second.chromosome[i] = p.second.chromosome[i];
		}
	}
	for (int i = 0; i < FEA; ++i) {
		if (children.first.weight > WEIGHT) {
			children.first.fitness = -1;
			break;
		}
		if (children.first.chromosome[i]) {
			children.first.fitness += listPop[i].v;
			children.first.weight += listPop[i].w;
		}
	}
	if (children.first.weight > WEIGHT) {
		children.first.fitness = -1;
	}
	for (int i = 0; i < FEA; ++i) {
		if (children.second.weight > WEIGHT) {
			children.second.fitness = -1;
			break;
		}
		if (children.second.chromosome[i]) {
			children.second.fitness += listPop[i].v;
			children.second.weight += listPop[i].w;
		}
	}
	if (children.second.weight > WEIGHT) {
		children.second.fitness = -1;
	}
	children.first.gen = children.second.gen = gen;
	return children;
}

// crossover operator - uniform
pair<chromosomes, chromosomes> crossoverUniform(pair<chromosomes, chromosomes> &p, int &gen) {
	pair<chromosomes, chromosomes> children;
	for (int i = 0; i < FEA; ++i) {
		if (coinToss()) {
			children.first.chromosome[i] = p.second.chromosome[i];
			children.second.chromosome[i] = p.first.chromosome[i];
		} else {
			children.first.chromosome[i] = p.first.chromosome[i];
			children.second.chromosome[i] = p.second.chromosome[i];
		}
	}
	for (int i = 0; i < FEA; ++i) {
		if (children.first.weight > WEIGHT) {
			children.first.fitness = -1;
			break;
		}
		if (children.first.chromosome[i]) {
			children.first.fitness += listPop[i].v;
			children.first.weight += listPop[i].w;
		}
	}
	if (children.first.weight > WEIGHT) {
		children.first.fitness = -1;
	}
	for (int i = 0; i < FEA; ++i) {
		if (children.second.weight > WEIGHT) {
			children.second.fitness = -1;
			break;
		}
		if (children.second.chromosome[i]) {
			children.second.fitness += listPop[i].v;
			children.second.weight += listPop[i].w;
		}
	}
	if (children.second.weight > WEIGHT) {
		children.second.fitness = -1;
	}

	children.first.gen = children.second.gen = gen;
	return children;
}

// mutation operator - single bit
void mutateSingleBit(pair<chromosomes, chromosomes> &p) {
	int random = rand() % 1000 + 1;
	if (random <= MUTRATIO * 10) {
		random = random % FEA;
		p.first.chromosome[random] ^= 1;
		if (p.first.chromosome[random] == 1) {
			if ((p.first.weight + listPop[random].w) > WEIGHT)
				p.first.chromosome[random] = 0;
			else {
				p.first.fitness += listPop[random].v;
				p.first.weight += listPop[random].w;
			}
		} else {
			p.first.fitness -= listPop[random].v;
			p.first.weight -= listPop[random].w;
		}
		p.second.chromosome[random] ^= 1;
		if (p.second.chromosome[random] == 1) {
			if ((p.second.weight + listPop[random].w) > WEIGHT)
				p.second.chromosome[random] = 0;
			else {
				p.second.fitness += listPop[random].v;
				p.second.weight += listPop[random].w;
			}
		} else {
			p.second.fitness -= listPop[random].v;
			p.second.weight -= listPop[random].w;
		}
	}
}

// wrapper function for three steps - selecting parent, cross-over, mutation
vector<chromosomes> scmWrapper(vector<chromosomes> &c, int &gen) {

	vector<chromosomes> children;

	// choice between steady state and genearation cross over
	int crossover = STEADY ? CROSSRATIO :  POP;

	for (int i = 0; i < crossover / 2; ++i) {
		//select()
		pair<chromosomes, chromosomes> parent = selectionTournament(c);

		//crossover()
		pair<chromosomes, chromosomes> child = crossoverUniform(parent, gen);

		//mutate()
		mutateSingleBit(child);
		if (child.first.fitness > 0)
			children.push_back(child.first);
		if (child.second.fitness > 0)
			children.push_back(child.second);
	}
	return children;
}

// survivor selection - age
void survivorSelectionAge(vector<chromosomes> &c, vector<chromosomes> &p) {

	// sort parent array according to fitness
	sort(c.begin(), c.end(), comparatorAge);

	// replace
	for (int i = POP - p.size(), k = 0; i < POP; ++i, k++) {
		c[i] = p[k];
	}

}

// survivor selection - fitness
void survivorSelectionFitness(vector<chromosomes> &c, vector<chromosomes> &p) {

	// sort parent array according to fitness
	sort(c.begin(), c.end(), comparatorFitness);

	// sort children array according to fitness
	sort(p.begin(), p.end(), comparatorFitness);

	// replace
	for (int i = POP - p.size(), k = 0; i < POP; ++i, k++) { // remove last k of a generation
		if (p[k].fitness < c[i].fitness)
			continue;
		c[i] = p[k];
	}

}

// terminating condition - generation number or iteration
bool terminateGenLimit(int &gen) {
	return gen > GENLIMIT;
}

// terminating condition - population
bool terminatePop(vector<chromosomes> &c) {
	unordered_map<int, int> m;
	int maxFitness = -1;
	for (int i = 0; i < POP; ++i) {
		if (c[i].fitness > maxFitness)
			maxFitness = c[i].fitness;
		if (m.find(c[i].fitness) == m.end())
			m[c[i].fitness] = 1;
		else
			m[c[i].fitness]++;
	}
	int max = 0, t;
	for (unordered_map<int, int> :: iterator it = m.begin(); it != m.end(); ++it) {
		if (it->second > max) {
			max = it->second;
			t = it->first;
		}
	}
	return maxFitness == t ? (max >= TERMINATION ? true : false) : false;
}

//  solution of given GA
int solutionGA(vector<chromosomes> &c) {
	int res = 0;
	for (int i = 0; i < POP; ++i) {
		if (res < c[i].fitness)
			res = c[i].fitness;
	}
	return res;
}

int main() {

	// declare chromosome
	vector<chromosomes> c(POP);


	// set item value weight -- later with a db
	setItem();

	// make initial population
	make_population(c);

	// calculate its fitness
	fitnessFunction(c);

	// variable to be used
	int gen = 0;
	int best;
	float avg;

	// files for best and avg graph CSV
	ofstream fa("average.csv");
	ofstream fb("best.csv");

	// main loop
	do {
		cout << "Generation " << gen + 1 << endl << endl << endl;

		// remove later
		// sort according to fitness
		//sort(c.begin(), c.end(),comparator);

		// select-crossover-mutate wrapper
		vector<chromosomes> children = scmWrapper(c, gen);

		// survivor selection
		survivorSelectionFitness(c, children);

		// to plot average
		avg = plotterAvg(c);
		fa << gen << "," << avg << endl;

		// to plot best
		best = plotterBest(c);
		fb << gen << "," << best << endl;

		// generation iterator
		gen++;

		// display current generation
		display(c);

		// terminating condition

		// population
		// if(terminatePop(c))
		// 	break;

		// generation limit
		// if(terminateGenLimit(gen))
		// 	break;

	} while (true);

	// closing files
	fa.close();
	fb.close();

	// solution
	cout << "Solution is " << solutionGA(c) << "g: " << gen - 1 << endl;

	return 0;
}