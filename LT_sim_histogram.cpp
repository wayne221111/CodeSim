#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <math.h>
#include "randomc.h"
#include "ConvoCode.h"
#include "LTCode.h"
#include "statistics.h"
#include <omp.h>
#include <algorithm>
#include <unistd.h>

using namespace std;
using namespace CodeSim;

int K;
long Run;
double Delta;
int STEPS;

long histo_bins[16] = {0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384};
const int N_histo_bins = 16;
vector< vector<double> > histoErrorCount;
vector< vector< vector< vector<int> > > > histoErrorPattern;

int Dsize;
int* Degree;
int g_seed = (int)time(0);
CRandomMersenne Rnd(g_seed);

double* D;
double* SD;

double* normolize(double* d){
        int i;
        double z = 0;
        for(i=0;i<Dsize;i++){
                if(d[i]<0) d[i] = -d[i];
                z = z + d[i];
        }
        for(i=0;i<Dsize;i++) d[i] = d[i]/z;
        return d;
}

int main(int argn, char **args){
        if(argn < 4) {
                cerr << "Usage: histogram.out LT_distibution_file convo_code_file interleaver_file [histogram_file]"<<endl;
                exit(1);
        }
        ifstream dist_file(args[1]);
        if(dist_file.fail()){
                cerr << "File can not be opened: " << args[1] <<endl;
                exit(1);
        }
        string histo_filename, pattern_filename;
        if(argn == 4) {
                histo_filename = args[1];
                histo_filename += "_histo.txt";
                pattern_filename = args[1];
                pattern_filename += "_pattern";
 }
        else    histo_filename = args[4];
        ofstream histo_file(histo_filename.c_str());
        if(histo_file.fail()) {
                cerr << "File can not be opened: " << histo_filename <<endl;
                exit(1);
        }
        ofstream pattern_file(pattern_filename.c_str());
        if(pattern_file.fail()) {
                cerr << "File can not be opened: " << pattern_filename <<endl;
                exit(1);
        }
        cout << "debug -1" << endl;
        nice(20);
        dist_file >> K;
        K = 80640;
        Run = 100;//(1000000000LL/K);
        dist_file >> Dsize;
        Degree = new int[Dsize];
        D = new double[Dsize];
        SD = new double[Dsize];
        for (int i=0; i<Dsize; i++) {
                dist_file >> Degree[i];
        }
        for (int i=0; i<Dsize; i++) {
                dist_file >> D[i];
        }
        STEPS = 10;
        Delta = 0.005;
        histoErrorCount.assign(STEPS, vector<double>() ); //[step][histobin]
        cout << "debug" << endl;
        histoErrorPattern.assign(STEPS, vector< vector< vector<int> > >(N_histo_bins, vector< vector<int> >() ) );
        {
                histo_file << "tab of degree\t";
                for (int i =0; i<Dsize; i++) {
                        histo_file << Degree[i] << '\t';
                }
                histo_file << "\ndistribution\t";
                for (int i =0; i<Dsize; i++) {
                        histo_file << D[i] << '\t';
                }
                histo_file << "\nEpsilons\t";
                for (int i = 0; i<STEPS; i++) {
                        histo_file << i*Delta << '\t';
                        histoErrorCount[i].assign(N_histo_bins,0);
                }
                histo_file << '\n';

                double histo_bins_ratio[N_histo_bins];

                if(histo_bins[N_histo_bins-1]<K)
                        histo_bins[N_histo_bins-1]=K;
                for (int i=0; i<N_histo_bins; i++) {
                        histo_bins_ratio[i] = histo_bins[i]/double(K);
                        if (histo_bins_ratio[i]>1) {
                                histo_bins_ratio[i]=1;
                        }
                }
cout<<"debug 1" << endl;
                int start_time = time(0);

        //      #pragma omp parallel for schedule(dynamic) num_threads(PARALLEL_THREADS)
                for(long run=0;run<Run;run++){
                        int seed;
        //              #pragma omp critical
                        seed = Rnd.BRandom();
cout <<"debug 2"<< endl;

for (int i = 0; i< STEPS; i++) {
                        Codeword<Bit> a;
                        a.reserve(80640);
                        for (int t=0; t<80640; t++) a.push_back(1);
                        Codeword<Byte> b = BitToByteCoverter::convert(a);
cout <<"debug 2.9" << endl;
                                LT_sim<Byte> sim(b.size(), b.size()*(1+Delta*i), Dsize, Degree, D, seed);
cout <<"debug 2.91" <<endl;
        Codeword<Byte> c = sim.encode(b);
cout << "debug 2.92" << endl;
c = sim.decode(c);
cout << "debug 2.93" << endl;
                                vector<int>     tmp;
cout << "debug 3" << endl;
                                for (int k=0; k<c.size(); k++) if(c[k].isErased() ) tmp.push_back(k);
                                double t = sim.failureRate();
                                for (int h=0; h<N_histo_bins; h++) {
                                        if (t<=histo_bins_ratio[h]+(1.0/K/10)) {
        //                                      #pragma omp atomic
                                                histoErrorCount[i][h]++;
                                                histoErrorPattern[i][h].push_back(tmp);
                                                break;
                                        }
                                }
                        }
                }
                //histo_file<<"Histogram";
                for (int i = 0; i<N_histo_bins; i++) {
                        histo_file << histo_bins_ratio[i]<<'\t';
                        for(int j = 0; j< STEPS; j++)
                                histo_file <<  histoErrorCount[j][i] / (double)Run<< '\t';
                        histo_file << '\n';
                }
                cout << "debug 4" << endl;
                for(int j = 0; j< STEPS; j++) {
                        pattern_file << "Epsilon: " << 1+Delta*j << endl;
                        for (int i = 0; i<N_histo_bins; i++) {
                                pattern_file << '\t' << histo_bins_ratio[i] << endl;
                                for(int k = 0; k<histoErrorPattern[j][i].size(); k++) {
                                        pattern_file << '\t' << '\t';
                                        for(int l = 0; l<histoErrorPattern[j][i][k].size(); l++)
                                                pattern_file << histoErrorPattern[j][i][k][l] << "  ";
                                        pattern_file << endl;
                                }
                        }
                }
                histo_file << "Time\t" << time(0) - start_time<< "\tK\t"<< K << "\tRun\t"<< Run <<endl;
        }
/*      ConvoCode cc(args[2]);
        int Layer = cc.getK();
        unsigned long L = 80000*6/7/Layer;
        Run = 10000//MAX_BIT / L;
        Permutator<Bit> inter(args[3], true);

        Codeword<Bit> a;
        a.reserve(Layer*L);
        for (int t=0; t<Layer*L; t++)   a.push_back(Rnd.IRandomX(0, 1));
        Codeword<Bit> b = cc.encode(a);
        b = inter.permutate(b);

        for (int s=0; s<STEPS; s++) {
                for(int j=1; j<N_histo_bins; j++) {
                        vector<unsigned long> total_err(Layer, 0);
                        double Epsilon = (1+s*Delta) -1;
                        vector<vector<bool> > histo_mask(Layer, vector<bool>(L+1,0));
                        vector<map<unsigned long, unsigned long> > histo(Layer, map<unsigned long, unsigned long>());
                        {
                                #pragma omp parallel for schedule(dynamic) num_threads(PARALLEL_THREADS)
                                for (int i=0; i<Run; i++) {
                                        vector<unsigned long> err(Layer, 0);
                                        Codeword<Byte> b1 = BitToByteCoverter::convert(b);
                                        int tmp = histoErrorPattern[s][j].size();
                                        tmp = Rnd.IRandomX(0, tmp-1);
                                        tmp = histoErrorPattern[s][j][tmp].size();
                                        for (int k=0; k<tmp; k++) {
                                                int index = histoErrorPattern[s][j][tmp][k];
                                                b1[index].setErased(True);
                                        }
                                        Codeword<Bit> c1 = BitToByteCoverter::revert(b1);
                                        c1 = inter.depermutate(c1);
                                        Codeword<Bit> c = cc.decode(c1);
                                        for (int k=0; k<c.size(); k++) {
                                                if (!(a[k] == c[k])) {
                                                        err[k%Layer]++;
                                                }
                                        }
                                        for (int k=0; k<Layer; k++) {
                                                #pragma omp atomic
                                                total_err[k]+=err[k];
                                        }
                                }
                        }
                        vector<double> avr_err_rate(Layer, 0);
                        for (int i=0; i<Layer; i++)
                                avr_err_rate[i] = total_err[i]/L/Run;
                }
        } */
        return 0;
}

