/*
 * File: main_Sedus_14_03_12.cpp
 * Author: Diego + Oriol + Marina
 * Date: 2014_03_12
 */

/*ARGUMENTS:
 * SimulationID R C Supertime TimeToFixation SC/WR/HS
 * If timeToFixation=0, a random fixation trajectory is simulated for each run (av. time = 4N),
 * otherwise, timeToFixation indicates the speed of fixation (time = timeToFixation * N) and is
 * linear trajectory towards fixation
 */


/*CHANGES:
 * printing has been reduced to the minimum and statistics for whole block have been eliminated
 */

/* PARTICULAR CHARACTERISTICS OF THIS MAIN
  */
//HOLA
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <math.h>
#include <cstdlib>
#include <time.h>
#include <ctime>
#include <algorithm>

#include <cstddef>
#include <array>
#include <vector>
#include <iterator>
#include <stdio.h>
#include <stdlib.h>
#include "sedus.h"
//#include "mainwindow.h"
#include <QObject>
#include <QString>
#include <unistd.h>
#include <QProgressBar>


using namespace std;

////////////////////////////////////////
//// SIMULATION PRINCIPAL VARIABLES ////
////////////////////////////////////////

int N; // Population size
int PROMETHEUS = 1000; // Number of generations for each genealogy
#define SUPERTIME //50 // Number of simulations per execution
int BLOCKLENGTH = 5000; // Block length
//#define BLOCKLENGTH 5000
int SAMPLE = 50; // Sample size
#define MUTTABLESIZE 1000000 // Maximum number of mutations (size of muttable)
#define B 3 // Maximum number of blocks per chromosome
#define BIGTIME 130 // BIGTIME * N is the total number of generations simulated per run
#define numOfBins 5 // Number of segments in which we divide each block (exclusively used for analysis)
#define maxNumOfHS 10  // Maximum number of crossover hotspots
#define numofsamples 1 // Number of different samples taken per run

int TIMELENGTH; // Total number of generations (including all phases)
int BURNIN; // Number of generations in phase I
int STRUCTURED; // Number of generations in phase II

float THETA = 0.001; // Population scaled mutation rate: THETA = 4 N mu
float R = 10; // Population scaled crossover rate: R = 4 N rho*BLOCKLENGTH
float C = 0.5; // Population scaled gene conversion rate: C = 4N*kappa*lambda

int numHS = 1; // Number of hotspots
int crossoverBegin[maxNumOfHS]; // Start point of crossover regions
int crossoverEnd[maxNumOfHS]; // End point of crossover regions


string letter = ""; // Simulation ID
string str;
string dir;

////////////////////////////////
//// DECLARATION OF CLASSES ////
////////////////////////////////

struct chrom { // Chromosomes

    int b; // Number of blocks
    int mpb[B]; // Number of mutations per block
    //int mutation[B][BLOCKLENGTH]; // Array of mutation positions (per block)
    std::vector<std::vector<int>> mutation;
    chrom(){
        mutation.resize(B);
        for(unsigned int i=0;i<mutation.size();i++){mutation[i].resize(BLOCKLENGTH);}
    }

};

struct mutation { // Mutations
//public:
    int position;
    int block;
    float frequency;
};

struct prev_pres {
    int prev;
    int pres;
};
struct fertility_info{
    int x;
    int y;
    bool recombinatrix;
};

typedef QVector<QVector<int>> qvint;
typedef QVector<QVector<double>> qvdouble;
///////////////////////////////
//// VARIABLES DECLARATION ////
///////////////////////////////

//// PARAMETERS OF EVENT ////
float mu; // Mutation rate per nucleotide and generation
float rho; // Crossover rate per nucleotide and generation
float meanTractLength = 100; // Mean Gene Conversion Tract Length (lambda)
float kappa; // Gene Conversion Initiation Rate
float numbptotalidentity = 0; // Length of the 100% identity tract
float IdentityInConvTract = 0; // Percent of identity required for conversion in all conversion tract
float donorRatio = 0.5; // Percentage of gene conversion events that occur from the original to the duplicated block

//// GENEALOGICAL MATRICES ////

std::vector<std::vector<int>> ancestry; // Matrix codifying the ancestor of each chromosome in each generation
std::vector<fertility_info> fertility_list;
std::vector<fertility_info> fertility_list_ini;
std::vector<std::vector<bool>> fertility;// Boolean Matrix codifying if a chromosome has descendants in the final generation of each era (each era is equivalent to PROMETHEUS generations) or not
std::vector<std::vector<bool>> fertility_ini;//Matrix fertility initialization
std::vector<std::vector<bool>> recombimatrix;// Boolean Matrix codifying if a chromosome comes from a recent recombination or not (has one or two parents)
std::vector<std::vector<int>> duplicontent;// Array indicating which chr carry the duplication (with present and previous lines)
//int ancestry[2 * N][PROMETHEUS];
//bool fertility[2 * N][PROMETHEUS];
//bool fertility_ini[2 * N][PROMETHEUS];
//bool recombimatrix[2 * N][PROMETHEUS];
//int duplicontent[2][2 * N];

//// GLOBAL VARIABLES AND STATISTICAL QUANTITIES ////
//chrom table[4 * N], *pointer[2][2 * N]; // GENOMIC INFORMATION
std::vector<chrom> table;
std::vector<std::vector<chrom*>> pointer;
int era, run; // t = generation inside each era, era = PROMETHEUS generations inside each phase (from 0 to BIGTIME-1), run = number of simulation runs (from 0 to SUPERTIME-1)
//int fixationTrajectory[20 * N + 1]; // Absolute frequency of the duplication in each generation of fixation process
std::vector<int> fixationTrajectory;
int timeToFixation;
int argc;
int correctArguments;
int sup;
int maxiterations;
int iterationspersupertime;
float perc;
float iterperc;
int phaseIiter;
int phaseIIiter;
int phaseIIIiter;

bool duFreq; // Duplication has occurred or not
//bool multihit[BLOCKLENGTH]; // Record the positions in which a mutation has occurred
std::vector<bool> multihit;
int duplicationFreq; // Absolute frequency of the duplication in the present generation

//// GLOBAL INTERNAL VARIABLES FOR FSL ////
int MutCount; // Total number of mutations segregating (the fixed ones are erased) in each moment
struct mutation muttable[MUTTABLESIZE], temporalmuttable[MUTTABLESIZE]; // Register of all the mutations

//// SAMPLING ////
//int sample[2 * N]; // Randomly sampled individuals
std::vector<int> sample;
int sampleN[] = {SAMPLE}; // Size of the samples (only one sample in this case)

//// FILES ////
ofstream profile;// N PROMETHEUS SUPERTIME BLOCKLENGTH SAMPLE MUTTABLESIZE B BIGTIME BURNIN/N STRUCTURED/N mu rho kappa meanTractLength
ofstream auxx; // NumOfFertileIndividuals, NumOfRecEvents, endTime of the Trajectory, NumOfMutEvents, NumOfMultihitCounts, NumOfConvEvents, NumOfFixationEvents (total and for each block)
ofstream samplefile[B + 2][numofsamples][2]; // For each block and the collapsed. 0 = pi; 1 = S. For each sample
ofstream mutationsNewFile[B]; // new mutation file, with ms-like format
ofstream SFS[B+2]; // site frequency spectra for each block + collapsed only for last era

///////////////////
//// FUNCTIONS ////
///////////////////

struct prev_pres phaseI();
int phaseII(int,int,int);
void phaseIII();

void open_files(); // Opening files
void close_files(); // Closing files

void genealogy(float, int); // rho, 0/1(non structured or structured) //// Filling genealogy matrices in each PROMETHEUS
void parentpicking(int[maxNumOfHS], int[maxNumOfHS], int, int, int,int,int); // crossoverBegin, crossoverEnd //// Create new generation from previous one (with recombination)

void duplication(int,int); // Create Duplication for eva (first duplicated chromosome)
void mutation(float, int, int); // For each fertile chromosome decide if a mutation happens and execute it if necessary
void conversion(float, int, int, float); // Only when duFreq is true. For each fertile chromosome decide if conversion happens and execute it if necessary

void statistics(int, bool, int); // Execution of all the statistic calculations (for each Era)
void FSL(int); // Count of All the Segregating, Fixed, Lost and Shared sites

void copychr(int, int, int, int); // Copy a chromosome (when there is no recombination)
int location(int, int, int, int); // Returns the location of a mutation or point in the mutation vector of a chromosome

void EraseFixedMutations(int, int, int); // Erase fixed mutations from chromosomes that have them (not necessary to consider them any more)

//float * SiteFrequencySpectrum(int, int, int); // Calculate pi and S values for a block
float * SiteFrequencySpectrumPrint(int, int, int, bool); // Calculate pi and S values for a block
//float * SiteFrequencySpectrumByBins(int, int, int); // Calculate pi and S values for a det block (by a det number of bins (if one = whople block))
float * SiteFrequencySpectrum_02(int, int, bool); // Calculate pi and S values for the collapsed block (0+2)
float * SiteFrequencySpectrum_02_Calling(int, int, bool); // Calculates pi and S for the collapsed sample doing the collapsed calling

void DivergenceForAll(int, int, int); // Average divergent positions between two blocks in the same chromosome
//float d_b_fromSample(int, int, int, int); // Average divergent positions between two blocks in different chromosomes

int DupliFreq(int, int, int); // Calculate Duplication frequency in the population

int muFrequencyIntWholePopAndSample(int, int, int, int); // Calculate absolute frequency of a given mutation
int muFrequencySampleIntDiscont(int, int, int, int, int, int); // Calculate absolute frequency of a given mutation but not in the whole population (sampling)
int muFrequencyCollapsedCallingFromSample(int, int, int); // Calculates absolute frequency of a given mutation collapsing both blocks from each individual

int SearchMutation(int, int, int); // Search a given mutation inside muttable (returns its position)
int tractpql(float); // Returns a given random tract length from the mean tract length

void SamplingIndividuals(int); // Sample corresponding number of individuals (register their number in sample[])

int GenerateFixationTrajectory(int, int); // Generates fixation trajectory of the duplication (fixed or not)
;
//void initialize_fertility_matrix(bool[2*N][PROMETHEUS]);
//void printMutationsFromSample(int, int, int); // Print Mutations of a sample in mutations file
void print_fertility();

float round(float, int);

int minim(int n1, int n2) {
    if (n1 > n2) { return n2;}
    else { return n1; }
}

int maxim(int n1, int n2) {
    if (n1 < n2) { return n2;}
    else { return n1; }
}

//////////////
//// MAIN ////
//////////////



qvdouble x, y;
qvdouble arry1;
qvdouble arry2;
qvdouble arry3;
qvdouble arry4;



void sedus::dowork() {
    // Reading command line arguments
    //bool fertility[2 * N][PROMETHEUS]
    /*
    using namespace std;
    //clock_t begin = clock();
    cout << "hola" << "\n";
    log->append("SUPERTIME \n");
    letter = argv[1];
    float argR; sscanf (argv[2],"%f",&argR); R = argR;
    rho = R / (4 * N * BLOCKLENGTH);
    float argC; sscanf (argv[3],"%f",&argC); C = argC;
    kappa = C / (4 * N * meanTractLength);
    int sup; sscanf (argv[4],"%d",&sup);
    #undef SUPERTIME
    #define SUPERTIME sup
    float timeToFixationFloat; sscanf (argv[5],"%f",&timeToFixationFloat);
    cout << timeToFixationFloat << "\n";
    int timeToFixation = (int) (timeToFixationFloat * N); // timeToFixation must be >= 1
    cout << timeToFixation << "\n";
    float argDonorRatio; sscanf (argv[6], "%f", &argDonorRatio);
    donorRatio = argDonorRatio;


    correctArguments = 0;
    if (string(argv[7]) == "SC"){
        if (argc == 8){
            crossoverBegin[0] = BLOCKLENGTH;
            crossoverEnd[0] = 2*BLOCKLENGTH;
            correctArguments = 1;
        }
    }
    else if (string(argv[7]) == "WR"){
        if (argc == 8){
            crossoverBegin[0] = 1;
            crossoverEnd[0] = B*BLOCKLENGTH;
            correctArguments = 1;
        }
    }
    else if (string(argv[7]) == "HS"){
        if (argc > 8){
            int num; sscanf (argv[7],"%d",&num);
            numHS = num;
            if (argc == (8+numHS*2)){
                correctArguments = 1;
                for (int HS = 0; HS < numHS; HS++) {
                    int Beg; int End;
                    sscanf (argv[8+HS*2],"%d",&Beg); sscanf (argv[9+HS*2],"%d",&End);
                    crossoverBegin[HS] = Beg;
                    crossoverEnd[HS] = End;
                }
            }
        }
    }
    */
    cout<<"N: "<<N<<endl;
    #undef SUPERTIME
    #define SUPERTIME sup
    perc=0;
    phaseIiter=BURNIN/PROMETHEUS;
    phaseIIiter=((BURNIN+STRUCTURED)/PROMETHEUS)-phaseIiter;
    phaseIIIiter=(TIMELENGTH/PROMETHEUS)-(phaseIiter+phaseIIiter);
    iterationspersupertime=phaseIiter+phaseIIiter+phaseIIIiter;
    maxiterations=SUPERTIME*(phaseIiter+phaseIIiter+phaseIIIiter);
    iterperc=static_cast<float>(100)/maxiterations;
    x.resize(4);
    y.resize(4);
    for(int i=0;i<4;i++){
        x[i].resize(iterationspersupertime);
        y[i].resize(iterationspersupertime);
    }
    for(int j=0;j<4;j++){for (int i=0; i<iterationspersupertime; ++i){x[j][i] = (i*PROMETHEUS)/1000;}}

    arry1.resize(SUPERTIME);
    arry2.resize(SUPERTIME);
    arry3.resize(SUPERTIME);
    arry4.resize(SUPERTIME);
    for(int i=0; i<arry1.length();i++){
        arry1[i].resize(iterationspersupertime);
        arry2[i].resize(iterationspersupertime);
        arry3[i].resize(iterationspersupertime);
        arry4[i].resize(iterationspersupertime);
    }




    if (correctArguments == 1){
        int i, j, h, o;// endTime;
        time_t seconds;
        time(&seconds);
        srand((unsigned int) seconds);

        open_files();

        profile << N << " " << PROMETHEUS << " " << SUPERTIME << " " << BLOCKLENGTH << " " << SAMPLE << " ";
        profile << MUTTABLESIZE << " " << B << " " << BIGTIME << " " << BURNIN / N << " " << STRUCTURED / N << " ";
        profile << mu << " " << rho << " " << kappa << " " << meanTractLength << "\n";

        for (j = 0; j < B; j++) {
            mutationsNewFile[j] << "ms " << SAMPLE << " " << SUPERTIME << " -s 5\n" << seconds << "\n";
        }

        //fertility initialization
        for (int i=0 ; i < 2*N ; i++){
                fertility_ini[i][PROMETHEUS-1] = true;
                for (int tt=0 ; tt < PROMETHEUS-1 ; tt++) {fertility_ini[i][tt] = false;}
         }

        // INITIALIZATION SUPERTIME
        for (run = 0; run < SUPERTIME; run++) {
            if(_abort)break;
            //cout << "SUPERTIME = " << run << "\n";
            QString esci=QString::fromStdString("Run "+std::to_string(run+1)+" of "+std::to_string(SUPERTIME));
            setLog(esci);
            //logger->append(esci);
            time(&seconds);
            srand((unsigned int) seconds);

            ////////////////////////////////////
            //  BUILD THE INITIAL POPULATION  //
            ////////////////////////////////////
            for (i = 0; i < 4 * N; i++) {
                table[i].b = 2; // Each chromosome begins with 2 blocks (b)
                for (j = 0; j < B; j++) {
                    table[i].mpb[j] = 0; // For each block Mutations per block = 0
                }
            }
            for (h = 0; h < 2; h++) {
                for (i = 0; i < 2 * N; i++) {
                    pointer[h][i] = &table[i + 2 * N * h]; // Inside pointer <- table values
                }
            }

            MutCount = 0;
            duFreq = false;
            for (j = 0; j < BLOCKLENGTH; j++) { multihit[j] = false;}

            //////////////////////
            //////// RUN /////////
            //////////////////////

            /*  PHASE I: BURN-IN  */
            cout << "PHASE I" << endl;
            setLog(QString::fromStdString("\tPHASE I"));
            prev_pres ret = phaseI();
            if(_abort)break;
            /* END PHASE I */
            //log->append("AA5\n");
            /*  PHASE II: STRUCTURED TRAJECTORY  */
            cout << "PHASE II" << endl;
            setLog(QString::fromStdString("\tPHASE II"));
            //endTime=phaseII(timeToFixation,ret.prev,ret.pres);
            phaseII(timeToFixation,ret.prev,ret.pres);
            if(_abort)break;
            /* END PHASE II */

            /*  PHASE III: RESOLUTION  */
            cout << "PHASE III" << endl;
            setLog(QString::fromStdString("\tPHASE III"));
            phaseIII();
            if(_abort)break;
            /* END PHASE III */

            for (j = 0; j < B+2; j++) {
                for (o = 0; o < numofsamples; o++) {samplefile[j][o][0] << "\n";samplefile[j][o][1] << "\n";}
            }



            qRegisterMetaType<qvdouble>("qvdouble");
            setChart(x,y);
            /*qvdouble x(250), y(250);
            for (int i=0; i<250; ++i)
              {
                x[i] = i;
                y[i] = run+exp(-i/150.0)*cos(i/10.0); // exponentially decaying cosine
                //y1[i] = exp(-i/150.0);             // exponential envelope
              }*/


         //   auxx << endTime <<"\n";

        } // END SUPERTIME
        close_files();
    }
    else { // INCORRECT ARGUMENTS
            cout << "INCORRECT ARGUMENTS\nUsage: SimulationID R C F/NF SC/WR/HS\nF/NF: Fixed or not fixed duplication fixation trajectory\nIn case of HS you should include the number of hostspots and each crossover begin point and crossover end point as arguments\n"<<endl;
    }
  // clock_t end = clock();
  //  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
   // cout << "time: " << elapsed_secs;
    setBar(100);
    setLog("Finished");
    mutex.lock();
   _working = false;
    mutex.unlock();
    emit finished();
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
//////  FUNCTIONS   /////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

struct sedus::prev_pres sedus::phaseI(){
    // 0 -> BURNIN/PROMETHEUS
    int prev = 0;
    int pres = 1;
    bool does_print = false;
    prev_pres ret;
    for (era = 0; era < (int) BURNIN / PROMETHEUS; era++) {
        if(_abort)break;
        //perc+=iterperc;
        //probar->setValue(int(perc+=iterperc));
        setBar(int(perc+=iterperc));
        genealogy(rho * BLOCKLENGTH, 0);// GENEALOGY (filling recombimatrix, ancestry and fertility matrices, determines all population genealogy)
        int prom=-1;
        prev = 0;
        pres = 1;
        for (std::vector<fertility_info>::iterator it = fertility_list.end()-1 ; it != fertility_list.begin()-1; --it){
            int i = (*it).x;
            int t = (*it).y;
            if(prom!=t){
                if (prev == 1) {prev = 0;pres = 1;} else {prev = 1;pres = 0;}
            }
            prom = t;
            parentpicking(crossoverBegin, crossoverEnd, numHS, prev,pres,i,t);
            mutation(mu, i, pres);
        }
        // CALCULATE THE STATISTICS
        statistics(pres, does_print, era);
    }
    ret.prev = prev;
    ret.pres = pres;
    return ret;
}

int sedus::phaseII(int timeToFixation,int prev, int pres){
    bool does_print = false;
    // Generating Fixation Trajectory in a diploid population
                int endTime = GenerateFixationTrajectory(STRUCTURED + 1, timeToFixation);
                cout << endTime << endl;
                // Picking the first chromosome with the duplication (eva)
                int eva = (int) (rand() % (2 * N));
                for (int i = 0; i < 2 * N; i++) {
                    duplicontent[0][i] = i;
                    duplicontent[1][i] = i;
                }
                duplicontent[1][0] = eva;
                duplicontent[1][eva] = 0;
                // Duplicate eva (create a new block with old mutations...)
                duplication(eva,pres);
                duFreq = true;

                // BURNIN/PROMETHEUS -> (BURNIN+STRUCTURED)/PROMETHEUS (30 -> 50)
                for (era = (int) BURNIN / PROMETHEUS; era < (int) (BURNIN + STRUCTURED) / PROMETHEUS; era++) {
                    if(_abort)break;
                    //perc+=iterperc;
                    //probar->setValue(int(perc+=iterperc));
                    setBar(int(perc+=iterperc));
                    // GENEALOGY (with recombination and taking into account that duplicated chr have duplicated ancestor)
                    genealogy(rho * BLOCKLENGTH, 1);
                    int prom=-1;
                    prev = 0;
                    pres = 1;
                    for (std::vector<fertility_info>::iterator it = fertility_list.end()-1 ; it != fertility_list.begin()-1; --it){
                            int i = (*it).x;
                            int t = (*it).y;
                            if(prom!=t){
                                    if (prev == 1) {prev = 0;pres = 1;} else {prev = 1;pres = 0;}
                            }
                            prom = t;
                            parentpicking(crossoverBegin, crossoverEnd, numHS,prev,pres,i,t);// PARENT PICKING (with recombination)
                            mutation(mu, i,pres);// MUTATION and CONVERSION (for each fertile chromosome)
                            conversion(kappa, i, pres, donorRatio);
                    }
                    // CALCULATE THE STATISTICS
                    statistics(pres, does_print,era);
                }
       return endTime;
}

void sedus::phaseIII(){
    double timefortotal=0.0;
    double timestatstotal=0.0;
    bool does_print = false;
    for (era = (int) (BURNIN + STRUCTURED) / PROMETHEUS; era < (int) TIMELENGTH / PROMETHEUS; era++) {
                    if(_abort)break;
                    //perc+=iterperc;
                    //probar->setValue(int(perc+=iterperc));
                    setBar(int(perc+=iterperc));
                    // GENEALOGY (all chr have the duplication, there are no two populations to take into account)
                    genealogy(rho * BLOCKLENGTH, 0);
                    int prom=-1;
                    int prev = 0;
                    int pres = 1;
                    clock_t tStart1 = clock();
                    for (std::vector<fertility_info>::iterator it = fertility_list.end()-1 ; it != fertility_list.begin()-1; --it){
                        int i = (*it).x;
                        int t = (*it).y;
                        if(prom!=t){
                            if (prev == 1) {prev = 0;pres = 1;} else {prev = 1;pres = 0;}
                        }
                        prom = t;
                        parentpicking(crossoverBegin, crossoverEnd, numHS, prev,pres,i,t);
                        mutation(mu, i,pres);
                        conversion(kappa, i,pres, donorRatio);
                    }
                    double timefor= (double)(clock() - tStart1)/CLOCKS_PER_SEC;
                    timefortotal=timefortotal+timefor;
                    // CALCULATE THE STATISTICS
                    if(era == ((int) (TIMELENGTH/PROMETHEUS)-1)){
                        does_print=true;
                    }
                    clock_t tStart2 = clock();
                    statistics(pres, does_print,era);
                    double timestats= (double)(clock() - tStart2)/CLOCKS_PER_SEC;
                    timestatstotal=timestatstotal+timestats;
    }
    printf("Time taken: %.2fs, %.2fs\n", timefortotal,timestatstotal);
}

//////////////////////////////
////  OPEN & CLOSE FILES  ////
//////////////////////////////

void open_files() { // Opens write-on files
    int j, o;
    stringstream ss;

    ss << "auxx_" << letter << ".dat" << endl;
    ss >> str;
    auxx.open(dir+str.c_str());
    ss.str("");

    ss << "profile_" << letter << ".dat" << endl;
    ss >> str;
    profile.open(dir+str.c_str());
    ss.str("");

    for (j = 0; j < B; j++) {
        ss << "mutations_new_" << j << "_" << letter << ".dat" << endl;
        ss >> str;
        mutationsNewFile[j].open(dir+str.c_str());
        ss.str("");
    }

    for (j = 0; j < B+2; j++) {
        ss << "SFS_" << j << "_" << letter << ".dat" << endl;
        ss >> str;
        SFS[j].open(dir+str.c_str());
        ss.str("");
    }

    for (o = 0; o < numofsamples; o++) {
        for (j = 0; j < B + 2; j++) {
            ss << "samplepi" << j << "[" << sampleN[o] << "]_" << letter << ".dat" << endl;
            ss >> str;
            samplefile[j][o][0].open(dir+str.c_str());
            ss.str("");

            ss << "sampleS" << j << "[" << sampleN[o] << "]_" << letter << ".dat" << endl;
            ss >> str;
            samplefile[j][o][1].open(dir+str.c_str());
            ss.str("");
        }
    }
}

void close_files() { // Closes write-on files
    int j, o;

    profile.close();
    auxx.close();
    for (j = 0; j < B; j++) { mutationsNewFile[j].close();}
    for (j = 0; j < B + 2; j++) {
        SFS[j].close();
        for (o = 0; o < numofsamples; o++) {
            samplefile[j][o][0].close();
            samplefile[j][o][1].close();
        }
    }
}

////////////////////////////////////////////
////  GENEALOGY & STRUCTURED GENEALOGY  ////
////////////////////////////////////////////

void genealogy(float probability, int strornot) { // Generates the genealogy based on the FixationTrajectory (with RECOMBINATION)
    // Determine recombination processes (with probability "probability")
    for (int tt=0 ; tt < PROMETHEUS ; tt++){
        for (int i=0 ; i < 2*N ; i++){
            recombimatrix[i][tt] = false;
        }
    }

    // STRUCTURED GENEALOGY (the population is subdivided in 2: one that carries the duplication, built according to trajectime, and the other not carrying the duplication
    if (strornot == 1){
        // Pick a parent from the corresponding duplicated/non-duplicated population
        int present=0;
        int previous=1;
        for (int tt=0 ; tt < PROMETHEUS ; tt++){
            // trajectime is the time in which we are in fixationTrajectory[] array (used to know the number of chr that carry the dup at each time)
            int trajectime = PROMETHEUS*(era-((int)BURNIN/PROMETHEUS))+tt+1;
            // Randomly mix all the chromosomes of the present generation
            for(int i=0 ; i < 2*N ; i++){
                 int val = (int) (rand()%(2*N));
                 int temp = duplicontent[present][i];
                 // Duplicontent is an array that indicates which chr carry the duplication (with present and previous lines)
                 duplicontent[present][i] = duplicontent[present][val];
                 duplicontent[present][val] = temp;
            }
            // For all the chr that have to have the duplication at this time...
            for(int i=0 ; i < fixationTrajectory[trajectime] ; i++){
                // Choose a father randomly (from the previous duplicated population)
                int val=(int) (rand()%(fixationTrajectory[trajectime-1]));
                ancestry[duplicontent[present][i]][tt] = duplicontent[previous][val];
            }
            // For all the chr that have not to have the duplication at this time...
            for(int i=fixationTrajectory[trajectime] ; i < 2*N ; i++){
                // Choose a father randomly (from the previous non-duplicated population)
                int val = (int) (rand()%(2*N-fixationTrajectory[trajectime-1])) + fixationTrajectory[trajectime-1];
                ancestry[duplicontent[present][i]][tt] = duplicontent[previous][val];
            }
            if (previous==1) {previous=0;  present=1;} else {previous=1;  present=0;}
        }
    }
    // NORMAL GENEALOGY (no 2 populations Duplicated/NoDuplicated)
    else if (strornot == 0){
        for (int tt=0 ; tt < PROMETHEUS ; tt++){
            for (int i=0 ; i < 2*N ; i++){
                int val=(int) (rand()%(2*N));
                ancestry[i][tt] = val; //   WITHOUT presLOSS OF GENERALITY AN INDIV. CAN MATE WITH HIMSELF
            }
        }
    }

    // TRACING BACK THE GENEALOGY AND BUILDING FERTILITY MATRIX
   // std::copy(std::begin(fertility_ini), std::end(fertility_ini), std::begin(fertility));
    //fertility initialization

    //memset(fertility, 0, sizeof(fertility[0][0]) * 2*N * PROMETHEUS);
   // for (int tt=0 ; tt < PROMETHEUS ; tt++) { for (int i=0 ; i < 2*N ; i++){
                    //	 cout << fertility[i][tt] << " ";} cout << "\n";} cout << "\n\n";

    //memcpy(fertility, fertility_ini, sizeof(fertility[0][0]) * 2*N * PROMETHEUS);
    //fertility.swap(fertility_ini);
    fertility = fertility_ini;
    //cout << fertility.size()<<endl;
    //for (int i=0 ; i < 2*N ; i++){
      //                  fertility[i][PROMETHEUS-1] = true;
     // }
                        //for (int tt=0 ; tt < PROMETHEUS ; tt++) { for (int i=0 ; i < 2*N ; i++){
                            //      	    	 cout << fertility[i][tt] << " ";} cout << "\n";} cout << "\n\n";



    fertility_list.clear();
    for (int i=0 ; i < 2*N ; i++){
                fertility_info fi;
                fi.x = i;
                fi.y = (PROMETHEUS-1);
                fertility_list.push_back(fi);
    }
   for (int tt=PROMETHEUS-1 ; tt > 0 ; tt--){
        for (int i=0 ; i < 2*N ; i++){
            if (fertility[i][tt] == true){

                if(fertility[ancestry[i][tt]][tt-1]==false){
                    fertility_info fi;
                    fi.x=(ancestry[i][tt]);
                    fi.y=(tt-1);
                    fertility_list.push_back(fi);
                }
                fertility[ancestry[i][tt]][tt-1] = true;
                float p = rand() / ((float) RAND_MAX + 1);// Determine recombination processes (with probability "probability")
                if (p < probability){

                    recombimatrix[i][tt] = true;
                    if(ancestry[i][tt]%2==0) {
                        if(fertility[ancestry[i][tt]+1][tt-1]==false){
                            fertility_info fi;
                            fi.x=(ancestry[i][tt]+1);
                            fi.y=(tt-1);
                            fertility_list.push_back(fi);
                        }
                        fertility[ancestry[i][tt]+1][tt-1]=true;
                    }
                    else {
                        if(fertility[ancestry[i][tt]-1][tt-1]==false){
                            fertility_info fi;
                            fi.x=(ancestry[i][tt]-1);
                            fi.y=(tt-1);
                            fertility_list.push_back(fi);
                        }
                        fertility[ancestry[i][tt]-1][tt-1]=true;
                    }
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////
////  PARENTPICKING, DUPLICATION, MUTATION, CONVERSION  ////
////////////////////////////////////////////////////////////

void parentpicking(int crossBegin[maxNumOfHS], int crossEnd[maxNumOfHS], int numCrossRegions, int prev, int pres, int i, int t) {
    int j, k, junctionBlock, defHS, HS;
    int father, partner, junction, vald0, valr0, childblocks, minblock, recTract, end, beg;
    chrom * chr;
    float p;
    double num;

    // PARENT-PICKING
            father = ancestry[i][t];
            // RECOMBINATION
            if (recombimatrix[i][t] == true) {
                childblocks = pointer[prev][father]->b;
                if (ancestry[i][t] % 2 == 0) {
                    partner = ancestry[i][t] + 1;
                } else {
                    partner = ancestry[i][t] - 1;
                }
                chr = pointer[pres][i];
                chr->b = childblocks;
                minblock = minim(pointer[prev][father]->b, pointer[prev][partner]->b);
                //maxblock = maxim(pointer[prev][father][0].b, pointer[prev][partner][0].b);
                p = rand() / ((float) RAND_MAX + 1);
                defHS = 0;
                for (HS = numCrossRegions; HS > 0; HS--) {
                    num = HS;
                    num /= numCrossRegions;
                    if (p < num){
                        defHS = HS-1;
                    }
                }
                end = crossEnd[defHS];
                beg = crossBegin[defHS];
                if (((minblock * BLOCKLENGTH) < end) && ((minblock * BLOCKLENGTH) > beg)) {
                    end = minblock * BLOCKLENGTH;
                }
                if (((minblock * BLOCKLENGTH) >= end) && ((minblock * BLOCKLENGTH) > beg)) {//LA SEGUNDA CONDICION ESTA IMPLCITA EN LA PRIMERA
                    // If neither father and partner have mutations in any block
                    if ((pointer[prev][father]->mpb[0] == 0) && (pointer[prev][father]->mpb[1] == 0) &&
                        (pointer[prev][father]->mpb[2] == 0) && (pointer[prev][partner]->mpb[0] == 0) &&
                        (pointer[prev][partner]->mpb[1] == 0) && (pointer[prev][partner]->mpb[2] == 0)) {
                        for (j = 0; j < chr->b; j++) {
                            chr->mpb[j] = 0;
                        }
                    }// If at least one of them have one or more mutations in one or more blocks
                    else {

                        recTract = end - beg;
                        junction = (int) (rand() % (recTract));
                        junction += beg;
                        junctionBlock = (int) (junction / BLOCKLENGTH); // block where junction fell
                        //       cout << junctionBlock << " ";
                        valr0 = location(junction - junctionBlock*BLOCKLENGTH, prev, partner, junctionBlock); // junction location in the block "junctionBlock" of partner
                        vald0 = location(junction - junctionBlock*BLOCKLENGTH, prev, father, junctionBlock); // junction location in the block "junctionBlock" of father

                        //COPY THE INFO FROM PARTNER (BLOCK 0)
                        for (j = 0; j < junctionBlock; j++) {
                            chr->mpb[j] = pointer[prev][partner]->mpb[j];
                            for (k = 0; k < pointer[prev][partner]->mpb[j]; k++) {
                                chr->mutation[j][k] = pointer[prev][partner]->mutation[j][k];
                            }
                        }

                        //COPY INFO IN JUNCTION BLOCK (BLOCK 1)
                        //FROM PARTNER
                        for (k = 0; k < valr0; k++) {
                            chr->mutation[junctionBlock][k] = pointer[prev][partner]->mutation[junctionBlock][k];
                        }
                        //FROM FATHER
                        for (k = 0; k < (pointer[prev][father]->mpb[junctionBlock] - vald0); k++) {
                            chr->mutation[junctionBlock][k + valr0] = pointer[prev][father]->mutation[junctionBlock][k + vald0];
                        }
                        //ESTABLISH MPB
                        chr->mpb[junctionBlock] = valr0 + (pointer[prev][father]->mpb[junctionBlock] - vald0);

                        //COPY INFO FROM BLOCK 2 IF PRESENT
                        for (j = junctionBlock + 1; j < chr->b; j++) {
                            if (j < pointer[prev][father]->b) { //j counts 0,1,2 but .b counts 1,2,3
                                chr->mpb[j] = pointer[prev][father]->mpb[j];
                                for (k = 0 ; k < pointer[prev][father]->mpb[j] ; k++) {
                                    chr->mutation[j][k] = pointer[prev][father]->mutation[j][k];
                                }
                            } else {
                                chr->mpb[j] = 0;
                            }
                        }
                    }
                } else if (((minblock * BLOCKLENGTH) < end) && ((minblock * BLOCKLENGTH) <= beg)) {//LA PRIMERA CONDICION ESTA IMPLICITA EN SEGUNDA
                     copychr(prev, father, pres, i);
                }
            }// NO RECOMBINATION
            else {
                copychr(prev, father, pres, i);
            }

}

void duplication(int i,int prev) {
    int k;
    int tempMutCount = 0;

    if (pointer[prev][i][0].b == 2) {
        pointer[prev][i][0].b++;
    }
    pointer[prev][i][0].mpb[2] = pointer[prev][i][0].mpb[0];
    for (k = 0; k < pointer[prev][i][0].mpb[0]; k++) {
        pointer[prev][i][0].mutation[2][k] = pointer[prev][i][0].mutation[0][k];
    }
    for (k = 0; k < MutCount; k++) {
        if (muttable[k].block == 0) {
            muttable[MutCount + tempMutCount].position = muttable[k].position;
            muttable[MutCount + tempMutCount].block = 2;
            tempMutCount++;
        }
    }
    MutCount += tempMutCount;
}


void mutation(float probability, int i, int pres) {
    int j, k, position, val, mutEvents;
    float p;
    chrom * chr;
    chr = pointer[pres][i];

    for (j = 0; j < chr->b; j++) {
        mutEvents=0;
        p = rand() / ((float) RAND_MAX + 1);
        if (p < (probability * BLOCKLENGTH)) { mutEvents++;
            if (p < pow((probability * BLOCKLENGTH), 2)) { mutEvents++;
                if (p < pow((probability * BLOCKLENGTH), 3)) { mutEvents++;
                    if (p < pow((probability * BLOCKLENGTH),4)) { mutEvents++;
                    }
                }
            }
        }
        for(int muts=0; muts < mutEvents; muts++){
            // Randomly choose one NEW mutation position
            do {
                position = (int) (rand() % (BLOCKLENGTH));
            } while (multihit[position] == true);
            multihit[position] = true;

            // If there are no previous mutations or if the new mutation fall behind all previous ones...,
            if (chr->mpb[j] == 0 || position > chr->mutation[j][chr->mpb[j] - 1]) {
                // Simply add the mutation
                chr->mutation[j][chr->mpb[j]] = position;
                chr->mpb[j]++;
            }// Otherwise, shift the mutations one position so that mutations appear ordered in the chr.mutation array
            else {
                val = location(position, pres, i, j);
                for (k = chr->mpb[j]; k > val; k--) {
                    chr->mutation[j][k] = chr->mutation[j][k - 1];
                }
                chr->mutation[j][val] = position;
                chr->mpb[j]++;
            }
            // ELABORATES MUTTABLE
            muttable[MutCount].position = position;
            muttable[MutCount].block = j;
            MutCount++;

            // IF MUTATION APPEARS IN BLOCK 0 OR 2, COPY THE MUTATION TO THE OTHER BLOCK
            if (j == 0 && duFreq == true) {
                muttable[MutCount].block = 2;
                muttable[MutCount].position = muttable[MutCount - 1].position;
                MutCount++;
            }
            if (j == 2) {
                muttable[MutCount].block = 0;
                muttable[MutCount].position = muttable[MutCount - 1].position;
                MutCount++;
            }
        }
    }

}

void conversion(float probability, int i, int pres, float donorRatio) {
    int k, donor, receptor, vald0, valdf, valr0, valrf, val, junction, tractlength, eqmut, otherk, differences;
    float p;
    chrom * chr;
    chr = pointer[pres][i];
    // If the chr has duplicated block
    if (chr[0].b == 3) {
        p = rand() / ((float) RAND_MAX + 1);

        // GENE CONVERSION
        if (p < (2 * probability * BLOCKLENGTH)) {
            // Determines which block will be the donor and which will be the receptor
            p = rand() / ((float) RAND_MAX + 1);
            if (p < donorRatio) {
                donor = 0;
                receptor = 2;
            } else {
                donor = 2;
                receptor = 0;
            }
            // Determines conversion initiation point
            junction = (int) (rand() % (BLOCKLENGTH));
            // Determines gene conversion tract length through function tractpql()
            tractlength = tractpql(meanTractLength);
            // If tractlength is an odd number, we must correct junction in order for gene conversion to be balanced
            if(tractlength%2 != 0){
                p = rand() / ((float) RAND_MAX + 1);
                if(p < 0.5){
                    junction += 1;
                }
            }

            // Tests for MEPS, 100% identity tract to allow conversion
            vald0 = location(junction - numbptotalidentity/2, pres, i, donor);
            valdf = location(junction + numbptotalidentity/2, pres, i, donor);
            valr0 = location(junction - numbptotalidentity/2, pres, i, receptor);
            valrf = location(junction + numbptotalidentity/2, pres, i, receptor);
            val = (valdf - vald0)-(valrf - valr0); // val is the number of positions to shift
            eqmut = 0;
            // If the sequences have the same number of mutations in the total identity region
            if (val == 0){
                for (k = 0; k < (valdf - vald0); k++) {
                    if (chr[0].mutation[receptor][k + valr0] == chr[0].mutation[donor][k + vald0]){
                        eqmut++;
                    }
                }
                // If the sequences have the same mutations in the total identity region
                if (eqmut == (valdf - vald0)){

                    // Determines location of the delimiting mutations for each block (donor and receptor)
                    vald0 = location(junction - 0.5 * tractlength, pres, i, donor);
                    valdf = location(junction + 0.5 * tractlength, pres, i, donor);
                    valr0 = location(junction - 0.5 * tractlength, pres, i, receptor);
                    valrf = location(junction + 0.5 * tractlength, pres, i, receptor);
                    val = (valdf - vald0)-(valrf - valr0); // val is the number of positions to shift
                    // Calculates differences between both sequences
                    differences = 0;
                    //Corrected bug; it used to say k<=
                    for (k = vald0; k < valdf; k++) {
                        otherk = location(chr[0].mutation[donor][k], pres, i, receptor);
                        if (chr[0].mutation[receptor][otherk] != chr[0].mutation[donor][k]){
                            differences++;
                        }
                    }
                    //Corrected bug; it used to say k<=
                    for (k = valr0; k < valrf; k++) {
                        otherk = location(chr[0].mutation[receptor][k], pres, i, donor);
                        if (chr[0].mutation[donor][otherk] != chr[0].mutation[receptor][k]){
                            differences++;
                        }
                    }
                    // If differences do not exceed the identity threshold for all the tractlength
                    if(differences <= tractlength*(1-(IdentityInConvTract/100))){

                        // Change receptor mutations in chr[0].mutation[receptor][] array
                        if (val > 0) {
                            for (k = chr[0].mpb[receptor] + val - 1; k >= valrf + val; k--) {
                                chr[0].mutation[receptor][k] = chr[0].mutation[receptor][k - val];
                            }
                        }
                        if (val < 0) {
                            for (k = valrf + val; k < chr[0].mpb[receptor] + val; k++) {
                                chr[0].mutation[receptor][k] = chr[0].mutation[receptor][k - val];
                            }
                        }
                        for (k = 0; k < (valdf - vald0); k++) {
                            chr[0].mutation[receptor][k + valr0] = chr[0].mutation[donor][k + vald0];
                        }
                        // Change receptor chr[0].mpb[] array
                        chr[0].mpb[receptor] += val;
                    }
                }
            }
        }
    }
}

//////////////////////
////  STATISTICS  ////
//////////////////////

void statistics(int prev, bool does_print, int era) {
    int j, o;
    float * resultsSample;

    FSL(prev); // Creating the summary vector fixedLostForAll

    // INFORMATION RECOVERED FROM SAMPLES
    for (o = 0; o < numofsamples; o++) {
        SamplingIndividuals(sampleN[o]);

        // CALCULATES THE SFS FOR THE SAMPLE BLOCK BY BLOCK
        for (j = 0; j < B; j++) {
            resultsSample = SiteFrequencySpectrumPrint(prev, j, sampleN[o], does_print);
            samplefile[j][o][0] << resultsSample[1] << " "; // the results array keeps (S,pi) but for the sake of tradition
            samplefile[j][o][1] << resultsSample[0] << " "; // we save it as (pi,S)
            double sum=0;
            switch (j)
            {
                case 0: {arry1[run][era]=resultsSample[1];for(int i=0;i<=run;i++){sum+=arry1[i][era];};y[0][era]=sum/(run+1);break;}
                case 1: {arry2[run][era]=resultsSample[1];for(int i=0;i<=run;i++){sum+=arry2[i][era];};y[1][era]=sum/(run+1);break;}
                case 2: {arry3[run][era]=resultsSample[1];for(int i=0;i<=run;i++){sum+=arry3[i][era];};y[2][era]=sum/(run+1);break;}
            }
        }

        resultsSample = SiteFrequencySpectrum_02(prev, sampleN[o], does_print); // Collapsed for samples
        samplefile[B][o][0] << resultsSample[1] << " "; // the results array keeps (S,pi) but for the sake of tradition
        samplefile[B][o][1] << resultsSample[0] << " "; // we save it as (pi,S)

        arry4[run][era]=resultsSample[1];
        double sum=0;
        for(int i=0;i<=run;i++){
            sum+=arry4[i][era];
        }
        y[3][era]=sum/(run+1);

        resultsSample = SiteFrequencySpectrum_02_Calling(prev, sampleN[o], does_print); // Collapsed for samples with calling
        samplefile[B+1][o][0] << resultsSample[1] << " "; // the results array keeps (S,pi) but for the sake of tradition
        samplefile[B+1][o][1] << resultsSample[0] << " "; // we save it as (pi,S)

    }

    // CALCULATES THE NUMBER OF PRIVATE & SHARED MUTATIONES BETWEEN BLOCKS 0 & 2
    DivergenceForAll(0, 2, prev);

}

/////////////////////////////////////////////////////////////////////
////  FSL (FIXED-SEGREGATING-LOST): EXTRACTS INFO FROM MUTTABLE  ////
/////////////////////////////////////////////////////////////////////

void FSL(int hh) {
    int m, otherm=0, mm;

    for (m = 0; m < MutCount; m++) {
        muttable[m].frequency = muFrequencyIntWholePopAndSample(hh, muttable[m].position, muttable[m].block, N);
        muttable[m].frequency = ((float) muttable[m].frequency) / (N * 2);
    }
    for (m = 0, mm = 0; m < MutCount; m++) {
        // SET MULTIHIT TO FALSE FOR MUTATIONS THAT HAVE BEEN LOST
        if (muttable[m].frequency == 0 && muttable[m].block == 1) {
            multihit[muttable[m].position] = false;
        }
        if (muttable[m].frequency == 0 && muttable[m].block == 0 && duFreq == false) {
            multihit[muttable[m].position] = false;
        }
        // IF MUTATION IS IN THE SINGLE-COPY BLOCK
        if (muttable[m].block == 1) {
            // SEGREGATING
            if (muttable[m].frequency > 0 && muttable[m].frequency < 1) {
                temporalmuttable[mm].block = muttable[m].block;
                temporalmuttable[mm].frequency = muttable[m].frequency;
                temporalmuttable[mm].position = muttable[m].position;
                mm++;
            }
            // FIXED
            if (muttable[m].frequency == 1) {
                EraseFixedMutations(muttable[m].position, muttable[m].block, hh);
            }

        }// IF THE MUTATION IS IN EITHER BLOCK 0 0R 2
        else {
            // IF THE DUPLICATION HAS NOT YET OCCURRED, TREATS BLOCK 0 AS SINGLE-COPY
            if (duFreq == false) {
                // SEGREGATING
                if (muttable[m].frequency > 0 && muttable[m].frequency < 1) {
                    temporalmuttable[mm].block = muttable[m].block;
                    temporalmuttable[mm].frequency = muttable[m].frequency;
                    temporalmuttable[mm].position = muttable[m].position;
                    mm++;
                }
                // FIXED
                if (muttable[m].frequency == 1) {
                    EraseFixedMutations(muttable[m].position, muttable[m].block, hh);
                }
            }// IF THE DUPLICATION HAS OCCURRED
            else {
                // CHECKS THE BLOCK IN WHICH IT HAS OCCURRED
                // AND LOOKS FOR THE POSITION OF THE MUTATION IN MUTTABLE FOR THE OTHER BLOCK
                if (muttable[m].block == 2) {
                    otherm = SearchMutation(0, muttable[m].position, MutCount);
                }
                if (muttable[m].block == 0) {
                    otherm = SearchMutation(2, muttable[m].position, MutCount);
                }
                // SEGREGATING
                if (muttable[m].frequency > 0 && muttable[m].frequency < 1) {
                    temporalmuttable[mm].block = muttable[m].block;
                    temporalmuttable[mm].frequency = muttable[m].frequency;
                    temporalmuttable[mm].position = muttable[m].position;
                    mm++;
                    // IN CASE THE MUTATION IS NOT PRESENT IN THE OTHER BLOCK, IT IS COPIED TO MUTTABLE ANYWAY
                    if (muttable[otherm].frequency == 0) {
                        temporalmuttable[mm].block = muttable[otherm].block;
                        temporalmuttable[mm].frequency = muttable[otherm].frequency;
                        temporalmuttable[mm].position = muttable[otherm].position;
                        mm++;
                    }
                }
                // FIXED
                if (muttable[m].frequency == 1) {
                    // IF THE DUPLICATION IS FIXED IN THE OTHER BLOCK
                    if (muttable[otherm].frequency == 1) {
                        EraseFixedMutations(muttable[m].position, muttable[m].block, hh);
                    }// IF THE DUPLICATION IS SEGREGATING OR NOT PRESENT IN THE OTHER BLOCK
                    else {
                        temporalmuttable[mm].block = muttable[m].block;
                        temporalmuttable[mm].frequency = muttable[m].frequency;
                        temporalmuttable[mm].position = muttable[m].position;
                        mm++;
                        // IN CASE THE MUTATION IS NOT PRESENT IN THE OTHER BLOCK, IT IS COPIED TO MUTTABLE ANYWAY
                        if (muttable[otherm].frequency == 0) {
                            temporalmuttable[mm].block = muttable[otherm].block;
                            temporalmuttable[mm].frequency = muttable[otherm].frequency;
                            temporalmuttable[mm].position = muttable[otherm].position;
                            mm++;
                        }
                    }
                }
                // LOST (TO CONSIDER IT LOST IT WOULD HAVE TO BE LOST IN BOTH BLOCKS)
                if (muttable[m].frequency == 0) {
                    if (muttable[otherm].frequency == 0) {
                        multihit[muttable[m].position] = false;
                    }
                }
            }
        }
    }
    MutCount = mm;
    for (m = 0; m < MutCount; m++) {
        muttable[m].block = temporalmuttable[m].block;
        muttable[m].frequency = temporalmuttable[m].frequency;
        muttable[m].position = temporalmuttable[m].position;
    }
}

/////////////////////////////////////////////////////
////  COPYCHR, LOCATION & ERASE FIXED MUTATIONS  ////
/////////////////////////////////////////////////////

// COPY A CHROMOSOME (WHEN THERE IS NO RECOMBINATION)
void copychr(int prev, int ind0, int pres, int ind1) { // (origin,end)
    chrom *c1;
    int j, k;
    c1 = pointer[pres][ind1];
    c1[0].b = pointer[prev][ind0][0].b;
    for (j = 0; j < pointer[prev][ind0][0].b; j++) {
        c1[0].mpb[j] = pointer[prev][ind0][0].mpb[j];
        for (k = 0; k < pointer[prev][ind0][0].mpb[j]; k++) {
            c1[0].mutation[j][k] = pointer[prev][ind0][0].mutation[j][k];
        }
    }
}

// LOCATE A GIVEN MUTATION OR POSITION INSIDE A mutation VECTOR OF A GIVEN CHROMOSOME
int location(int position, int h, int ind, int j) {// Found the position in the c[0].mutation array where "position" should be
    chrom *c;
    int k = 0;
    c = pointer[h][ind];
    bool found = false;
    while (found == false && k <= c->mpb[j]) {
        if (position <= c->mutation[j][k] && k < c->mpb[j]) {
            found = true;
        }
        k++;
    }
    return (k - 1);
}

// ERASES MUTATIONS
void EraseFixedMutations(int fixed, int block, int prev) {
    int ii, k;

    multihit[fixed] = false;
    for (ii = 0; ii < 2 * N; ii++) {
        pointer[prev][ii][0].mpb[block]--;
        for (k = location(fixed, prev, ii, block); k < pointer[prev][ii][0].mpb[block]; k++) {
            pointer[prev][ii][0].mutation[block][k] = pointer[prev][ii][0].mutation[block][k + 1];
        }
    }
}


///////////////////////////////////////////////////////////////
////  SITE FREQUENCY SPECTRUM, BY BINS SFS, COLLAPSED SFS  ////
///////////////////////////////////////////////////////////////

// Calculates Site Frequency Spectrum with Printing option for mutationNew
// also prints SFS for blocks 0, 1 and 2
// also calculates pi, S, eta, H
float * SiteFrequencySpectrumPrint(int h, int block, int n, bool does_print) {
    int s = n;
    int p, i, j, k, index, m, number, duplicationFreq;
    int xi[2 * s], list [2 * N];
    float value, mufreq;
    static float results[] = {0, 0};

    results[0] = 0;
    results[1] = 0;

        if(block == 2){
            //duplicationFreq = DupliFreq(prev, 2, n);
            duplicationFreq = DupliFreq(h, 2, n);
            auxx << duplicationFreq << " ";
            s = (int) duplicationFreq/2;
        }
        if (s == 0) { return results; }

        for (m = 0, number = 0; m < MutCount; m++) {
            if (block == muttable[m].block && muttable[m].frequency != 0) {
                    if (muttable[m].frequency != 1) {
                        mufreq = (float) muFrequencyIntWholePopAndSample(h, muttable[m].position, muttable[m].block,n);
                        mufreq = (mufreq)/(2*s);
                        if(mufreq != 0){
                            list[number] = muttable[m].position;
                            number++; // number = segregatingSites in sample
                        }
                    }
            }
        }
        results[0] = number;
       // number=5;
        if(does_print == true){
              mutationsNewFile[block] << "//\nsegsites: "<< number <<"\npositions: ";
        }

        if (number == 0) { return results; }

        int prototype[number];
        int mutationCounts[number];
        for (index = 0; index < 2 * s; index++) {
            xi[index] = 0;
        }

        // BUILD THE PROTOTYPE: AN ORDERED SEQUENCE WITH ALL POSSIBLE MUTATION POSITIONS
        for (m = 0, p = 0; m < number; m++) {
                if (p == 0) {
                        prototype[p] = list[m];
                        p++;
                }else {
                        j = 0;
                        while (list[m] > prototype[j] && j < p) { j++; }
                        if (j < p) {
                                for (k = p + 1; k > j; k--) {
                                        prototype[k] = prototype[k - 1];
                                }
                        }
                        prototype[j] = list[m];
                        p++;
                }
        }

        // PRINTS POSITIONS OF SEGREGATING SITES WITH FIVE DECIMAL POSITIONS
        if(does_print == true){
                     float rounded_mutation;
                     for(int nn=0; nn < number ; nn++){
                      rounded_mutation = (float) prototype[nn]/BLOCKLENGTH;
                      rounded_mutation = round(rounded_mutation,5);
                        // mutationsNewFile[block] << prototype[nn] << " ";
                      mutationsNewFile[block] << rounded_mutation << " ";
                     }
                     mutationsNewFile[block] << "\n";
               }

        // THIS FUNCTION PRINTS THE CONTENT OF MUTATIONSNEWFILE
        if(does_print == true){
            int pos;
                for (i = 0; i < 2 * s; i++) {
                        for (j = 0; j < number; j++) {
                            pos = prototype[j];
                            k = location(prototype[j], h, sample[i], block);
                            if (pointer[h][sample[i]][0].mutation[block][k] == pos && k < pointer[h][sample[i]][0].mpb[block]) {
                                mutationsNewFile[block] << "1";
                            }else{
                                mutationsNewFile[block] << "0";
                            }
                         }
                        mutationsNewFile[block] << "\n";
                }
        }

        //FIND THE FREQUENCY OF EACH MUTATION
        //MUTCOUNTS IS AN ARRAY THAT COUNTS THE NUMBER OF TIMES EACH MUTATION APPEARS
        for (j = 0; j < number; j++) {
                mutationCounts[j] = muFrequencyIntWholePopAndSample(h, prototype[j], block, n);
        }
        //FIND THE UNFOLDED SITE FREQUENCY SPECTRUM
        //xi[i] IS AN ARRAY THAT COUNTS THE NUMBER OF POLYMORPHIC SITES THAT HAVE i COPIES OF THE MUTATION
        for (i = 1; i < 2 * s; i++) {
                for (j = 0; j < number; j++) {
                        if (mutationCounts[j] == i) {
                                xi[i]++;
                        }
                }
        }

        if(does_print==true){
                     for (i = 1; i < 2 * s; i++) {
                        SFS[block] << xi[i] << " ";
                     }
                     SFS[block] << "\n";
        }

        // CALCULATE AVERAGE PAIRWISE DIFFERENCE
        for (i = 1, value = 0; i < 2 * s; i++) {
                value += (float) (i * (2 * s - i) * xi[i]);
        }
        value /= (s * (2 * s - 1));
        results[1] = value;
        return results;

}



// CALCULATES SFS COLLAPSING BLOCKS 0 & 2 FOR WHOLE POP AND FOR A SAMPLE
// PARAMETERS ARE h and SAMPLE SIZE nn and boolean does_print. IF SAMPLE SIZE IS N EVERYTHING IS TREATEAD
// FOR WHOLE POP, OTHERWISE, IT USES THE SAMPLE
// PERFORMS ALL CALCULATIONS TAKING ONLY THE FIRST HALF OF THE SAMPLE AND TAKING BLOCKS 0 & 2 FROM EACH CHROMOSOME

float * SiteFrequencySpectrum_02(int h, int nn, bool does_print) { // Collapsed
    int n, p, i, j, index, m, number, block, blockCount = nn;
    int mutationCountsTemp;
    int xi[2 * nn], list[2 * N];
    float value;
    static float results[] = {0, 0};

    for (m = 0, number = 0; m < MutCount; m++) {
        if (muttable[m].block == 0) {
            list[number] = muttable[m].position;
            number++; // number = segregatingSites in sample
        }
    }

    results[0] = number;
    results[1] = 0;
    if (number == 0) {
        return results;
    }
    int prototype[number];
    int mutationCounts[number];

    for (index = 0; index < 2 * nn; index++) {
        xi[index] = 0;
    }
    // BUILD THE PROTOTYPE: AN ORDERED SEQUENCE WITH ALL POSSIBLE MUTATION POSITIONS
    for (m = 0, p = 0; m < number; m++) {
        mutationCounts[m] = 0;
        if (p == 0) {
            prototype[p] = list[m];
            p++;
        } else {
            j = 0;
            while (list[m] > prototype[j] && j < p) {
                j++;
            }
            if (j < p) {
                for (n = p + 1; n > j; n--) {
                    prototype[n] = prototype[n - 1];
                }
            }
            prototype[j] = list[m];
            p++;
        }
    }
    //FIND THE FREQUENCY OF EACH MUTATION
    //MUTCOUNTS IS AN ARRAY THAT COUNTS THE NUMBER OF TIMES EACH MUTATION APPEARS
    for (j = 0; j < number; j++) {
            block = 0;
            mutationCountsTemp = muFrequencyIntWholePopAndSample(h, prototype[j], block, nn/2);
            mutationCounts[j] += mutationCountsTemp;
            block = 2;
            mutationCountsTemp = muFrequencyIntWholePopAndSample(h, prototype[j], block, nn/2);
            mutationCounts[j] += mutationCountsTemp;
    }
    for (index = 0; index < nn; index++) {
        if (pointer[h][index][0].b == 3) {
            blockCount++; // number of chr having the duplication in the sample
        }
    }
    //FIND THE UNFOLDED SITE FREQUENCY SPECTRUM
    //xi[i] IS AN ARRAY THAT COUNTS THE NUMBER OF POLYMORPHIC SITES THAT HAVE i COPIES OF THE MUTATION
    for (i = 1; i < blockCount; i++) {
        for (j = 0; j < number; j++) {
            if (mutationCounts[j] == i) {
                xi[i]++;
            }
        }
    }

    if(does_print==true){
         for (i = 1; i < 2 * nn; i++) {
            SFS[3] << xi[i] << " ";
        }
        SFS[3] << "\n";
    }

    // CALCULATE AVERAGE PAIRWISE DIFFERENCE
    for (i = 1, value = 0; i < blockCount; i++) {
        value += (float) (i * (blockCount - i) * xi[i]);
    }
    value /= (blockCount * (blockCount - 1) / 2);
    results[1] = value;
    return results;
}

float * SiteFrequencySpectrum_02_Calling(int h, int nn, bool does_print) { // Collapsed
    int n, p, i, j, index, m, number, blockCount = 2 * nn;
    int mutationCountsTemp;
    int xi[2 * nn], list[2 * N];
    float value;
    static float results[] = {0, 0};

    for (m = 0, number = 0; m < MutCount; m++) {
        if (muttable[m].block == 0) {
            list[number] = muttable[m].position;
            number++; // number = segregatingSites in sample
        }
    }

    results[0] = number;
    results[1] = 0;
    if (number == 0) {
        return results;
    }
    int prototype[number];
    int mutationCounts[number];

    for (index = 0; index < 2 * nn; index++) {
        xi[index] = 0;
    }
    // BUILD THE PROTOTYPE: AN ORDERED SEQUENCE WITH ALL POSSIBLE MUTATION POSITIONS
    for (m = 0, p = 0; m < number; m++) {
        mutationCounts[m] = 0;
        if (p == 0) {
            prototype[p] = list[m];
            p++;
        } else {
            j = 0;
            while (list[m] > prototype[j] && j < p) {
                j++;
            }
            if (j < p) {
                for (n = p + 1; n > j; n--) {
                    prototype[n] = prototype[n - 1];
                }
            }
            prototype[j] = list[m];
            p++;
        }
    }
    //FIND THE FREQUENCY OF EACH MUTATION
    //MUTCOUNTS IS AN ARRAY THAT COUNTS THE NUMBER OF TIMES EACH MUTATION APPEARS

    for (j = 0; j < number; j++) {
            mutationCountsTemp = muFrequencyCollapsedCallingFromSample(h, prototype[j], nn);
            mutationCounts[j] += mutationCountsTemp;
        }
    for (index = 0; index < 2 * nn; index++) {
        if (pointer[h][index][0].b == 3) {
            blockCount++; // number of chr having the duplication in the sample
        }
    }
    //FIND THE UNFOLDED SITE FREQUENCY SPECTRUM
    //xi[i] IS AN ARRAY THAT COUNTS THE NUMBER OF POLYMORPHIC SITES THAT HAVE i COPIES OF THE MUTATION
    for (i = 1; i < blockCount; i++) {
        for (j = 0; j < number; j++) {
            if (mutationCounts[j] == i) {
                xi[i]++;
            }
        }
    }

    if(does_print==true){
         for (i = 1; i < 2*nn; i++) {
            SFS[4] << xi[i] << " ";
        }
        SFS[4] << "\n";
    }

    // CALCULATE AVERAGE PAIRWISE DIFFERENCE
    for (i = 1, value = 0; i < blockCount; i++) {
        value += (float) (i * (blockCount - i) * xi[i]);
    }
    value /= (blockCount * (blockCount - 1) / 2);
    results[1] = value;
    return results;
}


////////////////////////////////////////////////////////////////
////  DIVERGENCE BETWEEN IN SAME AND DIFFERENT CHROMOSOMES  ////
////////////////////////////////////////////////////////////////

// CALCULATES DIVERGENCE BETWEEN TWO BLOCKS IN THE SAME CHROMOSOME
void DivergenceForAll(int blockA, int blockB, int h) {
    int i, k, j;
    int mutBlockA = 0;
    int mutBlockB = 0;
    int mutSharedAB = 0;
    int mutSharedABComp = 0;

    for (i = 0; i < 2 * N; i++) {
        for (k = 0; k < pointer[h][i][0].mpb[blockA]; k++) {
            j = location(pointer[h][i][0].mutation[blockA][k], h, i, blockB);
            if ((pointer[h][i][0].mutation[blockA][k] == pointer[h][i][0].mutation[blockB][j]) && (j < pointer[h][i][0].mpb[blockB])) {
                mutSharedAB++;
            } else {
                mutBlockA++;
            }
        }
        for (k = 0; k < pointer[h][i][0].mpb[blockB]; k++) {
            j = location(pointer[h][i][0].mutation[blockB][k], h, i, blockA);
            if ((pointer[h][i][0].mutation[blockB][k] == pointer[h][i][0].mutation[blockA][j]) && (j < pointer[h][i][0].mpb[blockA])) {
                mutSharedABComp++;
            } else {
                mutBlockB++;
            }
        }
    }
    if (mutSharedAB != mutSharedABComp) {
        cout << "ERROR en divergence for all\n" << mutSharedAB << " " << mutSharedABComp << endl;
    }
}


////////////////////////////
////  SIMPLE FUNCTIONS  ////
////////////////////////////

// CALCULATES ABSOLUTE FREQUENCY OF THE DUPLICATION
int DupliFreq(int h, int block, int n) {
    int i = 0, quantity = 0;

    if (n == N){
        for (i = 0; i < 2 * n; i++) {
                if (pointer[h][i][0].b == (block + 1)) {
                quantity++;
                }
        }
    }else{
        for (i = 0; i < 2 * n; i++) {
                if (pointer[h][sample[i]][0].b == (block + 1)) {
                quantity++;
                }
        }
    }
    return quantity;
}

// THIS FUNCTION SERVES TO SEARCH FOR THE ABSOLUTE FREQUENCY OF A GIVEN MUTATION IN THE WHOLE POPULATION OR A SAMPLE
int muFrequencyIntWholePopAndSample(int h, int pos, int block, int n) {
    int i, quantity, k;
    if (n == N) {
        for (i = 0, quantity = 0; i < 2 * n; i++) {
            k = location(pos, h, i, block);
            if (pointer[h][i][0].mutation[block][k] == pos && k < pointer[h][i][0].mpb[block]) {
                quantity++;
            }
        }
    } else {
        for (i = 0, quantity = 0; i < 2 * n; i++) {
            k = location(pos, h, sample[i], block);
            if (pointer[h][sample[i]][0].mutation[block][k] == pos && k < pointer[h][sample[i]][0].mpb[block]) {
                quantity++;
            }
        }
    }
    return quantity;
}

// THIS FUNCTION DOES THE CALLING FOR COLLAPSED DUPLICATIONS OF AN INDIVIDUAL FROM SAMPLE
int muFrequencyCollapsedCallingFromSample(int h, int pos, int n) {
    int i, j, quantity, collapsedQuantity, k, block;

       for (i = 0, collapsedQuantity=0; i < n; i++) {
            for(j = 0, quantity=0; j < 2; j++){
                 block=0;
                 k = location(pos, h, sample[(2*i)+j], block);
                 if (pointer[h][sample[(2*i)+j]][0].mutation[block][k] == pos && k < pointer[h][sample[(2*i)+j]][0].mpb[block]) {
                     quantity++;
                 }
                 block=2;
                 if(pointer[h][sample[(2*i)+j]][0].b > 2){
                     k = location(pos, h, sample[(2*i)+j], block);
                     if (pointer[h][sample[(2*i)+j]][0].mutation[block][k] == pos && k < pointer[h][sample[(2*i)+j]][0].mpb[block]) {
                         quantity++;
                     }
                 }
            }
            if(quantity == 4){ collapsedQuantity+=2;}
            else{
                if(quantity !=0){ collapsedQuantity++;}
            }
       }
       return collapsedQuantity;
}

//This function gives the quantity (integer) of a particular mutation found in a sample
//arguments are h (which chromosome array to search), pos (the integer indicating the mutation position),
//block (the block in which to search), n (sample size), init (the first individual chromosome to search
//in the sample), and jumps (the intervals with which to search within the sample)
int muFrequencySampleIntDiscont(int h, int pos, int block, int s, int init, int jumps) {
    int i, quantity, k;

    for (i = init, quantity = 0; i < 2 * s; i += jumps) {
        k = location(pos, h, sample[i], block);
        if (pointer[h][sample[i]][0].mutation[block][k] == pos && k < pointer[h][sample[i]][0].mpb[block]) {
            quantity++;
        }
    }

    return quantity;
}

// SEARCH A GIVEN MUTATION INSIDE MUTTABLE (RETURNS ITS POSITION)
int SearchMutation(int block, int pos, int mutcount) {
    int counter;
    bool found;

    counter = 0;
    found = false;
    while (found == false && counter < mutcount) {
        if (muttable[counter].position == pos && muttable[counter].block == block) {
            found = true;
        }
        counter++;
    }

    return (counter - 1);
}

// NUMERICAL FUNCTION: TRACTPQL (RETURN THE TRACT LENGTH)
int tractpql(float meanTL)
{
    // DISTRIBUTION TYPICAL PARAMETERS
    float q = (float) 1 / meanTL;
    float threshold = 0;
    float p;
    int x = 0;
    float power;

    p = rand() / ((float) RAND_MAX + 1);
    do {
        x++;
        power = ((float) x) - 0.5;
        threshold += q * pow(1 - q, power - 1);
    } while (p > threshold);

    return x;
}

// FUNCTION THAT CREATES A LIST OF n CHROMOSOMES TO SAMPLE. WORKS AWFUL WHEN n IS CLOSE TO 2*N
void SamplingIndividuals(int n) {
    int count, val;
    bool samplelist[2 * N];

    if (n != N) {
        for (count = 0; count < 2 * N; count++) {
            samplelist[count] = false;
        }
        count = 0;
        while (count < n) {
            val = (int) (rand() % (N));
            if (samplelist[val] == false) {
                samplelist[2 * val] = true;
                samplelist[2 * val + 1] = true;
                count++;
            }
        }
        for (val = 0, count = 0; val < 2 * N; val++) {
            if (samplelist[val] == true) {
                sample[count] = val;
                count++;
            }
        }
    } else {
        for (count = 0; count < 2 * N; count++) {
            sample[count] = count;
        }
    }
}

///////////////////////////////
////  FIXATION TRAJECTORY  ////
///////////////////////////////

int GenerateFixationTrajectory(int maxTime, int fixationTime) {
    int return_var;
    if (fixationTime != 0){
        // GENERATE A FIXATION TRAJECTORY WITH A CONSTANT INCREASE IN POPOULATION SIZE AND FIXED TO fixationTime
        int tt;
        float interval = (float) 2 * N / fixationTime;
        float fixationTrajectoryFloat[fixationTime];

        fixationTrajectoryFloat[0] = 1;
        fixationTrajectory[0] = 1;
        for (tt = 1; tt < fixationTime; tt++) {
            fixationTrajectoryFloat[tt] = (float) fixationTrajectoryFloat[tt - 1] + interval;
            fixationTrajectory[tt] = (int) (fixationTrajectoryFloat[tt]);
        }
        for (tt = fixationTime; tt < STRUCTURED + 1; tt++) {
            fixationTrajectory[tt] = 2 * N;
        }

        return_var = fixationTime;
    }
    else if (fixationTime == 0){
        // THE FIXATION TRAJECTORY IS GENERATED THROUGH THE ALGORITHM PROPOSED BY KIMURA
        // IN PNAS 77 (1), 1980, USING A PSEUDO-SAMPLING METHOD.
        // THIS IS A MODIFIED VERSION THAT ALWAYS MAINTAINS A DIPLOID POPULATION WITH THE DUPLICATION
        int endTime=0, j, success, tt;
        float u, p, varP, psv;
        float fixationTrajectoryFloat[maxTime]; // Probability of fixation of the duplication in each generation (DupInd/2N)(relative freq)

        // Initialize duplication
        for (tt = 0; tt < 2; tt++) {
            fixationTrajectoryFloat[tt] = (float) 1 / (2 * N);
            fixationTrajectory[tt] = (int) (fixationTrajectoryFloat[tt]*2 * N);
        }
        success = 0;
        tt = 1;

        // Calculate trajectory
        do {
            u = rand() / ((float) RAND_MAX + 1);
            p = fixationTrajectoryFloat[tt];
            varP = p * (1 - p) / (2 * N);
            psv = sqrt(3 * varP)*(2 * u - 1);
            fixationTrajectoryFloat[tt + 1] = p + psv;
            tt++;
            fixationTrajectory[tt] = (int) (fixationTrajectoryFloat[tt]*2 * N);

            // MODIFICATION TO ALWAYS HAVE A DIPLOID POPULATION WITH THE DUPLICATION
            if (fixationTrajectory[tt] % 2 != 0) {
                fixationTrajectory[tt]--;
            }

            // If duplication is already fixed
            if (fixationTrajectoryFloat[tt] >= 1) {
                success = 1;
                endTime = tt;
                for (j = endTime; j < maxTime; j++) {
                    fixationTrajectory[j] = 2 * N;
                }
            }

            // If maxTime was not enough or duplication has been lost, repeat the generation of the trajectory (control to ensure fixation)
            if ((fixationTrajectory[tt] < 1) || (tt == maxTime)) {
                tt = 1;
                fixationTrajectoryFloat[tt] = (float) 1 / (2 * N);
            }
        } while (success == 0);



        return_var = endTime;
    }
    return return_var;
}


void print_fertility(){
    ofstream out;
    out.open("fertility.out");
    for (std::vector<fertility_info>::iterator it = fertility_list.begin() ; it != fertility_list.end(); ++it){
            //printf ("%d|%d, ",(*it).x,(*it).y);
            cout << (*it).x << " " << (*it).y << endl;
        }
    // for (int tt=0 ; tt < PROMETHEUS ; tt++) { for (int i=0 ; i < 2*N ; i++){
         //  	    	 cout << fertility[i][tt] << " ";} cout << "\n";} cout << "\n\n";
}

float round(float number_to_round, int decimal_places) // pow() doesn't work with unsigned
{
    return int(number_to_round * pow(10.0, decimal_places) + .50001) /  pow(10.0, decimal_places);
}

//void initialize_fertility_matrix(bool (*fertility)[2*N][PROMETHEUS]){
//	static const bool z[2*N][PROMETHEUS] = { false };
//	memcpy(fertility, &z, sizeof z);
//	for (int i=0 ; i < 2*N ; i++){
//		fertility[i][PROMETHEUS-1] = true;
//	}
//
//}


void sedus::requestWork()
{
    mutex.lock();
    _working = true;
    _abort = false;
    mutex.unlock();

    emit workRequested();
}
void sedus::abort()
{
    mutex.lock();
    if (_working) {
           _abort = true;
       }
    mutex.unlock();
}
 void sedus::setLog(const QString &value){
     emit addLog(value);
 }

 void sedus::setBar(int value){
     emit addBar(value);
 }
 void sedus::setChart(const qvdouble &x,const qvdouble &y){
     emit addChart(x,y);
 }
sedus::sedus(parameters *params, QObject *parent):QObject(parent)
//sedus::sedus(QObject *parent):QObject(parent)
{
    _working = false;
    _abort = false;
    using namespace std;

    setLog("Initializing...");

    int argumentscount =8;
    //main parameters
    letter = params->id;
    dir = params->dir;

    //exec parameters
    #undef SUPERTIME
    #define SUPERTIME sup
    sup=params->exec.runs;
    SAMPLE=params->exec.sample_size;
    PROMETHEUS = params->exec.generations;

    //main parameters
    N = params->main.N;
    THETA = params->main.theta;
    BLOCKLENGTH = params->main.blocklenght;
    timeToFixation=params->main.fixation_linear;
    TIMELENGTH = params->main.total * N; // Total number of generations (including all phases)
    BURNIN = params->main.burnin * N; // Number of generations in phase I
    STRUCTURED = (int) 20 * N; // Number of generations in phase II

    //igc parameters
    meanTractLength=params->igc.lambda;
    numbptotalidentity=params->igc.MEPS;

    R = params->crossover.R;
    rho = R / (4 * N * BLOCKLENGTH);
    C= params->igc.C;
    kappa = C / (4 * N * meanTractLength);

    //VARIABLES with some dependency with those prior
    //chrom::setfirst(B);
    //for(unsigned int i=0;i<chrom.mutation.size();i++){chrom.mutation[i].resize(BLOCKLENGTH);}
    //chrom ch1();
    //chrom ch2();
    //std::vector<chrom> ch1v;
    //table=ch1v;
    //pointer=ch2;
    multihit.resize(BLOCKLENGTH);
    table.resize(4*N);
    pointer.resize(2);for(unsigned int i=0;i<pointer.size();i++){pointer[i].resize(2*N);}
    fixationTrajectory.resize(20 * N + 1);
    mu = THETA / (4 * N);
    sampleN[0] = {SAMPLE};
    sample.resize(2*N);
    ancestry.resize(2 * N);for(unsigned int i=0;i<ancestry.size();i++){ancestry[i].resize(PROMETHEUS);}
    fertility.resize(2 * N);for(unsigned int i=0;i<fertility.size();i++){fertility[i].resize(PROMETHEUS);}
    fertility_ini.resize(2 * N);for(unsigned int i=0;i<fertility_ini.size();i++){fertility_ini[i].resize(PROMETHEUS);}
    recombimatrix.resize(2 * N);for(unsigned int i=0;i<recombimatrix.size();i++){recombimatrix[i].resize(PROMETHEUS);}
    duplicontent.resize(2);for(unsigned int i=0;i<duplicontent.size();i++){duplicontent[i].resize(2*N);}

    /*ancestry[2 * N][PROMETHEUS];
    fertility[2 * N][PROMETHEUS]; // Boolean Matrix codifying if a chromosome has descendants in the final generation of each era (each era is equivalent to PROMETHEUS generations) or not
    fertility_ini[2 * N][PROMETHEUS]; //Matrix fertility initialization
    recombimatrix[2 * N][PROMETHEUS];
    */
    //float timeToFixationFloat; sscanf (argv[5],"%f",&timeToFixationFloat);
    /*float timeToFixationFloat=params->main.fixation_linear;
    cout << timeToFixationFloat << endl;
    timeToFixation = (int) (timeToFixationFloat * N); // timeToFixation must be >= 1
    cout << timeToFixation << endl;*/
    //float argDonorRatio; sscanf (argv[6], "%f", &argDonorRatio);
    float argDonorRatio=params->igc.donor;
    donorRatio = argDonorRatio;
    argc = argumentscount;
    correctArguments = 0;
    //if (params->crossover.isSC){
        //if (argc == 8){
            crossoverBegin[0] = BLOCKLENGTH;
            crossoverEnd[0] = 2*BLOCKLENGTH;
            correctArguments = 1;
        //}
    //}

    /*else if (string(argv[8]) == "WR"){
        if (argc == 8){
            crossoverBegin[0] = 1;
            crossoverEnd[0] = B*BLOCKLENGTH;
            correctArguments = 1;
        }
    }
    else if (string(argv[8]) == "HS"){
        if (argc > 8){
            //int num; sscanf (argv[7],"%d",&num);
            int num=atoi(argv[7].c_str());
            numHS = num;
            if (argc == (8+numHS*2)){
                correctArguments = 1;
                for (int HS = 0; HS < numHS; HS++) {
                    int Beg; int End;
                    //sscanf (argv[8+HS*2],"%d",&Beg);
                    //sscanf (argv[9+HS*2],"%d",&End);
                    Beg=atoi(argv[8+HS*2].c_str());
                    End=atoi(argv[9+HS*2].c_str());
                    crossoverBegin[HS] = Beg;
                    crossoverEnd[HS] = End;
                }
            }
        }
    }*/
}


