// KMC program
/*
#######################################
#                                     #
#-- Kinetic Monte Carlo Simulation --#
#------  Author: Devosmita Sen  --------#
#                                     #
#######################################

*/
#include<iostream>
#include<cstdio>
#include<vector>
#include<algorithm>
#include<cmath>
#include<ctime>
#include<omp.h>
#include<map>
#include<utility>
#include<queue>
#include<algorithm>
#include<fstream>
#include<string>
#include<stdexcept>

#define MATRIX_A 0x9908b0dfUL //for RNG
#define UPPER_MASK 0x80000000UL
#define LOWER_MASK 0x7fffffffUL
#define DIST_TYPE short
using namespace std;
// Global variables
/*********************************/
/* Simulation Parameters */
/*********************************/
// chemistry related parameters


double R2_func()
{
  double n;
  ifstream infile("R2.txt");
  if (infile>>n)
    return n;
  else
    throw std::runtime_error("Cannot read size from file");
}
const double Nbsq = R2_func();


size_t NRA_func()
{
  size_t n;
  ifstream infile("NRA.txt");
  if (infile>>n)
    return n;
  else
    throw std::runtime_error("Cannot read size from file");
}
const size_t NRA = NRA_func();


double c0_in_func() // input C0
{
  double n;
  ifstream infile("C0.txt");
  if (infile>>n)
    return n;
  else
    throw std::runtime_error("Cannot read size from file");
}
const double c0_in = c0_in_func();



const double ca0=c0_in*2*0.001; // can be overwritten w/ argument input into main
// this definition of conc is number of functional groups / volume- 0.780- is from paper- and the factor 2 helps to convert into KMC definition
// 0.001 factor to convert from mM to M
const double c_star = ca0*0.5 * (6.02214e23/1.0e24)*pow(Nbsq,1.5);
const size_t MAX=(2<<sizeof(DIST_TYPE)*8-1); // maximal loop size considered
// /* A4B2 settings (comment this line to use this setting)
const size_t fA=2,fB=4; //changed to A2B2 now
const char dA=1,dB=0;
//const size_t molsizeA=1,molsizeB=0;
size_t step;

class bond {       // The class
  public:             // Access specifier
    size_t JA;        // Junction A
    size_t JB;  //Junction B
};

// */
/* A4B4 settings (comment this line to use this setting)
const size_t fA=4,fB=4;
const char dA=1,dB=1;
const size_t mA=1,mB=1;
// */
// basic simulation settings (final conversion, system size, loop?, etc.)
double conversion=1.0; // final conversion
//size_t NRA=7500; // number of A precursor
// can be overwritten w/ argument input into main
size_t NRB=NRA*fA/fB; // number of B precursor // changes dewpending on NA
const size_t NA=NRA*fA,NB=NRB*fB;
size_t start = NA>NB ? NB : NA;
size_t i_rxn=start;// index for tracking conversion
const bool LOOPFREE=false; // loop free or not, affects probability and final conversion
bool DEBUG = false; // debug mode, if TRUE, RNG will use the TEST_SEED specified
size_t TEST_SEED = 100; //
// output/environment settings
char PATH[]="";
size_t suffix;
bool writeLoop = true;
bool writeMW = false;
bool writeDistAA = false;
bool logDPdist = false;
bool logNW = false;
bool write_conv_loop=true;
bool write_connectivity_data=true;
char prefixLoopFrac[]="lpfrac_", prefixMW[]="MW_", prefixMWdist[]="cumDist_", prefixDegree[]="deg_", prefixDist[]="Dist_", prefixconv_loop_step1[]="1_Conv_loop_",prefixconv_loop_step2[]="2_Conv_loop_",prefixconv_loop_step3[]="3_Conv_loop_",prefixconnect[]="network_";
// logging parameters
size_t MWfreq=10; // output/record frequency for DPw-ie. output after every 10 steps
size_t loopVecSize = 2; // max order of loop to log
//size_t nNW = 2; // number of times to output network topology// why 2 times?
//size_t nDPdist = 2; // number of times to log MW distributions // why 2 times?
// RNG related
unsigned long mt[624]; // someything to do with RNG??
int mti=625;
double fn[128], wn[128];
int kn[128];
void RNG_initialize(unsigned long);
unsigned long rand32();
size_t seed; 
/********************************/
/* Simulation Variables */
/********************************/

//double sumrA_check=fA*NRA+1;
// variables relating to KMC
double PAB,la; // KMC propensities
int assoc_loop; // store infor about loop formation/breaking in assoc case
size_t tmpJA,tmpJB; // temporary containers for junction numbers
vector<double> p; // probabilities of different reactions

vector<vector<unsigned short> > JunctionDist; // array for pairwise topological distance
vector<vector<unsigned short> > LoopSize ; // array for pairwise loopsize-> if LoopSize[JA][JB]=1-> this means that JA and JB are involved in formation of a loop of size 1
vector<double> Conv; // conversions at which the following values were logged
vector<double> full_reactedA_array; // number of A junctions fully reacted
vector<double> all_reactedA_fg; // number of all A functional groups (not junctions) which are reacted
vector<double> Conv_calc; // conversions calculated using formula: conv=1-urA.size()/NA
vector<double> Sum; // stores the value of sum at each step // for dissoc
vector<double> Sum_assoc_exchng;// for associative- depend of probability of unreacted chain end (A) and reacted junction (B)
vector<bond> all_bonds;

vector<vector<size_t> > neighA,neighB; // neighbors of A(B) precursors

vector<size_t> urA,urB; // unreacted A and B junctions, size changing
vector<size_t> reactedA(0);//,reactedB; // unreacted A and B junctions, size changing
vector<size_t> reactedB(0);

vector<double> sumA; // sum of probability relating to each A junction// dissoc
vector<double> sumA_assoc_exchng; // associative

double sum; // sum of propensities// dissociative
double sum_assoc_exchng; // assoc exchange
// variable related to rejection condition
double V,BoxSize; // V = volume of simulation box (nm^3)
// BoxSize = V^(1/3) in unit of nm
//vector<vector<double> > rA,rB; // spatial coordinates of A and B molecules

vector<size_t> mol; // molecule idx for each junction
// size = number of Junctions (NRA+NRB)
size_t largestMol; // index of largest molecule


//vector<double> loop0; //
vector<double> loop1; //
vector<double> loop2; //
vector<double> loop3; //
vector<double> loop4; //
/*vector<double> loop5; //
vector<double> loop6; //
vector<double> loop7; //
vector<double> loop8; //
vector<double> loop9; //
vector<double> loop10; //
*/
// Connectivity data for network

vector<size_t> node1; 
vector<size_t> node2; 
// all with size = NRA*fA reserved
// containers for logging number of loops for loop counting at end of KMC
vector<double> loop; // number of loops- if there is fA=2 or fB=2, then this will be 2*num_loops
vector<double> loopfrac; // loop fraction


/*********************/
/* Functions */
/*********************/
// KMC functions
void initialize();
bool KMCstep_forward(); // implements KMC step- forward

double KMCconv(double ); // outputs the conversion- i think this conversion is same as the actual conversion
void output(size_t);
// helper functions
// general helper functions
size_t dist(unsigned char);
int getJunctionDistAA(size_t,size_t); //part of the shortest path update step??
int getJunctionDistBB(size_t,size_t);
// steps in KMCstep()
void SelectJunct_forward(size_t &,size_t &,size_t &,size_t &);// pair selection algorithm 

void UpdateConnectivityData(size_t &JA,size_t &JB); // update connectivity data (B4 node 1 connected to B4 node 2)


void UpdateSum_forward(const size_t,const size_t,const size_t,const size_t); // propensity sum update

void UpdateLoop_forward(const size_t,const size_t);// loop count update

void CollectConnected(const size_t,const size_t,vector<size_t>&,vector<size_t>&,vector<size_t>&,vector<size_t>&,vector<size_t>&,vector<size_t>&);

void UpdateJuncDist_forward(const size_t,const size_t,const vector<size_t>&,const vector<size_t>&);

void UpdateJuncDist_common(const size_t,const size_t,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&);//part of the shortest path update step??
//void UpdateJuncDist_reverse(const size_t,const size_t,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&);//part of the shortest path update step??
void UpdateJuncDist_common_reverse(const size_t,const size_t,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&);//part of the shortest path update step??

void UpdateMol_forward(const size_t,const size_t,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&,const vector<size_t>&);

// used in KMCconv()
void updateWriteData(double,double,size_t);

// MAIN FUNCTION
int main(int argc,char* argv[])// possible input arguements are- cstar, NRA, suffix
{
  //freopen("log_file", "w", stdout ); // opens file named "filename" for output 
   
 
    // handle argument input
    /*if(argc<2) { }
    else {
    double c_star0 = c_star;
    c_star = atof(argv[1]);// atof- string to float
    ca0 = ca0 * (c_star/c_star0);
    }
    if(argc==3) {
    NRA = atoi(argv[2]);// atoi- string to integer
    NRB = atoi(argv[2])*fA/fB;
    } 
    else if(argc==4) {
    NRA = atoi(argv[2]);
    NRB = atoi(argv[2])*fA/fB;
    suffix = atoi(argv[3]); // this is a global variable, so it is not being used in main(), but still can be used in initialize() function 
    }
    */
    // log run time
    clock_t c;
    c=clock();
    // RUN KMC simulation
    initialize();
   
    //step 1
   // cout<<"step1\n";
    //cout<<"urA.size()"<<urA.size();
   // cout<<"reactedA.size()"<<reactedA.size();

    step=1;
    Conv.clear();
    full_reactedA_array.clear();
    all_reactedA_fg.clear();
    Conv_calc.clear();
    Sum.clear();
    Sum_assoc_exchng.clear();
    //loop0.clear();
    loop1.clear();
    loop2.clear();
    loop3.clear();
    loop4.clear();
    /*
    loop5.clear();
    loop6.clear();
    loop7.clear();
    loop8.clear();
    loop9.clear();
    loop10.clear();*/

    node1.clear();
    node2.clear();
    double finalConv1 = KMCconv(conversion);
    cout<<"KMCconv done\n";
    output(step);
    cout<<"Hi\n";
    //FILE * fp;
    //fp=fopen("molecule_size1.csv","a");

    cout<<"c_star\n";
    printf("%.8f\t",c_star);
    cout<<"\n";
    cout<<"loop_frac\n";
    for(size_t l=0;l<loop.size();++l) printf("%.8f\t",loopfrac[l]);
    cout<<"\n";
    cout <<"finalConv1 "<< finalConv1 << "\n";

    cout<< "BoxSize " << BoxSize << "\n";

    cout << "seed " <<seed << "\n";
    cout<< (double)(clock()-c)/CLOCKS_PER_SEC<<endl;


    system("pause");
    exit(0);



    system("pause");
    return 0;
}

void initialize()
{
    // initialize RNG
    if(DEBUG)
    seed = TEST_SEED;
    else
    seed = time(NULL)+clock()+suffix;
    RNG_initialize(seed);
    // initialize PAB, la
    PAB=(1.0e24/6.02214e23)*pow((3.0/(2.0*3.14159*Nbsq)),1.5);
    la=PAB/ca0;
    // initialize box size and rejection related parameters
    V = (double)NRA*fA / (6.02214e23*ca0) * 1.0e24; // volume of simulation box in units of nm^3
    BoxSize = pow(V,1.0/3); // Length of simulation box in unit of nm
    double MaxTrials = NRA*fA/10; // not used anywhere else!!
    for(size_t i=0;i<NRA;++i)
    /*
    rA.push_back(vector<double>(3,0.0));// vector with size 3 and all values as 0.0- ie. all NA functional groups have a position given by a vector of length 3
    for(size_t i=0;i<NRA;++i)
        for(size_t x=0;x<3;++x)// 3 coordinates- x,y,and z
            rA[i][x] = rand32()/4294967296.0*BoxSize;// randomly defining the position of each of the crosslinks- where did this formula come from??
    for(size_t i=0;i<NRB;++i) // similarly for B
        rB.push_back(vector<double>(3,0.0));
        for(size_t i=0;i<NRB;++i)
            for(size_t x=0;x<3;++x)
                rB[i][x] = rand32()/4294967296.0*BoxSize;

    */
                // initialize matrix p for determining probability
    p = vector<double>(MAX+1,0);
    if(!LOOPFREE) {
        for(size_t i=1;i<MAX;++i)
            p[i] = 1.0 + la*NRA*fA*pow(i,-1.5);
        p[MAX] = 1.0;
    }         
    p[0] = 1.0;
    // initialize JunctionDist array
    for(size_t i=0;i<NRA;++i){
        JunctionDist.push_back(vector<unsigned short>(NRB,0));
        LoopSize.push_back(vector<unsigned short>(NRB,0));
    }

    assoc_loop=-5;
    // initialize neighA neighB
    for(size_t i=0;i<NRA;++i)
        neighA.push_back(vector<size_t>());
    for(size_t i=0;i<NRB;++i)
        neighB.push_back(vector<size_t>());
    // initialize unreacted A and unreacted B vectors
    for(size_t i=0;i<NRA;++i) {
        urA.push_back(i);
    }
        for(size_t i=0;i<NRB;++i) {
            urB.push_back(i);
        }
        // initialize cumulative probability
        sum = NRA*fA * NRB*fB;
        //cout<<"sum initially: "<<sum<<"\n";
        sumA = vector<double>(NRA,fA*NRB*fB);
        sumA_assoc_exchng=vector<double>(NRA,0); // initially propensity for assoc rxn =0
        sum_assoc_exchng = 0; // since initially the system has no reacted B junctions
        //cout<<"sum initially: "<<sum<<"\n";
        sumA = vector<double>(NRA,fA*NRB*fB);

        mol = vector<size_t>(NRA+NRB,0);
        largestMol=0;
        for(size_t i=0;i<mol.size();++i) mol[i] = i;
            size_t NNN = NRA*fA>NRB*fB ? NRB*fB : NRA*fA;  // ie. min(NRA*fA,NRB*fB)

        //loop0.reserve(NNN);
        loop1.reserve(NNN);
        loop2.reserve(NNN);
        loop3.reserve(NNN);
        
        loop4.reserve(NNN);
        /*
        loop5.reserve(NNN);
        loop6.reserve(NNN);
        loop7.reserve(NNN);
        loop8.reserve(NNN);
        loop9.reserve(NNN);
        loop10.reserve(NNN);*/
        node1.reserve(NNN);
        node2.reserve(NNN);
        // initialize loop vector
        loop = vector<double>(loopVecSize,0);
        loopfrac = vector<double>(loop.size(),0); // loop.size() should be equal to loopVecSize
        // initialize post KMC calculation variable
}

void output(size_t step)
{
    cout<<"wrote to file\n";
    FILE *fp;
    // set file name
    char fn[50]; // no more than 50 characters for file name
    // write loop fraction
    if(writeLoop) {
        sprintf(fn,"%s%scs=%1.4fA%dB%d.csv",PATH,prefixLoopFrac,c_star,NRA,NRB);
        //cout<<fn;
        fp=fopen(fn,"a");
        fprintf(fp,"%.8f,",c_star);
        for(size_t l=0;l<loop.size();++l)
        fprintf(fp,"%.8f,",loopfrac[l]);
        fprintf(fp,"\n");
        fclose(fp);
    }
    // write MW
    if(writeMW) {
        sprintf(fn,"%s%scs=%1.4fA%dB%d_%02d.csv",PATH,prefixMW,c_star,NRA,NRB,suffix);
        fp=fopen(fn,"a");
        fprintf(fp,"Conv,DPw,DPwr,X,SolFrac,BranchFrac,Loop1,Loop2,Loop3,Loop4\n");

        fprintf(fp,"\n");
        fclose(fp);
    }

    if(write_conv_loop) {
        if(step==1)
            sprintf(fn,"%s%scs=%1.4fA%dB%d_%02d.csv",PATH,prefixconv_loop_step1,c_star,NRA,NRB,suffix);
        else if(step==2)
            sprintf(fn,"%s%scs=%1.4fA%dB%d_%02d.csv",PATH,prefixconv_loop_step2,c_star,NRA,NRB,suffix);
        else if(step==3)
            sprintf(fn,"%s%scs=%1.4fA%dB%d_%02d.csv",PATH,prefixconv_loop_step3,c_star,NRA,NRB,suffix);
        fp=fopen(fn,"a");
        //fprintf(fp,"Conv,Loop1,Conv_calc,Sum,loop2,loop3,loop4,loop5,loop6,loop7,loop8,loop9,loop10,loop_sum,loop0\n");
        fprintf(fp,"Conv,Loop1,Conv_calc,Sum,Sum_assoc_exchng,full_reactedA_array,all_reactedA_fg,loop2\n");

        //cout<<"Conv.size()"<<Conv.size()<<"\n";
        //cout<<"Conv_calc.size()"<<Conv_calc.size()<<"\n";
        //cout<<"loop1.size()"<<loop1.size()<<"\n";
        for(size_t l=0;l<loop1.size();++l){

            fprintf(fp,"%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%8f\n",Conv[l],loop1[l],Conv_calc[l],Sum[l],Sum_assoc_exchng[l],full_reactedA_array[l],all_reactedA_fg[l],loop2[l]);
        }
        //printf("Conv[0]",Conv[0],"\n");
        fprintf(fp,"\n");
        fclose(fp);
    }


    if(write_connectivity_data) {

        sprintf(fn,"%s%sKMC.txt",PATH,prefixconnect);//,c_star,NRA,NRB,suffix);
        fp=fopen(fn,"a");
        fprintf(fp,"Node1,Node2\n");
        cout<<"node1.size()"<<node1.size()<<"\n";
        for(size_t l=0;l<node1.size();++l){
        //printf("Hello inside loop\n");
        // fprintf(fp,"HEllo");
            fprintf(fp,"%zu,%zu\n",node1[l],node2[l]);
        }
        //printf("Conv[0]",Conv[0],"\n");
        fprintf(fp,"\n");
        fclose(fp);
    }

}

bool KMCstep_forward()
{
   /* for (size_t A=0;A<NRA;++A){
        for (size_t i=0;i<neighA[A].size();++i){
            size_t B=neighA[A][i];
            if(mol[A]!=mol[B+NRA]){
                cout<<"A"<<A<<"\n";
                cout<<"B"<<neighA[A][i]<<"\n";
                cout<<"mol[A]"<<mol[A]<<"\n";
                cout<<"mol[B]"<<mol[B+NRA]<<"\n";
                cout<<"neighA[A].size()"<<neighA[A].size()<<"\n";
                for (size_t j=0;j<neighA[A].size();++j){
                    cout<<"neighA[A][i]"<<neighA[A][j]<<"\n";}
                cout<<"neighB[B].size()"<<neighB[B].size()<<"\n";
                for (size_t j=0;j<neighB[B].size();++j){
                    cout<<"neighB[B][i]"<<neighB[B][j]<<"\n";}
               // CollectConnected_reverse()
        vector<size_t> AcA_new,BcA_new,AcB_new,BcB_new;
        vector<size_t> dAcA_new,dBcB_new;
    // Collect junctions connected to JA in AcA BcA; connected to JB in AcB BcB
        CollectConnected_reverse(A,B,AcA_new,AcB_new,BcA_new,BcB_new,dAcA_new,dBcB_new); // reverse done- just for testing- to ensure that sum is not updated anywhere
        //cout<<"FOR 144-133:\n";
        vector<size_t> AcA_n,BcA_n,AcB_n,BcB_n;
        vector<size_t> dAcA_n,dBcB_n;
    // Collect junctions connected to JA in AcA BcA; connected to JB in AcB BcB
        CollectConnected_reverse(A,B,AcA_n,AcB_n,BcA_n,BcB_n,dAcA_n,dBcB_n); // reverse done- just for testing- to ensure that sum is not updated anywhere
        
        //cout<<"AcBnew.size()"<<AcB_new.size()<<"\n";
        //for (size_t j=0;j<AcB_new.size();++j){
           // cout<<"AcBnew[i]"<<AcB_new[j]<<"\n";}
        //cout<<"neighB[B].size()"<<neighB[B].size()<<"\n";
        //cout<<"JunctionDist[A][B]"<<JunctionDist[A][B]<<"\n";
        //cout<<"Forward-problem in molecule update somewhere!- case 1";
        mol[A]=mol[B+NRA];
        //cout<<"AFTER CHANGE";
        //cout<<"mol[A]"<<mol[A]<<"\n";
        //cout<<"mol[B]"<<mol[B+NRA]<<"\n";
        //system("pause");
        
        //cout<<"AcBn.size()"<<AcB_n.size()<<"\n";
        //for (size_t j=0;j<AcB_n.size();++j){
            //cout<<"AcBn[i]"<<AcB_n[j]<<"\n";}
        //cout<<"neighB[B].size()"<<neighB[B].size()<<"\n";
        //cout<<"JunctionDist[A][B]"<<JunctionDist[A][B]<<"\n";
        //cout<<"Forward-problem in molecule update somewhere!- case 1";
        
        system("pause");
            }
        }
    }*/

    /*for (size_t B=0;B<NRB;++B){
        for (size_t i=0;i<neighB[B].size();++i){
            size_t A=neighB[B][i];
            if(mol[A]!=mol[B+NRA]){
                cout<<"B"<<B<<"\n";
                cout<<"A"<<neighB[B][i]<<"\n";
                cout<<"mol[A]"<<mol[A]<<"\n";
                cout<<"mol[B]"<<mol[B+NRA]<<"\n";
                cout<<"neighA[A].size()"<<neighA[A].size()<<"\n";
                cout<<"neighA[A][i]"<<neighA[A][i]<<"\n";
                cout<<"neighB[B].size()"<<neighB[B].size()<<"\n";
                cout<<"JunctionDist[A][B]"<<JunctionDist[A][B]<<"\n";
                
                cout<<"Forward-problem in molecule update somewhere!- case 2\n";
                mol[A]=mol[B+NRA];
                cout<<"AFTER CHANGE";
                cout<<"mol[A]"<<mol[A]<<"\n";
                cout<<"mol[B]"<<mol[B+NRA]<<"\n";
                system("pause");
            }
        }
    }*/
    // JA JB holds the selected junction for this step
    // idxA idxB holds the index of JA JB in urA urB
    size_t JA,JB,idxA,idxB;
    //size_t numTrials = 0; // not used anyuwhere!!
    // select the junction A and junction B to react
    SelectJunct_forward(JA,JB,idxA,idxB); // JA,JB,idxA,idxB are passed as pointers- so, these are not initialized earlier
    UpdateConnectivityData(JA,JB); // write to file- the junctions which react

    // Updates the sum of relative probabilities of unreacted A-B pairs,
    //by deleting the probabilities for the A and B reacting this step
    //cout<<"JA"<<JA<<"\n";
    //cout<<"JB"<<JB<<"\n";
    UpdateSum_forward(JA,JB,idxA,idxB);
    //cout<<urA.size()<<"\n";
    // Update Loop information
    UpdateLoop_forward(JA,JB);
    // Find all junctions connected to JA and JB
    // AcA holds all A junctions that is connected to JA
    // BcA holds all B junctions that is connected to JA
    // AcB holds all A junctions that is connected to JB
    // BcB holds all B junctions that is connected to JB
    vector<size_t> AcA,BcA,AcB,BcB;
    vector<size_t> dAcA,dBcB;
    // Collect junctions connected to JA in AcA BcA; connected to JB in AcB BcB
    CollectConnected(JA,JB,AcA,AcB,BcA,BcB,dAcA,dBcB);
    // Updates connectivity by updating neighbor list of JA and JB
    neighA[JA].push_back(JB);
    neighB[JB].push_back(JA);
    bond new_bond;
    new_bond.JA=JA;
    new_bond.JB=JB;
    all_bonds.push_back(new_bond);
    //cout<<"update neighbors for JA and JB: "<<JA<<"  "<<JB<<"\n";
    //cout<<"neighA[JA].size()"<<neighA[JA].size()<<"\n";
    //cout<<"neighB[JB].size()"<<neighB[JB].size()<<"\n";
    // update JunctionDist
    UpdateJuncDist_forward(JA,JB,BcA,AcB); // same function definition, but different inputs- need to decide which one is 
    // being used based on the inputs only
    UpdateJuncDist_common(JA,JB,AcA,BcB,dAcA,dBcB);
    // Update the molecule grouping info and molecule sizes
   // cout<<"before molecule update\n";
    //cout<<"mol[JA]"<<mol[JA]<<"\n";
    //cout<<"mol[JB]"<<mol[JB+NRA]<<"\n";
    //size_t mA1=mol[JA];
    //size_t mB1=mol[JB+NRA];
    UpdateMol_forward(JA,JB,AcA,AcB,BcA,BcB);
    //cout<<"after molecule update\n";
    //cout<<"mol[JA]"<<mol[JA]<<"\n";
    //cout<<"mol[JB]"<<mol[JB+NRA]<<"\n";
    //size_t mA2=mol[JA];
    //size_t mB2=mol[JB+NRA];
   /* if((mA1!=mB1) && (mA1==mA2) && (mB1==mB2)){
        cout<<"PROBLEM in moleucle update\n";
        system("pause");
    }*/
   

    /*cout<<"AcA[j]:";
        for(size_t j=0;j<AcA.size();++j) cout<<AcA[j]<<" ";
        cout<<endl;
    cout<<"BcA[j]:";
        for(size_t j=0;j<BcA.size();++j) cout<<BcA[j]<<" ";
        cout<<endl;
    cout<<"AcB[j]:";
    for(size_t j=0;j<AcB.size();++j) cout<<AcB[j]<<" ";
        cout<<endl;
    cout<<"BcB[j]:";
        for(size_t j=0;j<BcB.size();++j) cout<<BcB[j]<<" ";
        cout<<endl;
    cout<<"JunctionDist[JA][JB]"<<JunctionDist[JA][JB]<<endl;
    cout<<"mol[101]"<<mol[101]<<endl;
    cout<<"mol[43]"<<mol[43]<<endl;
        //system("pause");
    */
     /*for(size_t i=0;i<reactedA.size();++i){
        if(neighA[reactedA[i]].size()<fA){
            cout<<"in forward neighA[reactedA[i]].size<fA\n";
            cout<<"in forward reactedA[i]"<<reactedA[i]<<"\n";
            cout<<"neighA[reactedA[i]].size()"<<neighA[reactedA[i]].size()<<endl;
            system("pause");
        }
    }*/
    return true;
}











double KMCconv(double conversion)
{

    // do KMCstep until designated conversion
    //NA=NRA*fA,NB=NRB*fB;
    double End = (1-conversion)*start;
    if(LOOPFREE) End = 1 + start - ( start/fA+start/fB );
    double currConv=1.0;
    double currConv_calc=1.0;
    size_t totalbond = start - ceil(End);
    //size_t i_rxn=start;// index for tracking conversion
    size_t i_temp;
    size_t i_start=start;
    for( size_t i = i_start; i > End; i--) {
    /*
    double sumfA=0;//sum of all reacted B "functional groups"- useful for forward rxns
    double sumrA=0;// sum of all reacted A "functional groups"- useful for reverse rxns
    for(size_t k=0; k<reactedA.size(); ++k){
        size_t num_reactedA=neighA[reactedA[k]].size();
        sumrA=sumrA+num_reactedA;
    }
    for(size_t k=0; k<urA.size(); ++k){
        size_t num_reactedA=neighA[urA[k]].size();
        size_t num_unreactedA=fA-num_reactedA;
        sumfA=sumfA+num_unreactedA;
        sumrA=sumrA+num_reactedA;
    }
    double sumfB=0;
    double sumrB=0;// sum of all reacted B "functional groups"
    for(size_t k=0; k<reactedB.size(); ++k){
        size_t num_reactedB=neighB[reactedB[k]].size();
        sumrB=sumrB+num_reactedB;
    }
    for(size_t k=0; k<urB.size(); ++k){
        
        size_t num_reactedB=neighB[urB[k]].size();
        size_t num_unreactedB=fB-num_reactedB;
        sumfB=sumfB+num_unreactedB;
        sumrB=sumrB+num_reactedB;
    }

    if(sumrA!=sumrB) {
        cout<<"PROBLEM!! sumrA!=sumrB"<<"\n";
        cout<<"sumrA: "<<sumrA<<"\n";
        cout<<"sumrB "<<sumrB<<"\n";
        cout<<"sumfA: "<<sumfA<<"\n";
        cout<<"sumfB "<<sumfB<<"\n";
        system("pause");
    }
    */

    double sum_f=sum;//- THIS SHOULD BE EQUAL TO THE SUM (Sum of probabilities of unreacted junctions)- AND NOT MERELY THE NUMBER OF BONDS AVAILABLE FOR REACTION
   
    double sum_r=all_bonds.size();//sumrA; // // THIS SHOULD BE EQUAL TO THE NUMBER OF AVAILABLE BONDS!! - NOT THIS!!


    double P_f;
    double P_a;
    
    //if(invKeq==0.0) { // only dissoc forward and assoc exchng
        //P_f=((double) sum_f/(sum_f+sum_a)<1.0? (double) sum_f/(sum_f+sum_a) : 1.0);
        //min((double) sum_f/(sum_f+sum_a), 1.0);
       // cout<<"Pf"<<P_f;
        //system("pause");
        //P_a=1-P_f;
        P_f=1;
    //}
    
    //P_a=0;
    double rand_num=rand32()/4294967296.0;
    
    if(rand_num<P_f){
  

    if(KMCstep_forward()) {
        //cout<<"KMCstep_forward() done\n";
        i_rxn--;
        currConv = (double)(start-i_rxn+1)/(double)start;
            

    //currConv_calc= 1.0-(double) urA.size()/(NRA);
    double num_ur_A_junc=0;
    for (size_t j=0;j<urA.size();++j){
            num_ur_A_junc+=fA-neighA[urA[j]].size();
    }
    currConv_calc= 1.0-(double) num_ur_A_junc/((double)(NA));
   
    if(!((start-i+1)%MWfreq)) {//cout<<"writing to file at this i: Rxn1: "<<i<<"\n";
        updateWriteData(currConv,currConv_calc, totalbond);// basically means that do this step when (start-i+1)%MWfreq) is 0
    }
    if(sum < 1) {cout<<"breaking KMCConv -forward because sum<1\n";cout<<"sum:"<<sum<<"\n";break;}
    } 
    else { cout<<"KMCstep_forward is false"; break; }
    }


    else{ // reverse rxn
        //cout<<"reverse rxn\n";
    

    //doing reverse reaction
    
    }
  
    i_temp=i;
    
    }
    
    cout<<"total number of rxn steps: i_start-i_temp "<<i_start-i_temp<<endl;
    cout<<"step: "<<step<<endl;cout<<"i_temp: "<<i_temp<<endl;cout<<"i_start: "<<i_start<<endl;
    cout<<"currConv_calc"<<currConv_calc<<endl;
    // calculate loop fraction from loop[]
    cout<<"all rxn steps done\n";
    for(size_t l=0;l<loop.size();++l)
        loopfrac[l]=((double)loop[l])/totalbond;
        return currConv;
}

void UpdateConnectivityData(size_t &JA,size_t &JB) //TODO
{
char fn[50];
//sprintf(fn,"network.csv");//,PATH,"NW_",c_star,currConv,NRA,NRB,suffix);
//FILE *fp;
//fp=fopen(fn,"a");

if(neighA[JA].size()>0){// means there is at least one B4 connected to A2
    size_t n1=neighA[JA][0] ; // node1-there can be at most one B node connected to A
    size_t n2=JB; // node2- to be connected
    //fprintf(fp,"%d,%d\n",n1,n2);
    node1.push_back(n1);
    node2.push_back(n2);
    //cout<<"updated to list "<<n1<<" "<<n2<<"\n";
    //system("pause");
}
// if neighA[JA].size()==0 then this means that a free A node is just getting connected with B, and hence, no change in connectivity pattern of B
if(neighA[JA].size()==2){

    cout<<"PROBLEM in junction selection\n";
    system("pause");
}
/*
for(size_t i=0;i<neighA.size();++i) {
fprintf(fp,"A%d",i);
for(size_t j=0;j<neighA[i].size();++j) {
fprintf(fp,",B%d",neighA[i][j]);
}
}*/
}
void updateWriteData(double currConv,double currConv_calc, size_t totalbond)
{
    //size_t NA=NRA*fA,NB=NRB*fB;
    double totalWW=0,totalW=0;

    //loop0.push_back(loop[0]/totalbond);
    loop1.push_back(loop[1]/totalbond);
    loop2.push_back(loop[2]/totalbond);
    loop3.push_back(loop[3]/totalbond);
    loop4.push_back(loop[4]/totalbond);/*
    loop5.push_back(loop[5]/totalbond);
    loop6.push_back(loop[6]/totalbond);
    loop7.push_back(loop[7]/totalbond);
    loop8.push_back(loop[8]/totalbond);
    loop9.push_back(loop[9]/totalbond);
    loop10.push_back(loop[10]/totalbond);*/

    //cout<<"loop1 value"<<loop[1]/totalbond<<"\n";
    Conv.push_back(currConv);
    full_reactedA_array.push_back(reactedA.size());
    Conv_calc.push_back(currConv_calc);
    Sum.push_back(sum);
    Sum_assoc_exchng.push_back(sum_assoc_exchng);

    double num_react_A=0; // number of reacted groups
    for (size_t j=0;j<urA.size();++j){
            num_react_A+=neighA[urA[j]].size();
    }
    num_react_A+=fA*reactedA.size();
    all_reactedA_fg.push_back(num_react_A);
}

void SelectJunct_forward(size_t &JA,size_t &JB,size_t &idxA,size_t &idxB)
{
    double stop = rand32()/4294967296.0*sum; //rand32()/4294967296.0 is a uniformly distributed random number between 0 and 1
    //stop = 0.5*sum;// this is only for debugging- delete this line later!!
    double cump = 0;
    int JA_found=0;
    int JB_found=0;
    for(size_t i=0;i<urA.size();++i) {
        if(cump+sumA[urA[i]] >= stop) {
        // cout<<"JA done: sumA[urA[i]]="<<sumA[urA[i]]<<"\n";
        JA = urA[i];
        idxA = i;
        JA_found=1;
        break;
        } else {
        cump += sumA[urA[i]];
        //cout<<"sumA[urA[i]]"<<sumA[urA[i]]<<"\n";

        }
    }
    if(JA_found==0){
        cout<<"JA_found=0\n";
        cout<<"stop"<<stop<<"\n";
        cout<<"cump"<<cump<<"\n";
        cout<<"urA.size()"<<urA.size()<<"\n";
        system("pause");
    }
    for(size_t i=0;i<urB.size();++i) {
        JB = urB[i];
        double tmpP = p[dist(JunctionDist[JA][JB])] * (double)(fB - neighB[JB].size()) * (double)(fA -neighA[JA].size());// this is degeneracy
        if(cump+tmpP >= stop) {
        idxB = i;
        JB_found=1;
        break;
        } else {
        cump += tmpP;
        //cout<<"tmpP"<<tmpP<<"\n";
        }
    }
    if(JB_found==0){
        cout<<"JB_found=0\n";
        cout<<"stop"<<stop<<"\n";
        cout<<"cump"<<cump<<"\n";
        cout<<"sum"<<"\n";
        cout<<"urB.size()"<<urB.size()<<"\n";
        system("pause");
    }
}

void SelectJunct_assoc(size_t &JA,size_t &JB,size_t &JA_conn_JB, size_t &idxA,size_t &idxB)
{

    /*for(size_t i=0;i<reactedA.size();++i){
        if(neighA[reactedA[i]].size()<fA){
            cout<<"before assoc neighA[reactedA[i]].size<fA"<<endl;
            cout<<"before assoc reactedA[i]"<<reactedA[i]<<"\n";
            cout<<"reactedA[i]"<<reactedA[i]<<endl;
            system("pause");
        }
    }*/
    //cout<<"before junction selection\n";
    //cout<<"urA.size()"<<urA.size()<<"\n";

    double stop = (rand32()/4294967296.0)*sum_assoc_exchng; //rand32()/4294967296.0 is a uniformly distributed random number between 0 and 1
    //stop = 0.5*sum;// this is only for debugging- delete this line later!!
    double cump = 0;
    int JA_found=0;
    int JB_found=0;
    for(size_t i=0;i<urA.size();++i) {
        if(cump+sumA_assoc_exchng[urA[i]] >= stop) {
            // cout<<"JA done: sumA[urA[i]]="<<sumA[urA[i]]<<"\n";
            JA = urA[i];
            idxA = i;
            JA_found=1;
            //cout<<"found JA junction assoc\n";
            break;
    } 
    else {
        cump += sumA_assoc_exchng[urA[i]];
        //cout<<"sumA[urA[i]]"<<sumA[urA[i]]<<"\n";
        }
    }

    if(JA_found==0){
        cout<<"JA_found=0 in assoc\n";
        cout<<"stop"<<stop<<"\n";
        cout<<"cump"<<cump<<"\n";
        cout<<"urA.size()"<<urA.size()<<"\n";
        system("pause");
    }
    for(size_t i=0;i<NRB;++i) { // B junction to be selected should be reacted- because here we are doing associative exchange
        if(neighB[i].size()>0){
        JB =i;
        double tmpP = p[dist(JunctionDist[JA][JB])] * (double)(neighB[JB].size()) * (double)(fA -neighA[JA].size());// this is degeneracy
        if(cump+tmpP >= stop) {
            idxB = i;
            JB_found=1;
            //cout<<"found JB junction assoc\n";
            break;
        } 
        else {
            cump += tmpP;
        }
        }
    }
    if(JB_found==0){
        cout<<"JB_found=0 in assoc\n";
        cout<<"stop"<<stop<<"\n";
        cout<<"cump"<<cump<<"\n";
        cout<<"sum"<<"\n";
        cout<<"urB.size()"<<urB.size()<<"\n";
        system("pause");
    }



    double rand_num=rand32()/4294967296.0;
    //cout<<"found rand_num\n";
    //cout<<"neighB[JB].size() "<<neighB[JB].size()<<endl;
    //cout<<"rand_num*neighB[JB].size() "<<(size_t) rand_num*neighB[JB].size()<<endl;
    //system("pause");
    int JB_index=rand_num*neighB[JB].size();// which index of functional group of JB i am selecting for exchanging
    // currently considered to be randomly selected
    JA_conn_JB=neighB[JB][JB_index]; // this is the A group which was connected to JB before exchange( this one got kicked out due to exchange)
    //cout<<"found JA_conn_JB\n";


    size_t count=0;
    for(size_t i=0;i<neighA[JA_conn_JB].size();++i){
        if(neighA[JA_conn_JB][i]==JB) count++;
    }
    if(dist(JunctionDist[JA][JB])==1 && count==fA){//LoopSize[JA_conn_JB][JB]==1){ // if JA and JA_conn_JB belong to same JB and JA_conn_JB-JB form a loop
        assoc_loop=2;// intermolecular loop interchange
    }
    else if(JA==JA_conn_JB){ // if JA=JA_conn_JB- then nothing happens- similar to earlier case
        assoc_loop=3; // in this case- it is just one chain flipping its orientation- so no need to update anything (including distances, etc)
    }
    else if(dist(JunctionDist[JA][JB])==1){
        // // JA and JB were neighbors before exchange reaction 
        assoc_loop=1;
    }
    else{
        if(count<fA){//LoopSize[JA_conn_JB][JB]==0){
            assoc_loop=0; // no change in loop density
        }
        else if(count==fA){//LoopSize[JA_conn_JB][JB]==1){
            assoc_loop=-1;// loop density decreases- loop breaking
        }
        else{
            cout<<"PROBLEM in assoc_loop!!\n";
            system("pause");
        }
    }
    //cout<<"assoc_loop"<<assoc_loop<<"\n";
    //system("pause");
    if(assoc_loop!=3){ // no change needed if assoc_loop=3
    









    //cout<<"before urA update\n";
    //cout<<"urA.size()"<<urA.size()<<"\n";
    //cout<<"now updating urA and reactedA arrays\n";
    // updating urA and reactedA arrays
    if(neighA[JA].size()==fA-1) { //before the assoc rxn- JA was just one away from being fully reacted
        //cout<<"JA"<<JA<<endl;
        //cout<<"urA[idxA]"<<urA[idxA]<<endl;
        reactedA.push_back(JA);
        //cout<<"updating urA \n";
        urA.erase(urA.begin()+idxA);
        
        /*if(JA==185) {cout<<"assoc JA=185"<<endl;
        cout<<"neighA[185].size()"<<neighA[185].size()<<endl;
        //system("pause");
        }*/
    }
    //cout<<"after first urA update\n";
    //cout<<"urA.size()"<<urA.size()<<"\n";
    //cout<<"reactedA.size()"<<reactedA.size()<<"\n";

    //cout<<"updated arrays arrays\n";
    // now that JB is selected, i have to select which functional group in junction JB gets exchanged- and then accordingly i have to update the details about the chain associated with that functional group in crosslinker JB
    
    if(neighA[JA_conn_JB].size()==fA) { // was fully connected earlier
        reactedA.erase(std::remove(reactedA.begin(),reactedA.end(), JA_conn_JB),reactedA.end());
        //cout<<"reactedA.size()"<<reactedA.size()<<"\n";
       // if(JA_conn_JB==185){cout<<"185 deleted from reactedA"<<endl;system("pause");
        //}

        urA.push_back(JA_conn_JB);
    }
    /*cout<<"after second urA update\n";
    cout<<"urA.size()"<<urA.size()<<"\n";
    cout<<"reactedA.size()"<<reactedA.size()<<"\n";
    cout<<"after urA update\n";
    cout<<"urA.size()"<<urA.size()<<"\n";
    cout<<"assoc rxn\n";
    cout<<"JA:"<<JA<<" JB:"<<JB<<" JA_conn_JB:"<<JA_conn_JB<<"\n";*/

    /*cout<<"neighA[JA]\n";
    for(size_t i=0;i<neighA[JA].size();++i){cout<<neighA[JA][i]<<" ";}
    cout<<endl;

    cout<<"neighA[JA_conn_JB]\n";
    for(size_t i=0;i<neighA[JA_conn_JB].size();++i){cout<<neighA[JA_conn_JB][i]<<" ";}
    cout<<endl;

    cout<<"neighB[JB]\n";
    for(size_t i=0;i<neighB[JB].size();++i){cout<<neighB[JB][i]<<" ";}
    cout<<endl;*/

    // from the junctions JA and JB and JA_conn_JB selected- and their state before rxn- we can decide whether the rxn is loop forming, loop breaking, or none
    // if JA connected to JB- loop forming : assoc_loop=1
    // if JA not connected to JB: then 2 cases:
        // 1. LoopSize[JA_conn_JB][JB]=0- then there was no loop, hence no loop update: assoc_loop=0
        // 2. LoopSize[JA_conn_JB][JB]=1- then loop breaking: assoc_loop=-1

    /*for (size_t i=0;i<neighA[JA].size();++i){
        if(neighA[JA][i]==JB){
            // JA and JB were neighbors before exchange reaction 
            assoc_loop=1;
            break;
        }
    }*/

    //cout<<"updated urA and reactedA arrays\n";
    //cout<<"LoopSize[JA_conn_JB][JB]"<<LoopSize[JA_conn_JB][JB]<<endl;

    

    // updating neighbours

    neighA[JA].push_back(JB); // JA and JB get connected
   // if(JA==185){cout<<"after neigh update- neighA[185].size(): "<<neighA[JA].size()<<endl; system("pause");}
    for (size_t i=0; i<neighB[JB].size(); ++i){
        if(neighB[JB][i]==JA_conn_JB){
            neighB[JB][i]=JA; // JA replaces JA_conn_JB as the neighbors of JB in the exchange process
            break;// since we want to replace only one bond
        }
    }


   // cout<<"JA_conn_JB"<<JA_conn_JB<<endl;
    //cout<<"JB"<<JB<<endl;
    //cout<<"neighA[JA_conn_JB].size() before"<<neighA[JA_conn_JB].size()<<"\n";
    /*cout<<"urA elements:\n";
    for (size_t i=0; i<urA.size(); ++i){
        cout<<urA[i]<<" ";
    }*/

   /* cout<<"neighA[JA_conn_JB][i]:\n";//<<neighA[JA_conn_JB][i]<<endl;
    for (size_t i=0; i<neighA[JA_conn_JB].size(); ++i){
        cout<<neighA[JA_conn_JB][i]<<" ";
    }
    cout<<endl;
    cout<<"neighB[JB][i]:\n";//<<neighA[JA_conn_JB][i]<<endl;
    for (size_t i=0; i<neighB[JB].size(); ++i){
        cout<<neighB[JB][i]<<" ";
    }
    cout<<endl;
    */
    for (size_t i=0; i<neighA[JA_conn_JB].size(); ++i){
        if(neighA[JA_conn_JB][i]==JB){
            neighA[JA_conn_JB].erase(neighA[JA_conn_JB].begin()+i); // JA_conn_JB is kicked out of connection with JB
            /*if(JA_conn_JB==185){cout<<"neighA[185] deleted in assoc"<<endl; 
            cout<<"reactedA[i]:"<<endl;
            for(size_t i=0;i<reactedA.size();++i){
                cout<<reactedA[i]<<" "<<endl;
            }
            system("pause");
            }*/
            //cout<<"i"<<i<<"\n";
            break;
        }
    }
    //cout<<"neighA[JA_conn_JB].size() after"<<neighA[JA_conn_JB].size()<<"\n";
    //cout<<"neighA[JA_conn_JB][i]:\n";//<<neighA[JA_conn_JB][i]<<endl;

    /*for (size_t i=0; i<neighA[JA_conn_JB].size(); ++i){
        cout<<neighA[JA_conn_JB][i]<<" ";
    }
    cout<<endl;*/
    /*for (size_t i=0; i<neighA[JA_conn_JB].size(); ++i){
        if(neighA[JA_conn_JB][i]==JB){
            cout<<"problem with neighbour update!! in assoc\n";
            cout<<"i"<<i<<"\n";
            system("pause");
        }
    }*/
    // update bond info

    // form bond between JA and JB
    bond new_bond;
    new_bond.JA=JA;
    new_bond.JB=JB;
    all_bonds.push_back(new_bond); // for forward rxn
    // delete bond between JA_conn_JB and JB
    for (size_t i=0;i<all_bonds.size(); ++i){
        if(all_bonds[i].JA==JA_conn_JB && all_bonds[i].JB==JB){
            all_bonds.erase(all_bonds.begin()+i);
            break;
        }
    }
    }
    /*cout<<"after junction selection\n";
    cout<<"urA.size()"<<urA.size()<<"\n";

    cout<<"AFTER neigh UPDATE!!\n";
    cout<<"neighA[JA]\n";
    for(size_t i=0;i<neighA[JA].size();++i){cout<<neighA[JA][i]<<" ";}
    cout<<endl;

    cout<<"neighA[JA_conn_JB]\n";
    for(size_t i=0;i<neighA[JA_conn_JB].size();++i){cout<<neighA[JA_conn_JB][i]<<" ";}
    cout<<endl;

    cout<<"neighB[JB]\n";
    for(size_t i=0;i<neighB[JB].size();++i){cout<<neighB[JB][i]<<" ";}
    cout<<endl;*/

    /*for(size_t i=0;i<reactedA.size();++i){
        if(neighA[reactedA[i]].size()<fA){
            cout<<"in assoc neighA[reactedA[i]].size<fA"<<endl;
            cout<<"in assoc reactedA[i]"<<reactedA[i]<<"\n";
            cout<<"reactedA[i]"<<reactedA[i]<<endl;
            system("pause");
        }
    }*/
}



void UpdateSum_forward(const size_t JA,const size_t JB,const size_t idxA,const size_t idxB)
{
    double probAB=0;
    for(size_t i=0;i<urA.size();++i) {
        tmpJA = urA[i];
        probAB = p[dist(JunctionDist[tmpJA][JB])] * (double)(fA - neighA[tmpJA].size());
        sum -= probAB;
        sumA[tmpJA] -= probAB;
    }
    for(size_t i=0;i<urB.size();++i) {
        tmpJB = urB[i];
        probAB = p[dist(JunctionDist[JA][tmpJB])] * (double)(fB - neighB[tmpJB].size());
        sum -= probAB;
        sumA[JA] -= probAB;
    }
    probAB = p[dist(JunctionDist[JA][JB])];
    sum += probAB;
    sumA[JA] += probAB;
    // erase element for urA/urB
    // if size == fA-1 or fB-1, then all nodes have been reacted for the junction
    // neighbor list not updated here because it affects the calculation of original distance between junctions
    if(neighA[JA].size() == fA-1) {
        reactedA.push_back(JA);
        //cout<<"JA "<<JA<<" pushed to reactedA"<<endl;
        urA.erase(urA.begin()+idxA);
        //cout<<"JA "<<urA[idxA]<<" erased from unreactedA"<<endl;
        //if(JA==185) {cout<<"forward"<<endl;system("pause");}
    }
    if(neighB[JB].size() == fB-1) {
        reactedB.push_back(JB);
        urB.erase(urB.begin()+idxB);
    }
    
}


void UpdateLoop_forward(const size_t JA,const size_t JB)
{
    for(size_t i=0;i<loop.size();++i) {
        if(dist(JunctionDist[JA][JB]) == i) {
        // loop[0] is for new branch, not loop  --> that is why ++i is done
        // for system with fA or fB =2, loop[n] is for n/2 order loop
        // for other system, loop[n] is for nth order loop
        // loop[n] should be zero for odd n
            loop[i]++;
            LoopSize[JA][JB]=i;
            break;
        }
    }
}

void UpdateLoop_assoc(const size_t JA,const size_t JB,const size_t JA_conn_JB)//, const size_t assoc_loop)
{
    if(assoc_loop==1){
        for(size_t i=0;i<loop.size();++i) {
            if(dist(JunctionDist[JA][JB]) == i) {
            // loop[0] is for new branch, not loop  --> that is why ++i is done
            // for system with fA or fB =2, loop[n] is for n/2 order loop
            // for other system, loop[n] is for nth order loop
            // loop[n] should be zero for odd n
                loop[i]++;
                LoopSize[JA][JB]=i;
                break;
            }
        }
    }
    else if(assoc_loop==-1){
        if(LoopSize[JA_conn_JB][JB]==1){
            loop[LoopSize[JA_conn_JB][JB]]--;// LoopSize[JA][JB]- is the loopsize in which JA and JB were involved before breaking
    //cout<<"Updated loop fraction of size"<<LoopSize[JA][JB]<<"\n";
        }
        LoopSize[JA_conn_JB][JB]=0;// because after bond breaking, there is no more loop between JA_conn_JB and JB- hence LoopSize=0
        LoopSize[JA][JB]=0;
    }
    else if(assoc_loop==2){
        LoopSize[JA][JB]++; // should be equal to 1
        if(LoopSize[JA][JB]!=1){cout<<"PROBLEM in assoc_loop=2\n"; system("pause");}
        LoopSize[JA_conn_JB][JB]=0;
    }
}



void CollectConnected(const size_t JA,const size_t JB,vector<size_t> &AcA,vector<size_t> &AcB,
vector<size_t> &BcA,vector<size_t> &BcB,vector<size_t> &dAcA,vector<size_t> &dBcB)
{   
    //cout<<"in collectconnected_forward"<<endl;
   // system("pause");
    int tmpDist;
    // collecting A junctions connected to JA/JB
    //cout<<"JunctionDist[278][188]"<<JunctionDist[278][188]<<endl;
    for(size_t i=0;i<NRA;++i) {
        tmpJA = i;
        // collecting A junction connected to JB
        if(JunctionDist[tmpJA][JB])
            AcB.push_back(tmpJA);
        // collecting A junction connected to JA
        tmpDist = getJunctionDistAA(tmpJA,JA);
        if(tmpDist >= 0) {
            AcA.push_back(tmpJA);
            dAcA.push_back((size_t)tmpDist);
        }
    }
    // collecting B junctions connected to JA/JB
    for(size_t i=0;i<NRB;++i) {
        tmpJB = i;
        // collecting B junction connected to JA
        if(JunctionDist[JA][tmpJB])
        BcA.push_back(tmpJB);
        // collecting B junction connected to JB
        tmpDist = getJunctionDistBB(tmpJB,JB);
        if(tmpDist >= 0) {
            BcB.push_back(tmpJB);
            dBcB.push_back((size_t)tmpDist);
        }
    }
    //AcB.push_back(JA); // not actually needed 
   // BcA.push_back(JB);
    /*cout<<"BcB[j]:";
    for(size_t j=0;j<BcB.size();++j) cout<<BcB[j]<<" ";
        cout<<"BcA[j]:";
    for(size_t j=0;j<BcA.size();++j) cout<<BcA[j]<<" ";
        cout<<endl;
        cout<<"AcA[j]:";
    for(size_t j=0;j<AcA.size();++j) cout<<AcA[j]<<" ";
        cout<<endl;
        cout<<"AcB[j]:";
    for(size_t j=0;j<AcB.size();++j) cout<<AcB[j]<<" ";
        cout<<endl;
    cout<<"JunctionDist[JA][JB] before update"<<JunctionDist[JA][JB]<<endl;*/
    
    /*if((AcB.size()==0 && BcA.size()!=0)|| (BcA.size()==0 && AcB.size()!=0)){
        cout<<"Problem in collect_connected-forward- AcB-BcA"<<endl;
        
        system("pause");
    //system("pause");
    }
    if((BcB.size()>0 && AcB.size()==0)){
        cout<<"Problem in collect_connected-forward- BcB-AcB"<<endl;
        
        system("pause");
    //system("pause");
    }*/
    
    //system("pause");
}


void UpdateJuncDist_forward(const size_t JA,const size_t JB,const vector<size_t> &BcA,const vector<size_t> &AcB)
{
    size_t olddist,newdist;
    double probAB;
    for(size_t i=0;i<BcA.size();++i) {
        tmpJB = BcA[i];
        for(size_t j=0;j<AcB.size();++j) {
            tmpJA = AcB[j];
            olddist = JunctionDist[tmpJA][tmpJB];
            newdist = JunctionDist[tmpJA][JB] + JunctionDist[JA][tmpJB] + 1;
            if(olddist == 0) { // i.e. previously unconnected
            if(newdist > MAX) newdist = MAX;
            JunctionDist[tmpJA][tmpJB] = newdist;
            //if distance is updated, sum of relative probabilities must also be updated
            probAB = (p[dist(newdist)] - p[dist(olddist)]) * (double)(fA-neighA[tmpJA].size())*(double)(fB-neighB[tmpJB].size());
            if( (neighA[tmpJA].size()<fA) && (neighB[tmpJB].size()<fB) ) {
                sum += probAB;
                sumA[tmpJA] += probAB;
            }
            continue;
            }
            if(newdist < olddist) {
                JunctionDist[tmpJA][tmpJB] = newdist;
                //if distance is updated, sum of relative probabilities must also be updated
                probAB = (p[dist(newdist)] - p[dist(olddist)]) * (double)(fA-neighA[tmpJA].size())*(double)(fB-neighB[tmpJB].size());
                if( (neighA[tmpJA].size()<fA) && (neighB[tmpJB].size()<fB) ) {
                    sum += probAB;
                    sumA[tmpJA] += probAB;
                }
                continue;
            }
        }
    }
}


void UpdateJuncDist_common(const size_t JA,const size_t JB,const vector<size_t> &AcA,const vector<size_t> &BcB,
const vector<size_t> &dAcA,const vector<size_t> &dBcB)
{
    size_t olddist,newdist;
    double probAB;
    for(size_t i=0;i<AcA.size();++i) {
        tmpJA = AcA[i];
        for(size_t j=0;j<BcB.size();++j) {
            tmpJB = BcB[j];
            olddist = JunctionDist[tmpJA][tmpJB];
            newdist = dAcA[i] + dBcB[j] + 1;
            if(olddist == 0) { // i.e. previously unconnected
                if(newdist > MAX) newdist = MAX;
                JunctionDist[tmpJA][tmpJB] = newdist;
                //if distance is updated, sum of relative probabilities must also be updated
                probAB = (p[dist(newdist)] - p[dist(olddist)]) * (double)(fA-neighA[tmpJA].size())*(double)(fB-neighB[tmpJB].size());
                if( (neighA[tmpJA].size()<fA) && (neighB[tmpJB].size()<fB) ) {
                    sum += probAB;
                    sumA[tmpJA] += probAB;
                }
                continue;
            }
            if(newdist < olddist) {
                JunctionDist[tmpJA][tmpJB] = newdist;
                //if distance is updated, sum of relative probabilities must also be updated
                probAB = (p[dist(newdist)] - p[dist(olddist)]) * (double)(fA-neighA[tmpJA].size())*(double)(fB-neighB[tmpJB].size());
                if( (neighA[tmpJA].size()<fA) && (neighB[tmpJB].size()<fB) ) {
                    sum += probAB;
                    sumA[tmpJA] += probAB;
                }
                continue;
            }
        }
    }
}


void UpdateMol_forward(const size_t JA,const size_t JB,const vector<size_t> &AcA,const vector<size_t> &AcB,
const vector<size_t> &BcA,const vector<size_t> &BcB)
{
    size_t mA,mB,newMol,delMol;
    mA = mol[JA];// mol is molecule index for each junction(monomer) ), so mA is molecule index of JA, and mB is molecule index of JB (since this is B, NRA+JB)
    mB = mol[NRA+JB];
    //cout<<"mA "<<mA<<", mB "<<mB<<"\n";

    // new molecule number should be the smaller of the two molecule numbers
    if(mA != mB) { // if mA==mB then this is intramolecular rxn, no change is needed
    newMol = (mA<mB)?mA:mB;// min(mA,mB)
    delMol = (mA<mB)?mB:mA;// max(mA,mB)
    //cout<<"newMol "<<newMol<<"\n";
    //cout<<"delMol "<<delMol<<"\n";
    for(size_t j=0;j<AcA.size();++j) mol[AcA[j]] = newMol;
    for(size_t j=0;j<BcA.size();++j) mol[NRA+BcA[j]] = newMol;
    for(size_t j=0;j<AcB.size();++j) mol[AcB[j]] = newMol;
    for(size_t j=0;j<BcB.size();++j) mol[NRA+BcB[j]] = newMol;
    //molSize[newMol] = molSize[mA] + molSize[mB];
    //molSize[delMol] = 0;
    mol[JA]=newMol;
    mol[JB+NRA]=newMol;
    //if(molSize[newMol]> molSize[largestMol]) largestMol = newMol;

    }
    //if(JA==328 && mA==606) {cout<<"forward";system("pause");}

}

int getJunctionDistAA(size_t J1,size_t J2)
{ // return -1 if J1 J2 are not connected, otherwise return minimal distance between them
// find minimum distance between J1 and neighbour of J2
    if(J1 == J2) return 0;
    if(neighA[J2].size()<1) return -1;
    if(JunctionDist[J1][neighA[J2][0]] == 0) return -1;
    size_t minD = MAX,JB; // Equivalent to: int minD = MAX; int JB;
    for(size_t i=0;i<neighA[J2].size();++i) {
            JB = neighA[J2][i];
            if(JunctionDist[J1][JB] < minD)
                minD = JunctionDist[J1][JB];
    }
    return minD+1;
}
int getJunctionDistBB(size_t J1,size_t J2)
{
    if(J1 == J2) return 0;
    if(neighB[J2].size()<1) return -1;
    if(JunctionDist[neighB[J2][0]][J1] == 0) return -1;
    size_t minD = MAX,JA;
    for(size_t i=0;i<neighB[J2].size();++i) {
        JA = neighB[J2][i];
        if(JunctionDist[JA][J1] < minD)
            minD = JunctionDist[JA][J1];
    }
    return minD+1;
}
size_t dist(unsigned char JunctionDist)
{
// JunctionDist should always be odd, since only different types of junctions react
    if(JunctionDist == 0) return 0;
    else return (JunctionDist+1)/(2*(dA+dB));
}
//Initializes random number generator with seed
//RNG is Mersenne Twister MT19937 algorithm
void RNG_initialize(unsigned long seed)
{
    mt[0]= seed & 0xffffffffUL;
    for(mti=1; mti<624; mti++){
    mt[mti] = (1812433253UL*(mt[mti-1]^(mt[mti-1] >> 30)) + mti);
    mt[mti] &= 0xffffffffUL;
    /* for >32 bit machines */
    }
    double dn = 3.442619855899;
    int i;
    const double m1 = 2147483648.0;
    double q;
    double tn = 3.442619855899;
    const double vn = 9.91256303526217E-03;
    q = vn/exp(-0.5*dn*dn);
    kn[0] = (dn/q)*m1;
    kn[1] = 0;
    wn[0] = q/m1;
    wn[127] = dn/m1;
    fn[0] = 1.0;
    fn[127] = exp(-0.5*dn*dn);
    for(i = 126; i >= 1; i--){
    dn = sqrt(-2*log(vn/dn + exp(-0.5*dn*dn)));
    kn[i+1] = (dn/tn)*m1;
    tn = dn;
    fn[i] = exp(-0.5*dn*dn);
    wn[i] = dn/m1;
    }
}
//Returns a random long between 0 and 4294967295
unsigned long rand32()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A for x=0,1 */
    if(mti >= 624){
    int kk;
    for(kk=0;kk<227;kk++){
    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
    mt[kk] = mt[kk+397] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for(;kk<623;kk++){
    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
    mt[kk] = mt[kk-227] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y = (mt[623]&UPPER_MASK)|(mt[0]&LOWER_MASK);
    mt[623] = mt[396] ^ (y >> 1) ^ mag01[y & 0x1UL];
    mti = 0;
    }
    y = mt[mti++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    //y=0.5*4294967296.0;
    return y;
}
