#include <morph/HdfData.h>
#include <morph/Config.h>
#include <morph/Random.h>
using namespace std;
#include "Quine.h" // #include <morph/bn/quine.h>
#include "tools.h"


// GLOBALS
morph::RandUniform<double, std::mt19937>* rng;
unsigned int nL; // number of ports/links

vector<int> randPermute(int);

class node {
public:
    int N; // number of possible input combinations
    std::vector<int> G;
    std::vector<int> C; // compatible binding?
    int Cr;
    int x;  // initial state
    node(void){
        x = (rng->get()<0.5);
        N = pow(2,nL); // note could do N=1<<nL
        for(int i=0;i<N;i++){
            G.push_back(rng->get()<0.5);
        }

        for(int i=0;i<nL;i++){
            C.push_back(rng->get()<0.5);
        }
        Cr = rng->get()<0.5;
    }
};


class network{
public:
    int n, nC;
    std::vector<int> X;
    std::vector<std::vector<int> > K;
    std::vector<int> nodeA, linkA, nodeB;
    std::vector<node*> nodes;
    std::vector<long unsigned int> QQ;
    double cycLen;

    network(std::vector<node*> nodes, int nC){

        this->nodes = nodes;
        n = nodes.size();
        if(nL*n<nC){
            std::cout<<"request fewer connections"<<std::endl;
            return;
        }
        X.resize(n,0);
        for(int i=0;i<nL;i++){
            QQ.push_back(pow(2,i));
        }

        // SETUP THE NETWORK CONNECTIONS
        std::vector<int> linkOrder = randPermute(n*nL);
        {
            int j=0;
            int k=0;
            while((k<nC)&&(j<10000)){
                int nA = floor(linkOrder[k]/nL);
                int lA = linkOrder[k]%nL;
                int nB = floor(rng->get()*n);

                // NOTE THE C/Cr PAIRING ENSURES THAT MEMBERS OF A SPECIES CANNOT INTERACT
                //if((nA!=nB) && (nodes[nA]->C[lA] != nodes[nB]->Cr)){

                if((nodes[nA]->Cr != nodes[nB]->Cr)){ // CANNOT INTERACT WITH THEMSELVES
                    nodeA.push_back(nA);
                    linkA.push_back(lA);
                    nodeB.push_back(nB);
                    k++;
                } else {
                    j++;
                }
            }
            this->nC = k;
        }
    }

    void initStates(void){
        for(int i=0;i<n;i++){
            X[i]=nodes[i]->x;
        }
    }

    void step(void){
        // ONE STEP OF THE DYNAMICS (SYNCHRONOUS UPDATING)
        std::vector<int> Xp = X;
        std::vector<long unsigned int> Gindex(n,0);
        for(int i=0;i<nC;i++){
            Gindex[nodeA[i]] += Xp[nodeB[i]] * QQ[linkA[i]];
        }
        for(int i=0;i<n;i++){
            X[i] = nodes[i]->G[Gindex[i]];
        }
    }


    std::vector<double> settle(void){

        // SETTLE NETWORK DYNAMICS AND IDENTIFY REPEATED STATE
        std::vector<std::vector<int> > V;
        std::vector<int> Vtar;
        bool rep = false;
        initStates();
        while(!rep){
            step();
            for(int i=0;i<V.size();i++){
                bool allsame = true;
                for(int j=0;j<X.size();j++){
                    if(V[i][j]!=X[j]){
                        allsame = false;
                        break;
                    }
                }
                if(allsame){
                    Vtar = V[i];
                    rep = true;
                    break;
                }
            }
            V.push_back(X);
        }

        // CYCLE ROUND THE ATTRACTOR ONCE MORE
        std::vector<double> Xcycle(X.size(),0.0);
        bool cycling = true;
        int cycleLen = 1;
        while(cycling){
            step();
            for(int j=0;j<X.size();j++){
                Xcycle[j] += X[j];
            }
            bool allsame = true;
            for(int j=0;j<X.size();j++){
                if(Vtar[j] != X[j]){
                    allsame = false;
                    break;
                }
            }
            if(allsame){
                cycling = false;
            } else {
                cycleLen++;
            }
        }
        cycLen = (double)cycleLen;
        for(int j=0;j<X.size();j++){
            Xcycle[j] /= cycLen;
        }
        return Xcycle;
    }


};




int main(int argc, char **argv){

    if (argc < 3) { std::cerr << "\nUsage: ./evodevo configfile logdir seed(optional)"; return -1; }

    srand(time(NULL)); // note may not be different for simultaneously launched progs
    int seed = rand();
    if(argc==4){
        seed = std::stoi(argv[3]);
    }
    morph::RandUniform<double, std::mt19937> _rng(seed);
    rng = &_rng;

    std::string paramsfile (argv[1]);
    morph::Config conf(paramsfile);
    if (!conf.ready) { std::cerr << "Error setting up JSON config: " << conf.emsg << std::endl; return 1; }

    std::string logpath = argv[2];
    std::ofstream logfile;
    morph::Tools::createDir (logpath);
    { std::stringstream ss; ss << logpath << "/log.txt"; logfile.open(ss.str());}
    logfile<<"Hello."<<std::endl;

    // GET PARAMS FROM CONFIG FILE
    int popSize = conf.getInt("popSize", 1000);
    int netSize = conf.getInt("netSize", 10);
    int numGens = conf.getInt("numGens", 1000);
    int netLinks = conf.getInt("totalLinks", 50);
    nL = conf.getInt("maxLinksPerNode", 6);
    int slowMeasurePeriod = conf.getInt("slowMeasurePeriod",numGens+2);
    double mutationRate = conf.getFloat("mutationRate", 0.0);
    int mutationPeriod = conf.getInt("mutationPeriod", numGens+2);
    bool inheretStartState = conf.getBool("inheretStartState", false);
    int speciesToCount = conf.getInt("speciesToCount", 10);

    // STORAGE/HELPER VECTORS
    std::vector<double> lower(popSize);
    std::vector<double> upper(popSize);
    std::vector<double> meanF, maxF, minF, varF;
    std::vector<double> AmeanF, AmaxF, AminF, AvarF;
    std::vector<double> BmeanF, BmaxF, BminF, BvarF;
    std::vector<double> CmeanF, CmaxF, CminF, CvarF;
    std::vector<double> DmeanF, DmaxF, DminF, DvarF;
    std::vector<int> numOfSpecies;
    std::vector<std::vector<int> > allSpeciesCounts(numGens,std::vector<int>(speciesToCount,0));

    // CREATE INITIAL POPULATION
    std::vector<node> P(popSize);

    int Glen = P[0].G.size();
    int Clen = P[0].C.size();
    int GlenHalf = Glen/2;
    int ClenHalf = Clen/2;

    // EVOLUTIONARY LOOP
    for(int g=0;g<numGens;g++){

        std::vector<double> score(popSize,0.0);
        std::vector<double> metricA(popSize,0.0);
        std::vector<double> metricB(popSize,0.0);
        std::vector<double> metricC(popSize,0.0);

        std::vector<double> accum(popSize,0);
        std::vector<double> H(popSize,0);
        std::vector<int> trials(popSize,0);
        std::vector<double> cycleLen(popSize,0);

        // LOOP OVER INDIVIDUALS IN POPULATION
        for(int t=0;t<popSize;t++){

            // MAKE RANDOM NET
            std::vector<bool> available(popSize,true);
            available[t]=false;
            std::vector<int> I(1,t);
            int s=1;
            while(s<netSize){
                int i = floor(rng->get()*popSize);
                if(available[i]){
                    I.push_back(i);
                    available[i]=false;
                    s++;
                }
            }
            std::vector<node*> S;
            for(int i=0;i<I.size();i++){
                S.push_back(&P[I[i]]);
                trials[I[i]]++;
            }

            // EVALUATE RANDOM NET
            network Net(S,netLinks);
            std::vector<double> Xcycle = Net.settle();

            // ACCUMULATE ENTROPY OF STATES VISITED OVER TRIALS
            for(int i=0;i<netSize;i++){
                double h=0.0;
                if((0.0<Xcycle[i])&&(Xcycle[i]<1.0)){
                    h = - (Xcycle[i]*log2(Xcycle[i])+(1.0-Xcycle[i])*log2(1.0-Xcycle[i]));
                }
                H[I[i]] += h;
            }

            // ACCUMULATE METRICS VALUES OVER TRIALS
            for(int i=0;i<netSize;i++){
                accum[I[i]] += Xcycle[i];
                cycleLen[I[i]] += Net.cycLen;
            }
        }

        // AVERAGE OVER TRIALS
        for(int i=0;i<popSize;i++){
            double meanH = (double)H[i]/(double)trials[i];
            score[i] = 1.0+meanH; // (the 1 is necessary!)
        }

        // EXTRACT METRICS
        for(int i=0;i<popSize;i++){
            // A=genome bias
            metricA[i] = 0.0;
            if(inheretStartState){
                metricA[i] = P[i].x;
            }
            for(int j=0;j<P[i].N;j++){
                metricA[i] += P[i].G[j];
            }
            if(inheretStartState){
                metricA[i] /= (double)(P[i].N+1);
            } else {
                metricA[i] /= (double)(P[i].N);
            }
            // B=phenome bias
            metricB[i] = (double)accum[i]/(double)trials[i];
            // C=cycle length
            metricC[i] = (double)cycleLen[i]/(double)trials[i];
        }

        // RECOMBINATION
        double cumSum = 0.0;
        for(int i=0;i<popSize;i++){
            lower[i]=cumSum;
            cumSum += score[i];
            upper[i] = cumSum;
        }

        std::vector<node> P2 = P;
        for(int i=0;i<popSize;i++){
            int mum;
            {
                double r = rng->get()*cumSum;
                for(int j=0;j<popSize;j++){
                    if((lower[j]<=r)&&(r<upper[j])){
                        mum=j;
                        break;
                    }
                }
            }
            int dad;
            {
                double r = rng->get()*cumSum;
                for(int j=0;j<popSize;j++){
                    if((lower[j]<=r)&&(r<upper[j])){
                        dad=j;
                        break;
                    }
                }
            }

            P2[i].x = (rng->get()<0.5);
            P2[i].G = P[mum].G;
            P2[i].C = P[mum].C;
            P2[i].Cr = P[mum].Cr;

            bool order = (rng->get()<0.5);
            if(order){ // dad first
                for(int k=0;k<GlenHalf;k++){
                    P2[i].G[k] = P[dad].G[k];
                }
                for(int k=0;k<ClenHalf;k++){
                    P2[i].C[k] = P[dad].G[k];
                }
                P2[i].Cr = P[dad].Cr;
            } else { // mum first
                for(int k=GlenHalf;k<Glen;k++){
                    P2[i].G[k] = P[dad].G[k];
                }
                for(int k=ClenHalf;k<Clen;k++){
                    P2[i].C[k] = P[dad].G[k];
                }
            }
        }

        // MUTATION
        if(!((g+1)%mutationPeriod)){
            for(int i=0;i<popSize;i++){
                for(int j=0;j<Glen;j++){
                    if(rng->get()<mutationRate){
                        P2[i].G[j] = 1-P[i].G[j];
                    }
                }
            }
        }

        // INSTATE NEW POPULATION
        P = P2;


        // IDENTIFY AND ORDER UNIQUE SPECIES
        std::vector<std::vector<int> > speciesIDsorted;
        std::vector<int> speciesCountsSorted;
        {
            std::vector<std::vector<int> > speciesID(1,P[0].G);
            std::vector<int> speciesCounts(1,1);
            for(int i=1;i<popSize;i++){
                bool different = true;
                for(int j=0;j<speciesID.size();j++){
                    bool identical = true;
                    for(int k=0;k<P[0].N;k++){
                        if(P[i].G[k] != speciesID[j][k]){
                            identical=false;
                            break;
                        }
                    }
                    if(identical){
                        different=false;
                        speciesCounts[j]++;
                        break;
                    }
                }
                if(different){
                    speciesID.push_back(P[i].G);
                    speciesCounts.push_back(1);
                }
            }
            numOfSpecies.push_back(speciesID.size());

            // SORT THE UNIQUE SPECIES AND COUNTS
            std::vector<std::vector<int> > sIDcp = speciesID;
            std::vector<int> sCTcp = speciesCounts;
            while(sIDcp.size()){
                std::vector<int> minx = sIDcp[0];
                int ct = sCTcp[0];
                int minid = 0;
                for(int i=1;i<sIDcp.size();i++){
                    bool iIsSmaller = true;
                    for(int j=0;j<sIDcp[i].size();j++){
                        if((sIDcp[i][j]!=minx[j]) && (sIDcp[i][j]==0)){
                            iIsSmaller=false;
                            break;
                        }
                    }
                    if(iIsSmaller){
                        minx = sIDcp[i];
                        minid = i;
                        ct = sCTcp[i];
                    }
                }
                speciesIDsorted.push_back(minx);
                sIDcp.erase(sIDcp.begin()+minid);
                speciesCountsSorted.push_back(ct);
                sCTcp.erase(sCTcp.begin()+minid);
            }
            if(speciesToCount<speciesCountsSorted.size()){
                for(int i=0;i<speciesToCount;i++){
                    allSpeciesCounts[g][i] = speciesCountsSorted[i];
                }
            } else {
                for(int i=0;i<speciesCountsSorted.size();i++){
                    allSpeciesCounts[g][i] = speciesCountsSorted[i];
                }
            }
        }

        // DETERMINE GENOME COMPLEXITY (Quine-McCluskey algorithm)
        std::vector<double> complexityUnique(speciesIDsorted.size(),0.0);
        if((g+1)%slowMeasurePeriod==0){
            for(int i=0;i<speciesIDsorted.size();i++){
                morph::bn::Quine Q(nL);
                for (unsigned int j=0; j<speciesIDsorted[i].size(); ++j) { // combs of other genes
                    if (speciesIDsorted[i][j]) { // if its a 1
                        Q.addMinterm(j);
                    }
                }
                Q.go();
                complexityUnique[i] = Q.complexity();
            }
        }
        std::vector<double> complexity;
        for(int i=0;i<speciesIDsorted.size();i++){
            for(int j=0;j<speciesCountsSorted[i];j++){
                complexity.push_back(complexityUnique[i]);
            }
        }


        // CALCULATE SUMMARY METRICS
        {
            metrics M(score);
            std::vector<double> m = M.getAll();
            meanF.push_back(m[0]);
            varF.push_back(m[1]);
            maxF.push_back(m[2]);
            minF.push_back(m[3]);
        }

        {
            metrics M(metricA);
            std::vector<double> m = M.getAll();
            AmeanF.push_back(m[0]);
            AvarF.push_back(m[1]);
            AmaxF.push_back(m[2]);
            AminF.push_back(m[3]);
        }

        {
            metrics M(metricB);
            std::vector<double> m = M.getAll();
            BmeanF.push_back(m[0]);
            BvarF.push_back(m[1]);
            BmaxF.push_back(m[2]);
            BminF.push_back(m[3]);
        }

        {
            metrics M(metricC);
            std::vector<double> m = M.getAll();
            CmeanF.push_back(m[0]);
            CvarF.push_back(m[1]);
            CmaxF.push_back(m[2]);
            CminF.push_back(m[3]);
        }

        {
            metrics M(complexity);
            std::vector<double> m = M.getAll();
            DmeanF.push_back(m[0]);
            DvarF.push_back(m[1]);
            DmaxF.push_back(m[2]);
            DminF.push_back(m[3]);
        }

        std::cout<<"  "<<g<<"     \r"<<std::flush;
    }


    // SAVE VALUES TO H5 OUTPUT FILE
    std::stringstream fname;
    fname << logpath << "/out.h5";
    morph::HdfData data(fname.str());
    std::stringstream path;

    // fitness
    path.str(""); path.clear(); path << "/meanF";
    data.add_contained_vals (path.str().c_str(), meanF);
    path.str(""); path.clear(); path << "/maxF";
    data.add_contained_vals (path.str().c_str(), maxF);
    path.str(""); path.clear(); path << "/minF";
    data.add_contained_vals (path.str().c_str(), minF);
    path.str(""); path.clear(); path << "/varF";
    data.add_contained_vals (path.str().c_str(), varF);

    // genome bias
    path.str(""); path.clear(); path << "/AmeanF";
    data.add_contained_vals (path.str().c_str(), AmeanF);
    path.str(""); path.clear(); path << "/AmaxF";
    data.add_contained_vals (path.str().c_str(), AmaxF);
    path.str(""); path.clear(); path << "/AminF";
    data.add_contained_vals (path.str().c_str(), AminF);
    path.str(""); path.clear(); path << "/AvarF";
    data.add_contained_vals (path.str().c_str(), AvarF);

    // phenome bias
    path.str(""); path.clear(); path << "/BmeanF";
    data.add_contained_vals (path.str().c_str(), BmeanF);
    path.str(""); path.clear(); path << "/BmaxF";
    data.add_contained_vals (path.str().c_str(), BmaxF);
    path.str(""); path.clear(); path << "/BminF";
    data.add_contained_vals (path.str().c_str(), BminF);
    path.str(""); path.clear(); path << "/BvarF";
    data.add_contained_vals (path.str().c_str(), BvarF);

    // cycle length
    path.str(""); path.clear(); path << "/CmeanF";
    data.add_contained_vals (path.str().c_str(), CmeanF);
    path.str(""); path.clear(); path << "/CmaxF";
    data.add_contained_vals (path.str().c_str(), CmaxF);
    path.str(""); path.clear(); path << "/CminF";
    data.add_contained_vals (path.str().c_str(), CminF);
    path.str(""); path.clear(); path << "/CvarF";
    data.add_contained_vals (path.str().c_str(), CvarF);

    // complexity
    path.str(""); path.clear(); path << "/DmeanF";
    data.add_contained_vals (path.str().c_str(), DmeanF);
    path.str(""); path.clear(); path << "/DmaxF";
    data.add_contained_vals (path.str().c_str(), DmaxF);
    path.str(""); path.clear(); path << "/DminF";
    data.add_contained_vals (path.str().c_str(), DminF);
    path.str(""); path.clear(); path << "/DvarF";
    data.add_contained_vals (path.str().c_str(), DvarF);

    // final population genome
    std::vector<int> GG;
    for(int i=0;i<popSize;i++){
        for(int j=0;j<P[0].N;j++){
            GG.push_back(P[i].G[j]);
        }
    }
    path.str(""); path.clear(); path << "/G";
    data.add_contained_vals (path.str().c_str(), GG);

    // species information
    path.str(""); path.clear(); path << "/species";
    data.add_contained_vals (path.str().c_str(), numOfSpecies);

    std::vector<int> SC;
    for(int i=0;i<allSpeciesCounts.size();i++){
        for(int j=0;j<speciesToCount;j++){
            SC.push_back(allSpeciesCounts[i][j]);
        }
    }
    path.str(""); path.clear(); path << "/SC";
    data.add_contained_vals (path.str().c_str(), SC);

    std::cout<<"Species Counts: "<<std::endl;
    for(int j=0;j<speciesToCount;j++){
        int ct = allSpeciesCounts[allSpeciesCounts.size()-1][j];
        if(ct>0){
            std::cout<<ct<<",";
        }
    }
    std::cout<<std::endl;

    return 0;
}


vector<int> randPermute(int x){
    vector<int> X(x);
    for(int i=0;i<x;i++){
        X[i] = i;
    }
    vector<int> Y(x);
    for(int i=0;i<x;i++){
        int index = floor(rng->get()*X.size());
        Y[i] = X[index];
        X.erase (X.begin()+index);
    }
    return Y;
}















/*
 for(int i=0;i<n;i++){
 std::vector<int> linkOrder = randPermute(n);
 int k=0;
 for(int j=0;j<n;j++){
 if(!(i==linkOrder[j])){
 nodeA.push_back(i);
 linkA.push_back(k);
 nodeB.push_back(linkOrder[j]);
 k++;
 }
 }
 }
 this->nC = n*(n-1);
 }
 */
