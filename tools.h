class metrics {
public:
    int N;
    std::vector<double> x;
    double maxf, minf, meanf, varf;

    metrics(std::vector<double> x){
        this->x = x;
        N = x.size();
    }

    std::vector<double> getAll(void){
        getMax();
        getMin();
        getVar();
        std::vector<double> allf (1,meanf);
        allf.push_back(varf);
        allf.push_back(maxf);
        allf.push_back(minf);
        return allf;
    }

    double getMax(void){
        maxf = -1e9;
        for(int i=0;i<N;i++){
            if(x[i]>maxf){
                maxf = x[i];
            }
        }
        return maxf;
    }

    double getMin(void){
        minf = +1e9;
        for(int i=0;i<N;i++){
            if(x[i]<minf){
                minf = x[i];
            }
        }
        return minf;
    }

    double getMean(void){
        meanf = 0.0;
        for(int i=0;i<N;i++){
            meanf += x[i];
        }
        meanf /= (double)N;
        return meanf;
    }

    double getVar(void){
        getMean();
        varf = 0.0;
        for(int i=0;i<N;i++){
            varf += pow(meanf-x[i],2.0);
        }
        varf /= ((double)(N-1));
        return varf;
    }
};
