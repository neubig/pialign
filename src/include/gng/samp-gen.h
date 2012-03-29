#ifndef SAMPGEN_H__
#define SAMPGEN_H__

#include <cmath>
#include <cstdlib>
#include <vector>

inline int discreteSample(const double* vec, int size, double sum = -1) {
    if(size <= 0) throw std::runtime_error("Bad size-zero discrete sample");
    if(sum < 0) {
        sum = 0;
        for(int i = 0; i < size; i++)
            sum += vec[i];
    }
    sum *= (double)rand()/RAND_MAX;
    size -= 1; 
    int i = 0;
    for( ; i < size && ((sum -= vec[i]) > 0); i++);
    return i;
}

inline int discreteSample(const std::vector<double> & vec, double sum = -1) {
    return discreteSample(&vec[0],vec.size(),sum);
}

inline int discreteUniformSample(int size) {
    return rand() % size;
}

inline double addLogProbs(const std::vector<double> & probs) {
    const unsigned size = probs.size();
    double myMax = std::max(probs[0],probs[size-1]), norm=0;
    for(unsigned i = 0; i < probs.size(); i++)
        norm += exp(probs[i]-myMax);
    return log(norm)+myMax;
}
inline double addLogProbs(double a, double b) {
    double myMax = std::max(a,b);
#ifdef VITERBI_ON
    return myMax;
#else
    return log(exp(a-myMax)+exp(b-myMax))+myMax;
#endif
}

// distribution sampling functions
inline int bernoulliSample(double p) {
    return (rand() < p*RAND_MAX?1:0);
}
inline double exponSample(double l) {
    return -1*log(1-(double)rand()/RAND_MAX)/l;
}
inline double gammaSample(double a, double scale) {
    double b, c, e, u, v, w, y, x, z;
    if(a > 1) { // Best's XG method
        b = a-1;
        c = 3*a-.75;
        bool accept = false;
        do {
            u = (double)rand()/RAND_MAX;
            v = (double)rand()/RAND_MAX;
            w = u*(1-u);
            y = sqrt(c/w)*(u-.5);
            x = b+y;
            if(x >= 0) {
                z = 64*w*w*w*v*v;
                accept = (z <= 1-2*y*y/x || log(z) <= 2*(b*log(x/b)-y));
            }
        } while (!accept);
    } else { // Johnk's method
        do {
            u = (double)rand()/RAND_MAX;
            v = (double)rand()/RAND_MAX;
            x = pow(u,1/a);
            y = pow(v,1/(1-a));
        } while (x+y > 1);
        e = exponSample(1);
        x = e*x/(x+y);
    }
    return x / scale;
}
inline double betaSample(double a, double b) {
    double ga = gammaSample(a,1);
    double gb = gammaSample(b,1);
    return ga/(ga+gb);
}

inline double betaLogDensity(double x, double a, double b) {
    double ret =  log(tgamma(a+b)/tgamma(a)/tgamma(b)*pow(x,a-1)*pow(1-x,b-1));
    return ret;
}
inline double gammaLogDensity(double x, double k, double t) {
    double ret = log(pow(x,k-1)*exp(-1*x/t)/pow(t,k)/tgamma(k));
    return ret;
}

#endif
