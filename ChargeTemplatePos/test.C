#include "Minuit2/FCNBase.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include <vector>
#include "Minuit2/FCNGradAdapter.h"
#include <iostream>
#include "TROOT.h"
using namespace std;
class VertexRecLikelihoodFCN: public ROOT::Minuit2::FCNBase {
     public:
         VertexRecLikelihoodFCN() {}
         double operator() (const std::vector<double>& x) const{
             return a*x[0]*x[0];
         }
         double operator() (const double *x) const{
             std::vector<double> p(x,x+1);
             return (*this)(p);
         }
         double Up() const { return 0.5; }
     private:
         float a = 1;
};

VertexRecLikelihoodFCN* vtxllfcn = new VertexRecLikelihoodFCN();
ROOT::Minuit2::FCNBase* fMinuitFCN = new ROOT::Minuit2::FCNGradAdapter<ROOT::Minuit2::FCNGradAdapter<ROOT::Math::IMultiGradFunction>(vtxllfcn);
vector<double> a;
a.push_back(1.0);
vector<double> grad = fMinuitFCN->Gradient(a);
cout<<grad[0]<<endl;
