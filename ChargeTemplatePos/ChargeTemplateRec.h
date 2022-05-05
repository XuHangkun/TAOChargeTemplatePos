#ifndef ChargeTemplateRec_h
#define ChargeTemplateRec_h

#undef _POSIX_C_SOURCE // to remove warning message of redefinition
#undef _XOPEN_SOURCE // to remove warning message of redefinition

#include "TTree.h"
#include "SniperKernel/AlgBase.h"
#include "SniperKernel/AlgFactory.h"
#include "ChargeTemplatePos/ChargeTemplate.h"
#include "ChargeTemplatePos/TaoSiPM.h"
#include "ChargeTemplatePos/ElecEffects.h"

#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Minuit2/FCNBase.h"
#include <vector>
#include <string>

#define NSCAN 360

/*
 * ChargeTemplateRec
 */

class ChargeTemplateRec : public AlgBase
{
    public:
        ChargeTemplateRec(const std::string & name);
        ~ChargeTemplateRec();      

        double Chi2(double nhit,double x,double y,double z, double alpha_ge68);
        double LogPoisson(double obj,double exp);
        void CorrectCCVertex();

        bool initialize();
        bool execute();
        bool finalize();

        // update some state
        bool update();
        
        // Charge Center Reconstruction
        bool CalChargeCenter();

        // Use Minimizer for vertex reconstruction
        bool VertexMinimize();
        class VertexRecLikelihoodFCN: public ROOT::Minuit2::FCNBase {
            public:
                VertexRecLikelihoodFCN(ChargeTemplateRec* rec_alg) { m_alg = rec_alg; }
                double operator() (const std::vector<double>& x) const{
                    return m_alg->Chi2(x[0],x[1],x[2],x[3],x[4]);
                }
                double operator() (const double *x) const{
                    std::vector<double> p(x,x+5);
                    return (*this)(p);
                }
                double Up() const { return 0.5; }
            private:
                ChargeTemplateRec* m_alg;
        };

        // charge expectation
        float CalExpChargeHit(float radius,float theta, float alpha, float alpha_ge68);

    private:
        // input
        TaoSiPM* tao_sipm;
        ElecEffects* elec_effects;
        ChargeTemplate* charge_template;
        ChargeTemplate* charge_template_ge68;
        float CD_radius;

        //minimizer
        ROOT::Math::Minimizer* vtxllminimizer;
        ROOT::Math::Minimizer* vtxllminimizer_migrad;
        VertexRecLikelihoodFCN* vtxllfcn; 
        
        // result
        TTree* evt;
        int evtID;
        int evtType;
        float fGdLSEdep;
        float fGdLSEdepX;
        float fGdLSEdepY;
        float fGdLSEdepZ;
        float fGdLSEdepR;
        float fGdLSEdepTheta;
        float fGdLSEdepPhi;
        float fNSiPMHit;
        float fSiPMHits[SIPMNUM] = {0};
        std::vector<int> fSiPMHitID;
        float fSiPMDN[SIPMNUM] = {0};
        float fSiPMCT[SIPMNUM] = {0};
        float fSiPMCR[SIPMNUM] = {0};
        float fRecNHit;
        float fRecX;
        float fRecY;
        float fRecZ;
        float fRecR;
        float fRecTheta;
        float fRecPhi;
        float fCCRecX;
        float fCCRecY;
        float fCCRecZ;
        float fCCRecR;
        float fCCRecTheta;
        float fCCRecPhi;
        float fRecGammaTempRatio;
        float fChi2;
        float fEdm;

        // params to control the alg
        char* charge_template_file;
        bool close_dark_noise;
        bool close_inter_ct;
        bool close_charge_resolution;
        float cc_factor;
};

#endif
