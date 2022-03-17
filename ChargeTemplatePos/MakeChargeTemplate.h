#ifndef MakeChargeTemplate_h
#define MakeChargeTemplate_h

#include "SniperKernel/AlgBase.h"
#include "SniperKernel/AlgFactory.h"
#include "ChargeTemplatePos/ChargeTemplate.h"
#include "ChargeTemplatePos/TaoSiPM.h"
#include "ChargeTemplatePos/ElecEffects.h"

#include "TH1F.h"
#include <vector>

#define NSCAN 360

/*
 * MakeChargeTemplate
 */

class MakeChargeTemplate : public AlgBase
{
    public:
        MakeChargeTemplate(const std::string & name);
        ~MakeChargeTemplate();
        
        bool initialize();
        bool execute();
        bool finalize();

        // update some state
        bool update();
        
    private:
        // input
        TaoSiPM* tao_sipm;
        ElecEffects* elec_effects;

        // result
        TH1F* charge_temp;
        float num_counted_events[360] = {0};
        int evtID;
        int evtType;
        float fGdLSEdep;
        float fGdLSEdepX;
        float fGdLSEdepY;
        float fGdLSEdepZ;
        std::vector<float> fPrimParticleX;
        std::vector<float> fPrimParticleY;
        std::vector<float> fPrimParticleZ;
        float fNSiPMHit;
        float fSiPMHits[SIPMNUM] = {0};
        std::vector<int> fSiPMHitID;

        // params to control the alg
        bool close_dark_noise;
        bool close_inter_ct;
        bool close_charge_resolution;
        float template_radius;
};

#endif
