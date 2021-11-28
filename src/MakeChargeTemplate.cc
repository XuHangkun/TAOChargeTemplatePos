#include "ChargeTemplatePos/MakeChargeTemplate.h"
#include "ChargeTemplatePos/ChargeTemplate.h"
#include "ChargeTemplatePos/TaoSiPM.h"
#include "ChargeTemplatePos/Functions.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "Minuit2/FCNBase.h"
#include "TVector3.h"
#include <cmath>

#include "Event/SimHeader.h"
#include "Event/SimEvent.h"

#include "SniperKernel/AlgFactory.h"
#include "SniperKernel/SniperPtr.h"
#include "SniperKernel/SniperLog.h"
#include "RootWriter/RootWriter.h"
#include "BufferMemMgr/IDataMemMgr.h"
#include "EvtNavigator/NavBuffer.h"
#include "Geometry/SimGeomSvc.h"

#include <TGeoManager.h>

#include <boost/python.hpp>
#include <vector>
#include <iostream>

DECLARE_ALGORITHM(MakeChargeTemplate);

MakeChargeTemplate::MakeChargeTemplate(const std::string& name)
    : AlgBase(name)
{
    tao_sipm = new TaoSiPM();

    declProp("CloseDarkNoise", close_dark_noise = false);
    declProp("CloseInterCT", close_inter_ct = false);
    declProp("CloseChargeResolution", close_charge_resolution = false);
    declProp("TempRadius", template_radius = 0);

    // generate elec effects
    elec_effects = new ElecEffects(tao_sipm);


}

MakeChargeTemplate::~MakeChargeTemplate()
{
}

bool MakeChargeTemplate::initialize()
{
    // Parameters
    if ( close_dark_noise )
    { 
        elec_effects -> set_open_dark_noise(false); 
    } 
    if ( close_inter_ct )
    { 
        elec_effects -> set_open_inter_CT(false); 
    } 
    if ( close_charge_resolution )
    { 
        elec_effects -> set_open_charge_resolution(false); 
    } 
    std::cout << "Electronic Effects : " << std::endl;
    std::cout << "Open Dark Noise : "<<elec_effects->get_open_dark_noise() << std::endl;
    std::cout << "Open Internal Cross Noise : "<<elec_effects->get_open_inter_CT() << std::endl;
    std::cout << "Open Charge Resolution : "<<elec_effects->get_open_charge_resolution() << std::endl;

    // define charge template
    char temp_name[10];
    sprintf(temp_name, "r_%d", int(template_radius));
    charge_temp = new TH1F(temp_name, temp_name,360,0,180);

    // // = access the geometry =
    // SniperPtr<SimGeomSvc> simgeom_svc(getParent(), "SimGeomSvc");
    // // == check exist or not ==
    // if (simgeom_svc.invalid()) {
    //     LogError << "can't find SimGeomSvc" << std::endl;
    //     return false;
    // }
    // // == get the ROOT Geometry Manager ==
    // TGeoManager* geom = simgeom_svc->geom();
    // if (geom) {
    //     LogInfo << "get geometry geom: " << geom << std::endl;
    // }

    // = get RootWriter =
    SniperPtr<RootWriter> rootwriter(getParent(), "RootWriter");
    if (rootwriter.invalid()) {
        LogError << "Can't Find RootWriter. "
                 << std::endl;
        return false;
    }

    rootwriter->attach("ANASIMEVT/", charge_temp);

    return true;
}

bool MakeChargeTemplate::execute()
{

    SniperDataPtr<JM::NavBuffer> navBuf(getRoot(), "/Event");
    if (navBuf.invalid()) {
        return 0;
    }
    LogDebug << "navBuf: " << navBuf.data() << std::endl;

    JM::EvtNavigator* evt_nav = navBuf->curEvt();
    LogDebug << "evt_nav: " << evt_nav << std::endl;
    if (not evt_nav) {
        return 0;
    }

    Tao::SimHeader* sim_hdr = dynamic_cast<Tao::SimHeader*>(evt_nav->getHeader("/Event/Sim"));
    if (not sim_hdr) {
        return 0;
    }
    if (!sim_hdr->hasEvent()) {
        std::cout<<"no data is found, skip this event."<<std::endl;
        return true;
    }

    // == get event ==
    Tao::SimEvent* sim_event = dynamic_cast<Tao::SimEvent*>(sim_hdr->event());
    evtID = sim_event->EventID();
    evtType = sim_event->EventType();
    fGdLSEdep = sim_event->GdLSEdep();
    fGdLSEdepX = sim_event->GdLSEdepX();
    fGdLSEdepY = sim_event->GdLSEdepY();
    fGdLSEdepZ = sim_event->GdLSEdepZ();
    fNSiPMHit = sim_event->NSiPMHit();
    fSiPMHitID = sim_event->SiPMHitID();
    for(int i=0;i<fSiPMHitID.size();i++)
    {
        fSiPMHits[fSiPMHitID[i]]++;
    }

    // Add Electronic Effects
    fNSiPMHit = 0;
    for(int i=0;i < SIPMNUM;i++)
    {
        fSiPMHits[i] = elec_effects->AddElecEffects(fSiPMHits[i]);
        fNSiPMHit += fSiPMHits[i];
    }

    // Fill Histogram
    TVector3 edep_vec(fGdLSEdepX, fGdLSEdepY, fGdLSEdepZ);
    float edep_radius = edep_vec.Mag();
    if(edep_radius > template_radius - 1 && edep_radius < template_radius + 1)
    {
        for(int i=0; i<SIPMNUM;i++)
        {
            float theta = edep_vec.Angle(tao_sipm->get_vec(i) - edep_vec) * 180 / PI;
            int bin_index = charge_temp -> Fill(theta, fSiPMHits[i]*1.0/fNSiPMHit);
            num_counted_events[bin_index - 1] += 1;                        
        }
    }

    //update here
    update();
    return true;
}

bool MakeChargeTemplate::finalize()
{
    tao_sipm->finalize();
    for(int i=0; i< 360; i++)
    {
        if(num_counted_events[i] < 1)
        {
            continue;
        }
        float hits = charge_temp -> GetBinContent(i + 1);
        charge_temp -> SetBinContent(i + 1, hits/num_counted_events[i]);
        charge_temp -> SetBinError(i + 1, sqrt(hits/num_counted_events[i])/num_counted_events[i]);
    }
    return true;
}


bool MakeChargeTemplate::update()
{
    for(int i=0;i < tao_sipm->get_num();i++)
    {
        fSiPMHits[i] = 0;
    }
    return true;
}
