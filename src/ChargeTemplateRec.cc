#include "ChargeTemplatePos/ChargeTemplateRec.h"
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

DECLARE_ALGORITHM(ChargeTemplateRec);

ChargeTemplateRec::ChargeTemplateRec(const std::string& name)
    : AlgBase(name),evt(0)
{
    tao_sipm = new TaoSiPM();
    charge_template = new ChargeTemplate();
    
    declProp("open_dark_noise",open_dark_noise = false);
    declProp("cc_factor",cc_factor = 0.7035);

    // print
    LogDebug << "Open dark noise : "<<open_dark_noise<<std::endl;
    LogDebug << "Charge center alg. factor : "<<cc_factor<<std::endl;
}

ChargeTemplateRec::~ChargeTemplateRec()
{
}

bool ChargeTemplateRec::initialize()
{
    // = access the geometry =
    // SniperPtr<SimGeomSvc> simgeom_svc(getParent(), "SimGeomSvc");
    // // == check exist or not ==
    // if (simgeom_svc.invalid()) {
    //     LogError << "can't find SimGeomSvc" << std::endl;
    //     return false;
    // }
    // == get the ROOT Geometry Manager ==
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

    evt = rootwriter->bookTree("ANASIMEVT/myevt", "user defined data");
    evt->Branch("evtID", &evtID, "evtID/I");
    evt->Branch("evtType", &evtType, "evtType/I");
    evt->Branch("fNSiPMHit", &fNSiPMHit, "fNSiPMHit/f");
    evt->Branch("fGdLSEdep", &fGdLSEdep, "fGdLSEdep/f");
    evt->Branch("fGdLSEdepX", &fGdLSEdepX, "fGdLSEdepX/f");
    evt->Branch("fGdLSEdepY", &fGdLSEdepY, "fGdLSEdepY/f");
    evt->Branch("fGdLSEdepZ", &fGdLSEdepZ, "fGdLSEdepZ/f");
    evt->Branch("fRecNHit", &fRecNHit, "fRecNHit/f");
    evt->Branch("fRecX", &fRecX, "fRecX/f");
    evt->Branch("fRecY", &fRecY, "fRecY/f");
    evt->Branch("fRecZ", &fRecZ, "fRecZ/f");
    evt->Branch("fCCRecX", &fCCRecX, "fCCRecX/f");
    evt->Branch("fCCRecY", &fCCRecY, "fCCRecY/f");
    evt->Branch("fCCRecZ", &fCCRecZ, "fCCRecZ/f");
    evt->Branch("fChi2", &fChi2, "fChi2/f");

    // create minimizer
    vtxllfcn = new VertexRecLikelihoodFCN(this);
    vtxllminimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");

    return true;
}

bool ChargeTemplateRec::execute()
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
    std::vector<int> fSiPMHitID = sim_event->SiPMHitID();
    for(int i=0;i<fSiPMHitID.size();i++)
    {
        fSiPMHits[fSiPMHitID[i]]++;
    }

    // Add dark noise simply
    if(open_dark_noise){
        float d_p = tao_sipm->get_dark_noise_prob();
        for(int i=0;i < tao_sipm->get_num(); i++)
        {
            float rand = myrandom();
            if(rand < d_p)
            {
                fSiPMHits[i] ++;
                fNSiPMHit ++;
            }
        }
    }

    // charge center reconstruction
    CalChargeCenter();

    // start reconstruction
    VertexMinimize();
    // fRecX = fGdLSEdepX;
    // fRecY = fGdLSEdepY;
    // fRecZ = fGdLSEdepZ;
    // fChi2 = Chi2(fNSiPMHit,fGdLSEdepX,fGdLSEdepY,fGdLSEdepZ);
    
    // fill the event.
    evt->Fill();

    //update here
    update();
    return true;
}

bool ChargeTemplateRec::finalize()
{
    tao_sipm->finalize();
    charge_template->finalize();
    return true;
}

double ChargeTemplateRec::Chi2(double nhit,double vx,double vy,double vz)
{
    float total_chi2 = 0;
    float exp_dark_noise = tao_sipm->get_num() * tao_sipm->get_dark_noise_prob();

    // calculate some value that is needed.
    TVector3 v_vec = TVector3(vx,vy,vz);
    // v_vec.SetMag(vx);
    // v_vec.SetTheta(vy);
    // v_vec.SetPhi(vz);
    float vr = v_vec.Mag();

    int b_template_index = charge_template->get_template_index(vr); 
    if (b_template_index >= TEMPLATENUM - 1){
        b_template_index = TEMPLATENUM - 2;
    }

    float b_tmp_radius = charge_template->get_template_radius(b_template_index);
    float f_tmp_radius = charge_template->get_template_radius(b_template_index + 1);
    float fx = abs(pow(f_tmp_radius,1) - pow(vr,1));
    float bx = abs(pow(b_tmp_radius,1) - pow(vr,1));
    float b_weight = fx/(fx + bx);
    float f_weight = bx/(fx + bx);
    for(int i=0;i < tao_sipm->get_num(); i++)
    {
        float d2 = (tao_sipm->get_vec(i) - v_vec).Mag2();
        float b_d2 = (tao_sipm->get_vec(i) - v_vec*(b_tmp_radius/vr)).Mag2();
        float f_d2 = (tao_sipm->get_vec(i) - v_vec*(b_tmp_radius/vr)).Mag2();

        float b_angle = v_vec.Angle(tao_sipm->get_vec(i) - v_vec*(b_tmp_radius/vr));
        float cos_b_angle = std::cos(b_angle);
        float f_angle = v_vec.Angle(tao_sipm->get_vec(i) - v_vec*(f_tmp_radius/vr));
        float cos_f_angle = std::cos(f_angle);
        float exp_hit = b_weight * (b_d2/d2) * charge_template->Interpolate(b_template_index,cos_b_angle) + f_weight * (f_d2/d2) * charge_template->Interpolate(b_template_index + 1, cos_f_angle);
        if(open_dark_noise){
            exp_hit *= (nhit - exp_dark_noise);
            exp_hit += tao_sipm->get_dark_noise_prob();
        }else{
            exp_hit *= nhit;
        }
        total_chi2 += LogPoisson(fSiPMHits[i],exp_hit);
    }
    // std::cout << total_chi2 << " " << vx << " " << vy << " " << vz <<" "<< vr <<std::endl;
    return total_chi2;
}

bool ChargeTemplateRec::update()
{
    for(int i=0;i < tao_sipm->get_num();i++)
    {
        fSiPMHits[i] = 0;
    }
    return true;
}

bool ChargeTemplateRec::VertexMinimize()
{

    ROOT::Math::Functor vtxllf(*vtxllfcn,4);
    vtxllminimizer->SetFunction(vtxllf);
    vtxllminimizer->SetMaxFunctionCalls(1e5);
    vtxllminimizer->SetMaxIterations(1e4);
    vtxllminimizer->SetTolerance(0.001);
    vtxllminimizer->SetPrintLevel(1);
    
    TVector3 v_cc(fCCRecX,fCCRecY,fCCRecZ);
    float fCCRadius = v_cc.Mag();
    while (fCCRadius > 900)
    {
        v_cc *= (890./fCCRadius);
        fCCRadius = 890;
    } 

    vtxllminimizer->SetLimitedVariable(0,"hits",fNSiPMHit,0.01,0.5*fNSiPMHit,1.5*fNSiPMHit);
    // vtxllminimizer->SetVariable(1,"radius",v_cc.Mag(),0.01);
    // vtxllminimizer->SetVariable(2,"theta",v_cc.Theta(),0.01);
    // vtxllminimizer->SetVariable(3,"phi",v_cc.Phi(),0.01);
    vtxllminimizer->SetVariable(1,"x",v_cc.X(),0.01);
    vtxllminimizer->SetVariable(2,"y",v_cc.Y(),0.01);
    vtxllminimizer->SetVariable(3,"z",v_cc.Z(),0.01);

    int goodness = vtxllminimizer->Minimize();
    std::cout << "Vertex Minimize :: Goodness = " << goodness << std::endl;

    const double *xs = vtxllminimizer->X();
    fRecNHit = xs[0];
    fRecX    = xs[1];
    fRecY    = xs[2];
    fRecZ    = xs[3];
    fChi2    = vtxllminimizer->MinValue(); 
    return true;
}

bool ChargeTemplateRec::CalChargeCenter()
{
    TVector3 cc_vec(0,0,0);
    for(int i=0; i < SIPMNUM;i ++){
        cc_vec += fSiPMHits[i]*tao_sipm->get_vec(i);
    }
    cc_vec *= (1.0/cc_factor)*(1.0/fNSiPMHit);
    if(open_dark_noise)
    {
        float exp_dark_noise = tao_sipm->get_num()*tao_sipm->get_dark_noise_prob();
        float corr_factor = fNSiPMHit/(fNSiPMHit - exp_dark_noise);
        cc_vec *= corr_factor;
    }

    fCCRecX = cc_vec.X();
    fCCRecY = cc_vec.Y();
    fCCRecZ = cc_vec.Z();
    // fCCRecX = fGdLSEdepX;    // fCCRecY = fGdLSEdepY;
    // fCCRecY = fGdLSEdepY;
    // fCCRecZ = fGdLSEdepZ;
} 
 
double ChargeTemplateRec::LogPoisson(double obj,double exp_n)
{
   

    /*
    double up=TMath::Poisson(obj,exp_n);
    double low=TMath::Poisson(obj,obj);
    if((-2.0*log(up/low))>1.e30) { return 1; }
    return -2.0*log(up/low);
    */

    //return pow(exp_n-obj,2.0)/exp_n;
    
     
    double p=2*(exp_n-obj);
    if(obj>0){
        p+=2*obj*TMath::Log(obj/exp_n);
    }
    return p;
}

float ChargeTemplateRec::CalExpChargeHit(float radius, float theta, float alpha)
{
    float sipm_area = tao_sipm->get_sipm_area();
    float sipm_radius = tao_sipm->get_sipm_radius();
    float cos_theta = cos(theta*PI/180);
    float sin_theta = sin(theta*PI/180);
    float d = sqrt(sipm_radius*sipm_radius - radius*radius*sin_theta*sin_theta) - radius*cos_theta;
    float cos_theta_proj = (d*d + sipm_radius*sipm_radius - radius*radius)/(2*d*sipm_radius);
    float exp_value = alpha*cos_theta_proj*sipm_area/(4*PI*d*d);
    return exp_value;
}
