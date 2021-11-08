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
#include <iostream>

DECLARE_ALGORITHM(ChargeTemplateRec);

ChargeTemplateRec::ChargeTemplateRec(const std::string& name)
    : AlgBase(name),evt(0)
{
    tao_sipm = new TaoSiPM();
    charge_template = new ChargeTemplate();
    CD_radius = 900;
    
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
    evt->Branch("fSiPMHits", fSiPMHits, "fSiPMHits[4074]/I");
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
    evt->Branch("fDecayLength", &fDecayLength, "fDecayLength/f");
    evt->Branch("fChi2", &fChi2, "fChi2/f");
    evt->Branch("fEdm", &fEdm, "fEdm/f");
    evt->Branch("fScanX", fScanX, "fScanX[360]/F");
    evt->Branch("fScanXVal", fScanXVal, "fScanXVal[360]/F");
    evt->Branch("fScanY", fScanY, "fScanY[360]/F");
    evt->Branch("fScanYVal", fScanYVal, "fScanYVal[360]/F");
    evt->Branch("fScanZ", fScanZ, "fScanZ[360]/F");
    evt->Branch("fScanZVal", fScanZVal, "fScanZVal[360]/F");

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
    
    // Scan value
    // ScanLikelihood();
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

double ChargeTemplateRec::Chi2(
        double nhit,double vr,double vtheta,double vphi,double lambda)
{
    float total_chi2 = 0;
    float exp_dark_noise = tao_sipm->get_num() * tao_sipm->get_dark_noise_prob();

    // calculate some value that is needed.
    TVector3 v_vec = TVector3(0, 0 ,1);
    v_vec.SetMagThetaPhi(vr,vtheta,vphi);

    for(int i=0;i < tao_sipm->get_num(); i++)
    {

        float angle = v_vec.Angle(tao_sipm->get_vec(i) - v_vec);
        float exp_hit = CalExpChargeHit(vr, angle*180/PI, nhit, lambda);
        if(open_dark_noise){
            exp_hit += tao_sipm->get_dark_noise_prob();
        }else{
            exp_hit *= 1.0;
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

    ROOT::Math::Functor vtxllf(*vtxllfcn,5);
    vtxllminimizer->SetFunction(vtxllf);
    vtxllminimizer->SetMaxFunctionCalls(1e5);
    vtxllminimizer->SetMaxIterations(1e4);
    vtxllminimizer->SetTolerance(0.0001);
    vtxllminimizer->SetPrintLevel(1);
    
    TVector3 v_cc(fCCRecX,fCCRecY,fCCRecZ);
    float fCCRadius = v_cc.Mag();
    while (fCCRadius > 900)
    {
        v_cc *= (890./fCCRadius);
        fCCRadius = 890;
    } 
    
    float estimated_decay_length = 16.93*1000; // average absorption length in mm
    vtxllminimizer->SetVariable(0,"hits",fNSiPMHit*1.1,0.5);
    vtxllminimizer->SetVariable(1,"radius",v_cc.Mag(),0.01);
    vtxllminimizer->SetFixedVariable(2,"theta",v_cc.Theta());
    vtxllminimizer->SetFixedVariable(3,"phi",v_cc.Phi());
    vtxllminimizer->SetFixedVariable(4,"lambda",estimated_decay_length);

    int goodness = vtxllminimizer->Minimize();
    std::cout << "Vertex Minimize :: Goodness = " << goodness << std::endl;

    const double *xs = vtxllminimizer->X();
    TVector3 v_rec(0,0,1);
    v_rec.SetMagThetaPhi(xs[1],xs[2],xs[3]);
    fRecNHit = xs[0];
    fRecX    = v_rec.X();
    fRecY    = v_rec.Y();
    fRecZ    = v_rec.Z();
    fDecayLength    = xs[4];
    fChi2    = vtxllminimizer->MinValue(); 
    fEdm     = vtxllminimizer->Edm();
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
} 
 
double ChargeTemplateRec::LogPoisson(double obj,double exp_n)
{
    // simple chi2
    // double p = pow(exp_n - obj,2)/exp_n;
    
    // likelihood
    //double p = exp_n - obj * TMath::Log(exp_n);
    
    // likelihood ratio
    double p=2*(exp_n-obj);
    if(obj>0){
        p+=2*obj*TMath::Log(obj/exp_n);
    }
    return p;
}

float ChargeTemplateRec::CalExpChargeHit(float radius, float theta, float alpha, float lambda)
{
    float sipm_area = tao_sipm->get_sipm_area();
    float sipm_radius = tao_sipm->get_sipm_radius();
    float cos_theta = cos(theta*PI/180);
    float sin_theta = sin(theta*PI/180);
    float d = sqrt(sipm_radius*sipm_radius - radius*radius*sin_theta*sin_theta) - radius*cos_theta;
    float d_cd = sqrt(CD_radius*CD_radius - radius*radius*sin_theta*sin_theta) - radius*cos_theta;
    float cos_theta_proj = (d*d + sipm_radius*sipm_radius - radius*radius)/(2*d*sipm_radius);
    float exp_value = alpha*exp(-1.0*d_cd/lambda)*cos_theta_proj*sipm_area/(4*PI*d*d);
    return exp_value;
}

bool ChargeTemplateRec::ScanLikelihood()
{
    // Scan X first
    for(int i=0; i< NSCAN; i++){
        float x = (i - NSCAN/2)*900.0/(NSCAN/2);
        float y = fRecY; 
        float z = fRecZ;
        float nhit = fRecNHit;
        float lambda = fDecayLength; 
        float likelihood = Chi2(nhit,x,y,z,lambda);
        fScanX[i] = x;
        fScanXVal[i] = likelihood;
    }

    // Scan Y
    for(int i=0; i< NSCAN; i++){
        float y = (i - NSCAN/2)*900.0/(NSCAN/2);
        float x = fRecX; 
        float z = fRecZ;
        float nhit = fRecNHit;
        float lambda = fDecayLength; 
        float likelihood = Chi2(nhit,x,y,z,lambda);
        fScanY[i] = y;
        fScanYVal[i] = likelihood;
    }

    // Scan Z
    for(int i=0; i< NSCAN; i++){
        float z = (i - NSCAN/2)*900.0/(NSCAN/2);
        float y = fRecY; 
        float x = fRecX;
        float nhit = fRecNHit;
        float lambda = fDecayLength; 
        float likelihood = Chi2(nhit,x,y,z,lambda);
        fScanZ[i] = z;
        fScanZVal[i] = likelihood;
    }

    return true;
}
