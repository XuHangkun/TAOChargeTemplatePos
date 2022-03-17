#include "ChargeTemplatePos/TaoSiPM.h"
#include <fstream>
#include <iostream>
#include <string>
#include "TVector3.h"
#include <cmath>

using namespace std;

TaoSiPM::TaoSiPM()
{
    sipm_num = SIPMNUM;
    sipm_radius = 930.2;
    sipm_area = 2500;
    sipm_noise = 15;
    sipm_readout_window = 600*1.e-9;
    dark_noise_prob = sipm_area*sipm_noise*sipm_readout_window;
    inter_cross_talk_prob = 0.05;
    charge_resolution = 0.16;
    initialize();
}

TaoSiPM::~TaoSiPM()
{

}

bool TaoSiPM::initialize()
{
    string root_dir = getenv("CHARGETEMPLATEPOSROOT");
    string file = "/input/sipm_pos.txt";
    ifstream infile(root_dir + file);
    int index = 0;
    float theta,phi;
    for(int i = 0;i < SIPMNUM;i++)
    {
        infile >> index >> theta >> phi;
        // cout << "SiPM index : "<< index <<" theta : "<< theta <<" phi : "<<phi<<endl;
        sipm_theta[index] = theta;
        sipm_phi[index] = phi;
        sipm_vec[index] = TVector3(
                sipm_radius*sin(theta*3.141592/180)*cos(phi*3.141592/180),
                sipm_radius*sin(theta*3.141592/180)*sin(phi*3.141592/180),
                sipm_radius*cos(theta*3.141592/180)
                );
    }
    infile.close();
    cout << "TaoSiPM Initialization Finished !!!" << endl;
    return true;
}

bool TaoSiPM::finalize()
{
    return true;
}
