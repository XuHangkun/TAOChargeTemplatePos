#include "ChargeTemplatePos/ElecEffects.h"
#include "ChargeTemplatePos/Functions.h"
#include <cmath>
#include <random>
#include <iostream>
using namespace std;

ElecEffects::ElecEffects(TaoSiPM* a_tao_sipm)
{
    tao_sipm = a_tao_sipm;
    
    // some switch
    open_dark_noise = true;
    open_inter_CT = true;
    open_charge_resolution = true;

    // initialize
    Initialize();
}

ElecEffects::~ElecEffects()
{

}

bool ElecEffects::Initialize()
{
    for(int i=0;i<APPROX_NUM;i++)
    {
        if(i == 0){
            dark_noise_ps[i] += poisson(tao_sipm->get_dark_noise_prob(),i);
            cross_talk_ps[i] += ct_prob(tao_sipm->get_inter_cross_talk_prob(),i + 1);
        }else{
            dark_noise_ps[i] = dark_noise_ps[i-1] + poisson(tao_sipm->get_dark_noise_prob(),i);
            cross_talk_ps[i] = cross_talk_ps[i-1] + ct_prob(tao_sipm->get_inter_cross_talk_prob(),i + 1);
        }
    }
    
    // print something
    cout << "Dark Noise Prob. " << dark_noise_ps[0] << "\t" << dark_noise_ps[1] << "\t" << dark_noise_ps[2] << endl;
    cout << "Cross Talk Prob. " << cross_talk_ps[0] << "\t" << cross_talk_ps[1] << "\t" << cross_talk_ps[2] << endl;

    return true;
}

float ElecEffects::AddDarkNoise(float hit)
{
    if(!open_dark_noise)
    {
        return 0;
    }
    // number of dark noise satisfy poisson distribution
    // m = 0.15
    // p(n|0.15) = (m)^n * exp(-n) / (n!)
    // we will use approximation, only consider 0,1,2,3
    float random = myrandom();
    int dk_hit = 0;
    for(int i=0;i < APPROX_NUM; i++)
    {
        if(i == APPROX_NUM - 1 && random > dark_noise_ps[i])
        {
            dk_hit = i + 1;
            break;
        }
        if(random < dark_noise_ps[i])
        {
            dk_hit = i;
            break;
        }
    }
    return dk_hit;
}

float ElecEffects::AddInterCT(float hit)
{
    if(!open_inter_CT)
    {
        return 0;
    }
    // number of dark noise satisfy poisson distribution
    // we will use approximation, only consider 0,1,2,3
    float total_ct_hit = 0;
    for(int j = 0; j < hit; j++)
    {
        float random = myrandom();
        int ct_hit = 0;
        for(int i=0;i < APPROX_NUM; i++)
        {
            if(i == APPROX_NUM - 1 && random > cross_talk_ps[i])
            {
                ct_hit = i + 1;
                break;
            }
            if(random < cross_talk_ps[i])
            {
                ct_hit = i;
                break;
            }
        }
        total_ct_hit += ct_hit;
    }
    return total_ct_hit;
}

float ElecEffects::AddChargeResolution(float hit)
{
    if(!open_charge_resolution)
    {
        return 0;
    }
    if(hit < 0.1)
    {
        return 0;
    }
    float total_cr_hit = 0;
    std::normal_distribution<double> charge_resolution_dist(0,1*tao_sipm->get_charge_resolution());
    total_cr_hit = hit*charge_resolution_dist(generator);
    return total_cr_hit;
}

float ElecEffects::AddElecEffects(float hit)
{
    float total_hit = hit;
    float ct_hit = 0;
    float dk_hit = 0;
    float cr_hit = 0;
    if(open_inter_CT)
    {
        // Add Cross talk
        ct_hit = AddInterCT(hit);
        total_hit += ct_hit;
    }
    if(open_dark_noise)
    {
        // Add Dark Noise
        dk_hit = AddDarkNoise(hit);
        total_hit += dk_hit;
    }
    if(open_charge_resolution)
    {
        // Add Charge Resolution
        cr_hit = AddChargeResolution(total_hit);
        total_hit += cr_hit;
    }
    // cout << "true hit : " << hit << endl;
    // cout << "cross talk hit : " << ct_hit << endl;
    // cout << "dark noise hit : " << dk_hit << endl;
    // cout << "charge resolution hit : " << cr_hit << endl;
    return total_hit;
}

float ElecEffects::poisson(float m, int n)
{
    float value = pow(m,n)*exp(-1.0*m);
    for(int i=1;i<=n;i++)
    {
        value /= (i*1.0);
    }
    return value;
}

float ElecEffects::ct_prob(float c, int n)
{
    float m = c*n;
    float value = pow(m,n - 1)*exp(-1.0*m);
    for(int i=1;i<=n;i++)
    {
        value /= (i*1.0);
    }
    return value;
}
