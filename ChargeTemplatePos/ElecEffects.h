#ifndef ELECEFFECTS_h
#define ELECEFFECTS_h
#include <random>
#include "TaoSiPM.h"
#define APPROX_NUM 5

class ElecEffects
{

    public:
        ElecEffects(TaoSiPM* tao_sipm);
        ~ElecEffects();

        bool Initialize();
        
        float AddDarkNoise(float hit);
        float AddInterCT(float hit);
        float AddChargeResolution(float hit);
        float AddElecEffects(float hit);

        // usefull function
        float poisson(float m, int n);
        float ct_prob(float c, int n);

        // get and set
        void set_open_dark_noise(bool v) { open_dark_noise = v; }
        bool get_open_dark_noise() { return open_dark_noise; }
        void set_open_inter_CT(bool v) { open_inter_CT = v; }
        bool get_open_inter_CT() { return open_inter_CT; }
        void set_open_charge_resolution(bool v) { open_charge_resolution = v; }
        bool get_open_charge_resolution() { return open_charge_resolution; }

    private:
        TaoSiPM* tao_sipm;
        bool open_dark_noise;
        bool open_inter_CT;
        bool open_charge_resolution;

        // accumulate param
        float dark_noise_ps[APPROX_NUM] = {0};
        float cross_talk_ps[APPROX_NUM] = {0};
        // charge resolution
        std::default_random_engine generator;
};

#endif
