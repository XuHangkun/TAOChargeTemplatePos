#ifndef TaoSiPM_h
#define TaoSiPM_h

#define SIPMNUM 4074
#include "TVector3.h"

/*
 * TaoSiPM
 */

class TaoSiPM
{
    public:
        TaoSiPM();
        ~TaoSiPM();

        bool initialize();
        bool finalize();

        int get_num() { return sipm_num; }
        float get_sipm_area() { return sipm_area; }
        float get_sipm_radius() { return sipm_radius; }
        float get_dark_noise_prob() { return dark_noise_prob; }
        float get_theta(const int index) { return sipm_theta[index]; }
        float get_phi(const int index) { return sipm_phi[index]; }
        TVector3 get_vec(const int index) { return sipm_vec[index]; }

    private:
        //params of sipm
        int sipm_num;
        float sipm_area;    // mm^2
        float sipm_noise;   // Hz/mm^2
        float sipm_readout_window; //s
        float dark_noise_prob;
        float sipm_radius;
        float sipm_theta[SIPMNUM] = {0};
        float sipm_phi[SIPMNUM] = {0};
        TVector3 sipm_vec[SIPMNUM];
};

#endif
