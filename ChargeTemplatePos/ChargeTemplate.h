#ifndef ChargeTemplate_h
#define ChargeTemplate_h

#include "TH1F.h"
#include "TFile.h"
#define TEMPLATENUM 101

/*
 * ChargeTemplate
 */

class ChargeTemplate
{
    public:
        ChargeTemplate();
        ~ChargeTemplate();

        bool initialize();
        bool finalize();

        int correct_index(int index);

        int get_template_index(float radius);
        float get_template_radius(int index);
        TH1F* get_template(int index);
        float Interpolate(int index,float theta);

        float CalExpChargeHit(float radius, float theta, float alpha);
        float cal_sipm_proj(float radius, float theta);
        float cal_sipm_distance(float radius, float theta);

    private:
        int tmp_num;
        float cd_radius;
        float sipm_radius;
        float tmp_radius[TEMPLATENUM];
        TH1F* tmp[TEMPLATENUM];
        TFile* tmp_file;
};

#endif
