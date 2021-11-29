#include "ChargeTemplatePos/ChargeTemplate.h"
#include "TFile.h"
#include "math.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include "ChargeTemplatePos/Functions.h"

using namespace std;

ChargeTemplate::ChargeTemplate()
{
    tmp_num = TEMPLATENUM;
    cd_radius = 900;
    sipm_radius = 900;
    delta_r3 = pow(cd_radius,3)*1.0/(tmp_num - 1);
    initialize();
}

ChargeTemplate::~ChargeTemplate()
{
}

bool ChargeTemplate::initialize()
{
    string root_dir = getenv("CHARGETEMPLATEPOSROOT");
    ifstream info_file((root_dir + "/input/templates.txt").c_str());
    string file_name = "/input/old_templates.root";
    tmp_file = new TFile((root_dir + file_name).c_str());
    for(int i = 0;i < tmp_num ;i++)
    {
        double radius = 0;
        char tmp_name[30];
        info_file >> radius >> tmp_name;   
        tmp_radius[i] = radius;
        // sprintf(tmp_name,"r_%d.0",int(radius));
        tmp[i] = (TH1F*) tmp_file->Get(tmp_name);
        // cout<< tmp_name << tmp[i]->Interpolate(90) << endl;
    }
    info_file.close();
    cout << "ChargeTemplate Initialization Finished !!!" << endl;
    return true;
}

bool ChargeTemplate::finalize()
{
    tmp_file->Close();
    return true;
}

int ChargeTemplate::get_template_index(float radius)
{
    for(int i=0; i < tmp_num - 1; i++)
    {
        if(tmp_radius[i] <= radius and radius < tmp_radius[i+1])
        {
            return i;
        }
    }
    return tmp_num - 1;
}

TH1F* ChargeTemplate::get_template(int index)
{
    int c_index = correct_index(index);
    return tmp[c_index];
}

float ChargeTemplate::Interpolate(int index, float theta)
{
    float angle = theta; 
    int c_index = correct_index(index);
    return tmp[c_index]->Interpolate(angle);
}

float ChargeTemplate::get_template_radius(int index)
{
    int c_index = correct_index(index);
    return tmp_radius[c_index];
}
int ChargeTemplate::correct_index(int index)
{
    if(index < 0)
    {
        return 0;
    }else if(index >= TEMPLATENUM)
    {
        return TEMPLATENUM - 1;
    }
    return index;
}

float ChargeTemplate::cal_sipm_proj(float radius, float theta)
{
    float cos_theta = cos(theta*PI/180);
    float sin_theta = sin(theta*PI/180);
    float d = sqrt(sipm_radius*sipm_radius + radius*radius*sin_theta*sin_theta) - radius*cos_theta;
    float cos_theta_proj = (d*d + sipm_radius*sipm_radius - radius*radius)/(2*d*sipm_radius);
    return cos_theta_proj;
}

float ChargeTemplate::cal_sipm_distance(float radius, float theta)
{
    float cos_theta = cos(theta*PI/180);
    float sin_theta = sin(theta*PI/180);
    float d = sqrt(sipm_radius*sipm_radius + radius*radius*sin_theta*sin_theta) - radius*cos_theta;
    return d;
}

float ChargeTemplate::CalExpChargeHit(float radius, float theta, float alpha)
{
    int index = get_template_index(radius);
    // before charge template information
    float b_tmp_radius = get_template_radius(index);
    float b_correct_factor = cal_sipm_proj(b_tmp_radius, theta)/pow(cal_sipm_distance(b_tmp_radius, theta),2);
    TH1F* b_temp = get_template(index);
    float b_temp_hit = b_temp -> Interpolate(theta);
    // after charge template information
    float f_tmp_radius = get_template_radius(index + 1);
    float f_correct_factor = cal_sipm_proj(f_tmp_radius, theta)/pow(cal_sipm_distance(f_tmp_radius, theta),2);
    TH1F* f_temp = get_template(index + 1);
    float f_temp_hit = f_temp -> Interpolate(theta);

    // calculate weight
    float b_weight = abs(pow(radius,3) - pow(f_tmp_radius,3))/(abs(pow(radius,3) - pow(b_tmp_radius,3)) + abs(pow(radius,3) - pow(f_tmp_radius,3)));

    // cal ..
    float correct_factor = cal_sipm_proj(radius, theta)/pow(cal_sipm_distance(radius, theta),2);
    float exp_hit = correct_factor*(b_weight*b_temp_hit/b_correct_factor + (1 - b_weight)*f_temp_hit/f_correct_factor);
    
    return alpha*exp_hit;
}
