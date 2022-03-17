#include "ChargeTemplatePos/ChargeTemplate.h"
#include "TFile.h"
#include "math.h"
#include "TString.h"
#include "TVector3.h"
#include <iostream>
#include <fstream>
#include <string>
#include "ChargeTemplatePos/Functions.h"

using namespace std;

ChargeTemplate::ChargeTemplate(string charge_tmp_file)
{
    tmp_num = TEMPLATENUM;
    cd_radius = 900;
    sipm_radius = 930.2;
    max_tmp_radius = 0;
    tmp_numbers = 0;
    charge_template_file = charge_tmp_file;
    initialize();
}

ChargeTemplate::~ChargeTemplate()
{
}

bool ChargeTemplate::initialize()
{
    string root_dir = getenv("CHARGETEMPLATEPOSROOT");
    ifstream info_file((root_dir + "/input/" + charge_template_file + ".txt").c_str());
    string file_name = "/input/" + charge_template_file +".root";
    tmp_file = new TFile((root_dir + file_name).c_str());

    double radius = 0;
    char tmp_name[30];
    while(info_file >> radius >> tmp_name)
    { 
        tmp_radius.push_back(radius);
        TH1F* hist = (TH1F*) tmp_file->Get(tmp_name);
        for(int i = 0;i <180; i++)
        {
            hist -> Interpolate(i);
        }
        tmp.push_back(hist);
        cout << "Templare radius : "<< radius <<" name : "<< tmp_name <<endl;
        tmp_numbers ++;
        max_tmp_radius = radius;
    }
    cout << "Max template radius : " << max_tmp_radius << "\nTemplate numbers : "<< tmp_numbers << endl;
    info_file.close();

    cout << "ChargeTemplate Initialization Finished !!!" << endl;
    return true;
}

bool ChargeTemplate::finalize()
{
    tmp_file->Close();
    return true;
}

TH1F* ChargeTemplate::get_template(int index)
{
    return tmp[index];
}

float ChargeTemplate::get_template_radius(int index)
{
    return tmp_radius[index];
}

float ChargeTemplate::cal_sipm_proj(float radius, float sipm_distance)
{
    float cos_theta_proj = (sipm_distance*sipm_distance + sipm_radius*sipm_radius - radius*radius)/(2*sipm_distance*sipm_radius);
    return cos_theta_proj;
}

float ChargeTemplate::cal_sipm_distance(float radius, float theta)
{
    float cos_theta = cos(theta*PI/180);
    float sin_theta = sin(theta*PI/180);
    float d = sqrt(sipm_radius*sipm_radius - radius*radius*sin_theta*sin_theta) - radius*cos_theta;
    return d;
}

float ChargeTemplate::LinearInterpolation(float radius, float x0, float y0, float x1, float y1)
{
    float value = 0;
    if (fabs(x0 - x1) < 1.e-2)
    {
        value = (y0 + y1)/2.0;
    }else{
        value = y0 + (radius - x0)*(y1 - y0)/(x1 - x0);
    }
    return value;
}

int ChargeTemplate::FindBeforeIndex(float radius, int low, int high)
{
    // for(int i = 0; i < tmp_numbers - 1; i++)
    // {
    //     if((tmp_radius[i] <= radius) & (radius < tmp_radius[i+1]))
    //     {
    //         return i;
    //     }
    // }
    // return 0;
    if (low == high){
        if(radius < tmp_radius[low]){
            return max(low - 1, 0);
        }else{
            return low;
        }
    }else if(high < low){
        return max(low - 1, 0);
    }
    int mid = int((low + high)/2);
    if(radius >= tmp_radius[mid]){
        return FindBeforeIndex(radius, mid + 1, high);
    }else{
        return FindBeforeIndex(radius, low, mid - 1);
    }

}

float ChargeTemplate::CalExpChargeHit(float radius, float theta)
{

    int bindex = FindBeforeIndex(radius,0,tmp_numbers-1);
    int findex = bindex + 1;
    if(radius >= max_tmp_radius)
    {
        findex = tmp_numbers - 1;
        bindex = findex - 1;
    }
    float cos_theta = cos(theta*PI/180);
    float sin_theta = sin(theta*PI/180);
    float sipm_distance = cal_sipm_distance(radius, theta);

    // before charge template information
    float b_tmp_radius = get_template_radius(bindex);
    float b_sipm_distance = cal_sipm_distance(b_tmp_radius, theta);
    float b_correct_factor = cal_sipm_proj(b_tmp_radius, b_sipm_distance)/pow(b_sipm_distance,2);
    TH1F* b_temp = get_template(bindex);
    float b_temp_hit = b_temp -> Interpolate(theta)/b_correct_factor;

    // after charge template information
    float f_tmp_radius = get_template_radius(findex);
    float f_sipm_distance = cal_sipm_distance(f_tmp_radius, theta);
    float f_correct_factor = cal_sipm_proj(f_tmp_radius, f_sipm_distance)/pow(f_sipm_distance,2);
    TH1F* f_temp = get_template(findex);
    float f_temp_hit = f_temp -> Interpolate(theta) / f_correct_factor;
    
    // get linear interpolation
    float exp_hit = LinearInterpolation(pow(radius,1), pow(b_tmp_radius,1), b_temp_hit, pow(f_tmp_radius,1), f_temp_hit);
    
    // cal ..
    float correct_factor = cal_sipm_proj(radius, sipm_distance)/pow(sipm_distance,2);
    exp_hit = correct_factor * exp_hit;
    return exp_hit;
}
