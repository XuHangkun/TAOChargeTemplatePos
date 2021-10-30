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
    string file_name = "/input/templates.root";
    tmp_file = new TFile((root_dir + file_name).c_str());
    for(int i = 0;i < tmp_num ;i++)
    {
        double radius = 0;
        char tmp_name[30];
        info_file >> radius >> tmp_name;   
        cout << radius<< " " << tmp_name << endl;
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

float ChargeTemplate::Interpolate(int index, float cos_theta)
{
    float angle = acos(cos_theta)*180/PI; 
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
