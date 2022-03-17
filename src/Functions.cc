#include "ChargeTemplatePos/Functions.h"
#include <cstdlib>

using namespace std;

float myrandom()
{
    return rand()%(100000)*1.0/100000;
}

int max(int x, int y)
{
    if(x >= y)
    {
        return x;
    }
    return y;
}
