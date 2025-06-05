#include "geffe_generator.hpp"

using namespace std;

geffe_generator::geffe_generator():r_30(30, 0x32800000),r_31(31, 0x48000000),r_32(32, 0xF5000000){};

bool geffe_generator::clock_function(bool x, bool y, bool s){
    return ((s&&x)^((!s)&&y));
}

void geffe_generator::set_register(uint32_t reg, uint8_t index){
    switch(index)
    {
    case 0:
        if((reg>>30)!=0){
            throw invalid_argument("Invalid register");
        }
        r_30.set_register(reg);
        break;
    case 1:
        if((reg>>31)!=0){
            throw invalid_argument("Invalid register");
        }
        r_31.set_register(reg);
        break;
    case 2:
        r_32.set_register(reg);
        break;
    default:
        throw invalid_argument("Invalid index");
        break;
    }
}

bool geffe_generator::clock(){
    bool x = r_30.fast_clock(), y = r_31.fast_clock(), s = r_32.fast_clock();
    return clock_function(x, y, s);
}

void geffe_generator::generate_gamma(vector<uint8_t>& gamma, int size){
    if(gamma.size()!=size){
        gamma.resize(size);
    }
    for(int i = 0; i < size; i++){
        gamma[i] = static_cast<uint8_t>(clock());
    }
}