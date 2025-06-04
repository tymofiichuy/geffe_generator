#include "lfsr.hpp"
#include <intrin.h>

using namespace std;

lfsr::lfsr(int len, uint32_t poly):length(len),polynom(poly),length_set(true){}

void lfsr::set_length(int len){
    if(len<1||len>64){
        throw invalid_argument("Invalid length");
    }
    length = len;
    length_set = true;
}

int lfsr::read_length(){
    if(length_set){
        return length;
    }
    else{
        throw runtime_error("Length isn't set");
    }
}

void lfsr::set_register(uint32_t reg){
    l_register = reg;
}

void lfsr::set_poly(uint32_t poly){
    polynom = poly;
}

bool lfsr::clock(){
    if(!length_set){
        throw runtime_error("Length isn't set");
    }
    uint32_t bit = 0x1;
    bool next = 0, res;
    for(int i = 0; i < length; i++){
        next ^= ((l_register&bit)!=0)&&((polynom&bit)!=0);
        bit <<= 1;
    }
    res = ((l_register&(0x1<<(length-1)))!=0);
    l_register <<= 1;
    l_register ^= static_cast<uint32_t>(next);
    return res;
}

bool lfsr::fast_clock(){
    uint32_t mask;
    bool res;
    res = ((l_register>>(length-1))&1);
    mask = l_register&polynom;
    l_register = (l_register<<1)^(__popcnt(mask)&0x1);
    return res;
}