#pragma once

#include "lfsr.hpp"
#include <vector>

class geffe_generator{
private:
    lfsr r_30, r_31, r_32;
    bool clock_function(bool x, bool y, bool s);
public:
    void set_register(uint32_t reg, uint8_t index);
    bool clock();

    void generate_gamma(std::vector<uint8_t>& gamma, int size);

    geffe_generator();
};