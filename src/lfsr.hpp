#pragma once

#include <iostream>

class lfsr{
private:
    uint32_t l_register = 0;
    uint32_t polynom;
    int length;
    bool length_set = false;
public:
    void set_length(int len);
    int read_length();

    void set_register(uint32_t reg);
    void set_poly(uint32_t poly);
    bool clock();
    bool fast_clock();

    lfsr(int len = 0, uint32_t poly = 0);
};