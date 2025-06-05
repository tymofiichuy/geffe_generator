#pragma once

#include "geffe_generator.hpp"

#include <string>
#include <bitset>
#include <ppl.h>
#include <concurrent_vector.h>

class register_recovery{
private:
    float alpha_quantile;
    float beta_quantile;
    float criterion;
    int population;
    bool parameters_set = false;

    std::bitset<320> gamma_template;
    concurrency::concurrent_vector<uint32_t> L1_candidates;
    concurrency::concurrent_vector<uint32_t> L2_candidates;
public:
    void set_quantiles(float alpha_q, float beta_q);
    void set_critical_set();
    void set_gamma_template(std::string& sequence);

    // void prepare_file(const std::string& out_file);
    bool recognize(std::bitset<320>& gamma);
    void recover_L1();
    void recover_L2();
    void recover_L3();

    register_recovery(float alpha_q = 0, float beta_q = 0);
};