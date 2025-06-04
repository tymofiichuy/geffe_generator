#pragma once

#include <string>

template<typename T> class register_recovery{
private:
    float alpha_quantile;
    float beta_quantile;
public:
    void set_quantiles(float alpha_q, float beta_q);
    void prepare_file(const std::string& out_file);
    bool recognize();
    void recover_register(T& generator, int register_index, int register_lenght);

    register_recovery(float alpha_q, float beta_q);
};