#include "register_recovery.hpp"

#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

register_recovery::register_recovery(float alpha_q, float beta_q):alpha_quantile(alpha_q),beta_quantile(beta_q){}

// alpha - 2.336
// beta - 6.009, 6.121
void register_recovery::set_quantiles(float alpha_q, float beta_q){
    alpha_quantile = alpha_q;
    beta_quantile = beta_q;
    parameters_set = false;
}

void register_recovery::set_critical_set(){
    population = static_cast<int>(pow((beta_quantile*sqrt(0.25)+alpha_quantile*sqrt(0.1875))/(0.25), 2));
    criterion = static_cast<float>(population*(0.25)+alpha_quantile*sqrt(population*(0.25)));

    parameters_set = true;
}

void register_recovery::set_gamma_template(std::string& sequence){
    if(!parameters_set){
        throw runtime_error("Undefined population size");
    }
    if(sequence.size() != population){
        throw invalid_argument("Invalid template size");
    }
    if(!std::all_of(sequence.begin(), sequence.end(), [](char c) { return c == '0' || c == '1'; })) {
        throw std::invalid_argument("Invalid template");
    }
    gamma_template.reset();
    for(int i = 0; i < population; i++){
        if(sequence[i] == '1'){
            gamma_template.set(i);
        }
    }
}

void register_recovery::set_full_gamma_template(std::string& sequence){
    if(!std::all_of(sequence.begin(), sequence.end(), [](char c) { return c == '0' || c == '1'; })) {
        throw std::invalid_argument("Invalid template");
    }
    full_gamma_template.reset();
    for(int i = 0; i < population; i++){
        if(sequence[i] == '1'){
            full_gamma_template.set(i);
        }
    }
}

// void register_recovery::prepare_file(const string& out_file){
//     ofstream out(out_file, ios::binary);
//     if(!out){
//         throw runtime_error("Unable to open file");
//     }

//     uint64_t fuse = 0;
//     for(uint32_t i = 0; true; i++){
//         out.write(reinterpret_cast<char*>(&i), sizeof(i));
//         if(i == UINT32_MAX){
//             break;
//         }

//         fuse++;
//         if(fuse == 0x100000000){
//             throw runtime_error("Fuse lock");
//             break;
//         }
//     }
// }

bool register_recovery::recognize(std::bitset<320>& gamma){
    return ((gamma^gamma_template).count() < criterion);
}

void register_recovery::recover_L1(){
    uint32_t chunk_num = 0x40000000/16384;
    concurrency::parallel_for(uint32_t(0), chunk_num, [&](uint32_t chunk_index){
        bitset<320> gamma;
        lfsr l_register(30, 0x32800000);
        for(uint32_t i = chunk_index*16384; i < (chunk_index+1)*16384; i++){
            gamma.reset();           
            l_register.set_register(i);
            for(int j = 0; j < population; j++){
                gamma[j] = l_register.fast_clock();
            }
            if(recognize(gamma)){
                L1_candidates.push_back(i);
            }
        }
    }
    );
    cout << "L1 candidates are found.\n";
    // for(int i = 0; i < static_cast<int>(L1_candidates.size()); i++){
    //     cout << L1_candidates[i] << "\n";
    // }
}

void register_recovery::recover_L2(){
    uint32_t chunk_num = 0x80000000/32768;
    concurrency::parallel_for(uint32_t(0), chunk_num, [&](uint32_t chunk_index){
        bitset<320> gamma;
        lfsr l_register(31, 0x48000000);
        for(uint32_t i = chunk_index*32768; i < (chunk_index+1)*32768; i++){
            gamma.reset();           
            l_register.set_register(i);
            for(int j = 0; j < population; j++){
                gamma[j] = l_register.fast_clock();
            }
            if(recognize(gamma)){
                L2_candidates.push_back(i);
            }
        }
    }
    );
    cout << "L2 candidates are found.\n";
    // for(int i = 0; i < static_cast<int>(L2_candidates.size()); i++){
    //     cout << L2_candidates[i] << "\n";
    // }
    // cout << L2_candidates.size();
}

void register_recovery::recover_L3(){
    const uint32_t gamma_template_32 = 0x82AB0478;    
    for(int i = 0; i < static_cast<int>(L2_candidates.size()); i++){
        uint32_t mask_1, mask_0;
        for(int j = 0; j < static_cast<int>(L1_candidates.size()); j++){
            if(((~(L1_candidates[j]^L2_candidates[i]))&(L2_candidates[i]^gamma_template_32))!=0){
                continue;
            }

            unsigned int m_1_count, m_0_count;
            //1 if 1 in the mask is set, 0 in the other case
            mask_1 = (L1_candidates[j]^L2_candidates[i])&(L2_candidates[i]^gamma_template_32);
            m_1_count = __popcnt(mask_1);
            //0 if 0 in the mask is set, 1 in the other case
            mask_0 = (~(L1_candidates[j]^L2_candidates[i]))|(L2_candidates[i]^gamma_template_32);
            m_0_count = __popcnt(mask_0);
            
            uint64_t chunk_num = 0x100000000ull/131072;
            atomic_bool term = false;
            concurrency::parallel_for(uint64_t(0), chunk_num, [&](uint64_t chunk_index){
                if(!term.load()){
                    geffe_generator gn;
                    uint32_t sample;
                    bool bit;
                    for(uint64_t k = chunk_index*131072; k < (chunk_index+1)*131072; k++){
                        gn.set_register(0, L1_candidates[j]);
                        gn.set_register(1, L2_candidates[i]);                       
                        gn.set_register(2, static_cast<uint32_t>(k));
                        sample = 0;
                        for(uint8_t s = 0; s < 32; s++){
                            bit = gn.clock();
                            if(bit){
                                sample |= 1<<s;
                            }
                        }
                        if(__popcnt(sample&mask_1)!=m_1_count){
                            continue;
                        }
                        if(__popcnt(sample|mask_0)!=m_0_count){
                            continue;
                        }
                        gn.set_register(0, L1_candidates[j]);
                        gn.set_register(1, L2_candidates[i]);
                        gn.set_register(2, static_cast<uint32_t>(k));
                        bool flag = true;
                        for(int s = 0; s < static_cast<int>(full_gamma_template.size()); s++){
                            if(!(gn.clock()&full_gamma_template[s])){
                                flag = false;
                                break;
                            }
                        }
                        if(flag){
                            cout << L1_candidates[j] << " " << L2_candidates[i] << " " << static_cast<uint32_t>(k) << "\n";
                            term.store(true);
                            return;
                        }
                    }                    
                }
            }
            );
        }
    }
}