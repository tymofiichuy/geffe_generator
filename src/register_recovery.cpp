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
    for(int i = 0; i < static_cast<int>(sequence.size()); i++){
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
    cout << "L1 candidates are found: " << L1_candidates.size() << " items.\n";
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
    cout << "L2 candidates are found: " << L2_candidates.size() << " items.\n";
}

void register_recovery::recover_L3(){
    int counter_320 = 0, counter_32 = 0;
    const uint32_t gamma_template_32 = 0x82AB0478;
    bitset<320> gamma_template_320("10000010101010110000010001111000001000110010000011101101000100101010110111000110010000010001011001110100001000101011011001000011000011010101011000000111110011010001111100110110100010011110010111101111110101000100110010000101011011010111100101100011000011101101010000100100101100110111011111011100101101000001110011110111");
    uint32_t mask_1, mask_0, not_mask_0;
    uint32_t L1_sample_32, L2_sample_32;
    bitset<320> L1_sample_320, L2_sample_320;
    lfsr L1_lfsr(30, 0x32800000), L2_lfsr(31, 0x48000000);    
    for(int i = 0; i < static_cast<int>(L2_candidates.size()); i++){
        for(int j = 0; j < static_cast<int>(L1_candidates.size()); j++){
            L1_lfsr.set_register(L1_candidates[j]);
            L2_lfsr.set_register(L2_candidates[i]);
            L1_sample_32 = 0;
            L2_sample_32 = 0;
            L1_sample_320.reset();
            L2_sample_320.reset();
            for(int l = 319; l >= 0; l--){
                if(L1_lfsr.fast_clock()){
                    L1_sample_320.set(l);
                }
                if(L2_lfsr.fast_clock()){
                    L2_sample_320.set(l);
                } 
            }
            for(int l = 319, bit = 31; l > 287; l--, bit--){
                if(L1_sample_320[l]==1){
                    L1_sample_32 |= (1<<bit);
                }
                if(L2_sample_320[l]==1){
                    L2_sample_32 |= (1<<bit);
                }
            }
            if(((~(L1_sample_320^L2_sample_320))&(L2_sample_320^gamma_template_320)).any()){
                continue;
            }  
            unsigned int m_1_count, m_0_count;
            //1 if 1 in the mask is set, 0 in the other case
            mask_1 = (L1_sample_32^L2_sample_32)&(L2_sample_32^gamma_template_32);
            m_1_count = __popcnt(mask_1);
            //0 if 0 in the mask is set, 1 in the other case
            mask_0 = (~(L1_sample_32^L2_sample_32))|(L2_sample_32^gamma_template_32);
            m_0_count = __popcnt(mask_0);
            not_mask_0 = ~mask_0; 
            
            uint64_t sh = (static_cast<uint64_t>(1)<<__popcnt(~(not_mask_0|mask_1)));
            uint64_t chunk_size = 131072;
            uint64_t chunk_num = sh/chunk_size;
            while(chunk_num < 16){
                chunk_size >>= 1;                
                chunk_num = sh/chunk_size;
            }

            atomic_bool term;
            term.store(false);  
            concurrency::parallel_for(uint64_t(0), chunk_num, [&](uint64_t chunk_index){
                if(term.load()==false){
                    geffe_generator gn;
                    uint32_t L3_candidate;
                    for(uint64_t k = chunk_index*chunk_size; k < (chunk_index+1)*chunk_size; k++){
                        uint32_t iteration = static_cast<uint32_t>(k);
                        int counter = 0;
                        L3_candidate = 0;
                        for(int t = 31; t >= 0; t--){
                            if((mask_1&(1<<t))!=0){
                                L3_candidate ^= (1<<t);
                            }
                            else if((not_mask_0&(1<<t))!=0){
                                continue;
                            }
                            else{
                                if((iteration&(1<<counter))!=0){
                                    L3_candidate ^= (1<<t);
                                }
                                counter++;
                            }
                        }
                        gn.set_register(L1_candidates[j], 0);
                        gn.set_register(L2_candidates[i], 1);                       
                        gn.set_register(L3_candidate, 2);
                        bool flag = true;
                        for(int s = 0; s < static_cast<int>(full_gamma_template.size()); s++){
                            if(gn.clock()!=full_gamma_template[s]){
                                flag = false;
                                break;
                            }
                        }  
                        if(flag){
                            cout << full_gamma_template;
                            cout << "L1: " << L1_candidates[j] << ", L2: " << L2_candidates[i] << ", L3: " << L3_candidate << "\n";
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