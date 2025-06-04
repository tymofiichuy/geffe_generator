#include "register_recovery.hpp"
#include <fstream>
#include <cmath>

using namespace std;

template<typename T> register_recovery<T>::register_recovery(float alpha_q, float beta_q):alpha_quantile(alpha_q),beta_quantile(beta_q){}

template<typename T> void register_recovery<T>::set_quantiles(float alpha_q, float beta_q){
    alpha_quantile = alpha_q;
    beta_quantile = beta_q;
}

template<typename T> void register_recovery<T>::prepare_file(const string& out_file){
    ofstream out(out_file, ios::binary);
    if(!out){
        throw runtime_error("Unable to open file");
    }

    uint64_t fuse = 0;
    for(uint32_t i = 0; true; i++){
        out.write(reinterpret_cast<char*>(&i), sizeof(i));
        if(i == UINT32_MAX){
            break;
        }

        fuse++;
        if(fuse == 0x100000000){
            throw runtime_error("Fuse lock");
            break;
        }
    }
}

template<typename T> bool register_recovery<T>::recognize(){
    
}