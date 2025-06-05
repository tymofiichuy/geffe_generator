#include "register_recovery.hpp"

using namespace std;

int main(){
    lfsr l_register(31, 0x48000000);
    l_register.set_register(1741154825);
    for(int i = 0; i < 265; i++){
        cout << l_register.fast_clock();
    }
    return 0;
}