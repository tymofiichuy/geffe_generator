#include "register_recovery.hpp"

using namespace std;

int main(){
    lfsr l1_register(31, 0x32800000);
    l1_register.set_register(537051710);
    lfsr l2_register(31, 0x48000000);
    l2_register.set_register(1477968685);
    for(int i = 0; i < 265; i++){
        cout << l1_register.fast_clock();
    }
    cout << "\n";
    for(int i = 0; i < 265; i++){
        cout << l2_register.fast_clock();
    }
    return 0;
}