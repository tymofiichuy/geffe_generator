#include "geffe_generator.hpp"
#include <string>

using namespace std;

int main(int argc, char *argv[]){
    if(argc != 2){
        return 1;
    }
    geffe_generator gf;
    gf.set_register(0x2AAAAAAA, 0);
    gf.set_register(0x2AAAAAAA, 1);
    gf.set_register(0xCCCCCCCC, 2);

    for(int i = 0; i < stoi(argv[1]); i++){
        cout << static_cast<int>(gf.clock());
    }

    return 0;
}