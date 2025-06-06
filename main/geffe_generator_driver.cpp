#include "geffe_generator.hpp"
#include <string>

using namespace std;

//806014269 55649069 2352825186
int main(int argc, char *argv[]){
    if(argc != 2){
        return 1;
    }
    geffe_generator gf;
    gf.set_register(806014269, 0);
    gf.set_register(55649069, 1);
    gf.set_register(2352825186, 2);

    for(int i = 0; i < stoi(argv[1]); i++){
        cout << static_cast<int>(gf.clock());
    }

    return 0;
}