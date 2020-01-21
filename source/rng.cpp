//
// Created by ilaria on 2020-01-09.
//

#include "rng.h"

/*To do list:
 * insert all the distribution I need here DONE
 * add the seeding procedure to the main DONE
 * remove rng from MC_parameters struct DONE
 * insert rn::uniform_integer_box() ecc in all the functions extracting random numbers (montecarlo, parallel tempering..) DONE
 * (not urgent) print the seed used in a file in the correct folder
 * */
namespace rn{
    pcg32 rng;
    void seed(long n){
        if(n>=0) {
            auto given_seed = (unsigned long) n;
            std::cout << "Seeding : " << given_seed << std::endl;
            std::seed_seq seq{given_seed};
            rng.seed(seq);
            std::srand(rng());
        }
        else{
            std::cout << "Seeding : std::random_device" << std::endl;
            pcg_extras::seed_seq_from<std::random_device> seed_source;
            rng.seed(seed_source);
            std::srand(rng());
        }
    }


    int uniform_integer_box(const int min, const int max){
        std::uniform_int_distribution<>  rand_int(std::min(min,max),std::max(min,max));
        return rand_int(rng);
    }

    double uniform_real_box(const double min, const double max){
        std::uniform_real_distribution<>  rand_real(std::min(min,max),std::max(min,max));
        return rand_real(rng);
    }


}
