//
// Created by alext on 4/21/2022.
//

#ifndef COMPUTE_CV_GAMMA_Z_SPECTRO_ATOM_H
#define COMPUTE_CV_GAMMA_Z_SPECTRO_ATOM_H
#include <vector>
#include <fstream>
#include <cmath>
#include <iostream>
#include <sstream>
#include <exception>
#include "Constants.h"
struct atom_data {
    double gn;
    double el_levl;
};
class Spectro_Atom {
private:
    std::vector<atom_data> data;
public:
   Spectro_Atom(int n_max,std::string file_name);
   Spectro_Atom(const Spectro_Atom& other);

    atom_data Get_data_at(int n) const;
    const size_t Get_data_size()const;



};


#endif //COMPUTE_CV_GAMMA_Z_SPECTRO_ATOM_H
