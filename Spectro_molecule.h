//
// Created by alext on 4/21/2022.
//

#ifndef COMPUTE_CV_GAMMA_Z_SPECTRO_MOLECULE_H
#define COMPUTE_CV_GAMMA_Z_SPECTRO_MOLECULE_H

#include <vector>
#include <fstream>
#include <cmath>
#include <iostream>
#include <sstream>
#include <exception>
#include "Constants.h"

struct molecu_data {
    double Te;
    double gn;
    double omega_e;
    double omega_e_xe;
    double omega_e_ye;
    double omega_e_ze;
    double omega_e_ke;
    double E_diss;
    double Be;
    double alpha_e;
    double gamma_e;
    double delta_e;
    double xi_e;
    double De;
    double beta_e;
    double g_e;
    double Hv0;
    double Hv1;
    double r_e;
};

class Spectro_molecule {
private:
    std::vector<molecu_data> data;
public:
    Spectro_molecule(int n_max, std::string file_name);
    Spectro_molecule(const Spectro_molecule& other);

    molecu_data Get_data_at(int n) const;

    const size_t Get_data_size()const;


};


#endif //COMPUTE_CV_GAMMA_Z_SPECTRO_MOLECULE_H
