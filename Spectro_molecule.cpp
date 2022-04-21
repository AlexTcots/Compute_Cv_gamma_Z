//
// Created by alext on 4/21/2022.
//

#include "Spectro_molecule.h"

Spectro_molecule::Spectro_molecule(int n_max, std::string file_name) {
    data.reserve(n_max);
    std::ifstream ifs(file_name);
    if (!ifs.is_open()) {
        throw std::invalid_argument("fail to open file " + file_name);
    } else {
        std::string line;
        for (int i = 0; i < n_max; ++i) {
            getline(ifs, line);
            std::istringstream iss(line);
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
            iss >> Te
                >> gn
                >> omega_e
                >> omega_e_xe
                >> omega_e_ye
                >> omega_e_ze
                >> omega_e_ke
                >> E_diss
                >> Be
                >> alpha_e
                >> gamma_e
                >> delta_e
                >> xi_e
                >> De
                >> beta_e
                >> g_e
                >> Hv0
                >> Hv1
                >> r_e;
            data.push_back(
                    {Te, gn, omega_e, omega_e_xe, omega_e_ye, omega_e_ze, omega_e_ke, E_diss, Be, alpha_e, gamma_e,
                     delta_e, xi_e, De, beta_e, g_e, Hv0, Hv1, r_e});
        }
    }
}

molecu_data Spectro_molecule::Get_data_at(int n) const {
    return data.at(n);
}

const size_t Spectro_molecule::Get_data_size() const {
    return data.size();
}

Spectro_molecule::Spectro_molecule(const Spectro_molecule &other) {
    data=other.data;

}


