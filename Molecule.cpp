//
// Created by alext on 4/22/2022.
//

#include "Molecule.h"

Molecule::Molecule(double m, const Spectro_molecule &molecule_data)
        : Particle(), mass(m), Molecule_data(molecule_data) {
    std::cout << "Molecule constructed!" << '\n';
}

double Molecule::Compute_en_ele(int n_elev) {
    return Molecule_data.Get_data_at(n_elev).Te * 100 * h_plc * c_light;
}

double Molecule::Compute_en_vib(int n_elev, int i_vib) {
    double omegae_n = Molecule_data.Get_data_at(n_elev).omega_e;
    double omegae_xe_n = Molecule_data.Get_data_at(n_elev).omega_e_xe;
    return omegae_n * (i_vib + 0.5) - omegae_xe_n * pow(i_vib + 0.5, 2.0);
}

double Molecule::Compute_en_rot(int n_elev, int i_vib, int j_rot) {
    double B_ne = Molecule_data.Get_data_at(n_elev).Be;
    double alpha_ne = Molecule_data.Get_data_at(n_elev).alpha_e;
    double D_ne = Molecule_data.Get_data_at(n_elev).De;
    double beta_ne = Molecule_data.Get_data_at(n_elev).beta_e;

    double B_ni = B_ne - alpha_ne * (i_vib + 0.5);
    double D_ni = D_ne - beta_ne * (i_vib + 0.5);

    return B_ni * j_rot * (j_rot + 1) - D_ni * pow(j_rot, 2.0) * pow(j_rot + 1, 2.0);
}

bool Molecule::LessEdiss(int n_elev, int i_vib, int j_rot) {
    double E_diss = Molecule_data.Get_data_at(n_elev).E_diss * 100 * h_plc * c_light;
    double E_nij = Compute_en_ele(n_elev) + Compute_en_vib(n_elev, i_vib) +
                   Compute_en_rot(n_elev, i_vib, j_rot);
    return E_nij < E_diss;
}

void Molecule::Compute_Ztrans(const double &Temp) {
    double ztrans = pow(2*Pi*mass*BOLTZ* Temp/(h_plc*h_plc),3.0/2);
    Set_Z_trans(ztrans);

}


