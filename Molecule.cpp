//
// Created by alext on 4/22/2022.
//

#include "Molecule.h"

Molecule::Molecule(double m, const Spectro_molecule &molecule_data)
        : Particle(), mass(m), Molecule_data(molecule_data) {
    std::cout << "Molecule constructed!" << '\n';
}

double Molecule::Compute_en_ele(int n_elev) {
    return Molecule_data.Get_data_at(n_elev).Te * 1e2 * h_plc * c_light;
}

double Molecule::Compute_en_vib(int n_elev, int i_vib) {
    double omegae_n = Molecule_data.Get_data_at(n_elev).omega_e;
    double omegae_xe_n = Molecule_data.Get_data_at(n_elev).omega_e_xe;
    return (1e2 * omegae_n * (i_vib + 0.5) - 1e2 * omegae_xe_n * pow(i_vib + 0.5, 2.0)) * h_plc * c_light;
}

double Molecule::Compute_en_rot(int n_elev, int i_vib, int j_rot) {
    double B_ne =  100*Molecule_data.Get_data_at(n_elev).Be;
    double alpha_ne = Molecule_data.Get_data_at(n_elev).alpha_e;
    //double D_ne = Molecule_data.Get_data_at(n_elev).De;
    //double beta_ne = Molecule_data.Get_data_at(n_elev).beta_e;

    double B_ni = B_ne - alpha_ne * (i_vib + 0.5);
    //double D_ni = D_ne - beta_ne * (i_vib + 0.5);

    //return (B_ni * j_rot * (j_rot + 1) - D_ni * pow(j_rot, 2.0) * pow(j_rot + 1, 2.0))*h_plc*c_light;
    return B_ni * j_rot * (j_rot + 1)* h_plc * c_light;
}

double Molecule::Compute_e_nij(int n_elev, int i_vib, int j_rot) {
    return Compute_en_ele(n_elev) + Compute_en_vib(n_elev, i_vib) +
           Compute_en_rot(n_elev, i_vib, j_rot);
}

double Molecule::Compute_e_diss(int n_elev) {
    double e_diss = Molecule_data.Get_data_at(n_elev).Te + fabs(Molecule_data.Get_data_at(n_elev).E_diss);

    return e_diss * 100 * h_plc * c_light;
}

bool Molecule::LessEdiss(int n_elev, int i_vib, int j_rot) {
    double E_diss = Compute_e_diss(n_elev);
    double E_nij = Compute_e_nij(n_elev, i_vib, j_rot);
    return E_nij < E_diss;
}

void Molecule::Compute_Ztrans(const double &Temp) {
    double ztrans = pow(2 * Pi * mass * BOLTZ * Temp / (h_plc * h_plc), 3.0 / 2);
    Set_Z_trans(ztrans);

}

void Molecule::Compute_Zint(const double &Temp) {
    size_t nmax = Molecule_data.Get_data_size();
    double sum = 0;
    double sum_j = 0;
    double sum_i = 0;
    double sum_n = 0;
    double e_nij = 0;
    int i_vib = 0;
    int j_rot = 0;
    double beta = -1*BOLTZ*Temp;

    for (int n_elevl = 0; n_elevl < nmax; ++n_elevl) {
        double e_dissn = Compute_e_diss(n_elevl);
        e_nij = Compute_e_nij(n_elevl, i_vib, j_rot);
        while (e_nij < e_dissn) {// loop for vib
            while (e_nij < e_dissn) { // loop for rot

                sum_j += (2 * j_rot + 1) * exp(Compute_en_rot(n_elevl, i_vib, j_rot) / beta);
                j_rot++;
                //std::cout << " n i j " << n_elevl << ' ' << i_vib << ' ' << j_rot << '\n';
                e_nij = Compute_e_nij(n_elevl, i_vib, j_rot);
            }
            sum_i += exp(Compute_en_vib(n_elevl, i_vib) / beta) * sum_j;
            sum_j = 0;
            i_vib++;
            j_rot = 0;
            e_nij = Compute_e_nij(n_elevl, i_vib, j_rot);
        }
        sum_n += Molecule_data.Get_data_at(n_elevl).gn * exp(Compute_en_ele(n_elevl) / beta) * sum_i;
        sum_i = 0;
        i_vib = 0;

    }

    /*
    for (int n_elevl = 0; n_elevl < nmax; ++n_elevl) {
        double e_dissn = Compute_e_diss(n_elevl);
        while (e_nij < e_dissn) {
            while (e_nij < e_dissn) {
                e_nij = Compute_e_nij(n_elevl, i_vib, j_rot);
                sum += Molecule_data.Get_data_at(n_elevl).gn * (2 * j_rot + 1) * exp(e_nij / (-1 * BOLTZ * Temp));
                j_rot++;
                //std::cout << " n i j " << n_elevl << ' ' << i_vib << ' ' << j_rot << '\n';

            }
            i_vib++;
            j_rot = 0;
            e_nij = Compute_e_nij(n_elevl, i_vib, j_rot);
        }
        i_vib = 0;
        e_nij = Compute_e_nij(n_elevl,i_vib,j_rot);

    }
     */

    Set_Z_int(sum_n/2);
    // Z_tatol
    Compute_Ztotal();

}

void Molecule::Compute_etrans(const double &Temp) {
    double etrans = 1.5 * BOLTZ * Temp / mass;
    Set_etrans(etrans);
}

void Molecule::Compute_eint(const double &Temp) {
    size_t nmax = Molecule_data.Get_data_size();
    double sum_j = 0;
    double sum_i = 0;
    double sum_n = 0;
    double e_nij = 0;
    double e_dissn = 0;
    int i_vib = 0;
    int j_rot = 0;
    double beta = -1 * BOLTZ * Temp;

    for (int n_elevl = 0; n_elevl < nmax; ++n_elevl) {
        e_dissn = Compute_e_diss(n_elevl);
        e_nij = Compute_e_nij(n_elevl, i_vib, j_rot);
        while (e_nij < e_dissn) {// loop for vib
            while (e_nij < e_dissn) { // loop for rot
                double erot = Compute_en_rot(n_elevl, i_vib, j_rot);
                sum_j += (2 * j_rot + 1) * erot * exp(erot / beta);
                j_rot++;
               //std::cout << " n i j " << n_elevl << ' ' << i_vib << ' ' << j_rot <<' '
               //<<"sum_j = "<< sum_j<<' '<<" sum_i = "<< sum_i<<" sum_n = "<< sum_n<< '\n';
                e_nij = Compute_e_nij(n_elevl, i_vib, j_rot);
            }
            double evib = Compute_en_vib(n_elevl, i_vib);
            sum_i += evib * exp(evib / beta) * sum_j;

            sum_j = 0;
            i_vib++;
            j_rot = 0;
            e_nij = Compute_e_nij(n_elevl, i_vib, j_rot);
        }
        double eele = Compute_en_ele(n_elevl);
        sum_n += Molecule_data.Get_data_at(n_elevl).gn * eele * exp(eele / beta) * sum_i;

        sum_i = 0;
        i_vib = 0;
        j_rot = 0;

    }
    Set_eint(sum_n / (mass * Get_Zint()));
    // etotal
    Compute_etotal();

}

void Molecule::Compute_Cv_trans(const double &Temp) {
    double cv_trans = 1.5 * BOLTZ / mass;
    Set_Cv_trans(cv_trans);

}

void Molecule::Compute_Cv_int(const double &Temp) {

    size_t nmax = Molecule_data.Get_data_size();
    double sum_jl = 0;
    double sum_il = 0;
    double sum_nl = 0;
    double sum_jr = 0;
    double sum_ir = 0;
    double sum_nr = 0;
    double e_nij = 0;
    int i_vib = 0;
    int j_rot = 0;
    double beta = -1*BOLTZ*Temp;

    for (int n_elevl = 0; n_elevl < nmax; ++n_elevl) {
        double e_dissn = Compute_e_diss(n_elevl);
        e_nij = Compute_e_nij(n_elevl, i_vib, j_rot);
        while (e_nij < e_dissn) {// loop for vib
            while (e_nij < e_dissn) { // loop for rot
                double erot = Compute_en_rot(n_elevl, i_vib, j_rot);
                sum_jl += pow(erot / (BOLTZ * Temp), 2.0) * (2 * j_rot + 1) * exp(erot / beta);
                sum_jr += (erot / (BOLTZ * Temp)) * (2 * j_rot + 1) * exp(erot / beta);
                j_rot++;
                //std::cout << " n i j " << n_elevl << ' ' << i_vib << ' ' << j_rot << '\n';
                e_nij = Compute_e_nij(n_elevl, i_vib, j_rot);
            }
            double evib = Compute_en_vib(n_elevl, i_vib);
            sum_il += pow(evib / (BOLTZ * Temp), 2.0) * exp(evib / beta) * sum_jl;
            sum_ir += (evib / (BOLTZ * Temp)) * exp(evib / beta) * sum_jr;
            sum_jl = 0;
            sum_jr = 0;
            i_vib++;
            j_rot = 0;
            e_nij = Compute_e_nij(n_elevl, i_vib, j_rot);
        }
        double eele = Compute_en_ele(n_elevl);
        sum_nl += pow(eele / (BOLTZ * Temp), 2.0) * Molecule_data.Get_data_at(n_elevl).gn * eele *
                  exp(eele / beta) * sum_il;
        sum_nr += (eele / (BOLTZ * Temp)) * Molecule_data.Get_data_at(n_elevl).gn * eele *
                  exp(eele / beta) * sum_ir;
        sum_il = 0;
        sum_ir = 0;
        i_vib = 0;

    }

    double cv_int = (BOLTZ / mass) * (sum_nl / Get_Zint() - pow(sum_nr / Get_Zint(), 2.0));
    Set_Cv_int(cv_int);
    // Cv_total
    Compute_Cv_total();

}

double Molecule::Get_Ztrans() const {
    return Particle::Get_Ztrans();
}

double Molecule::Get_Zint() const {
    return Particle::Get_Zint();
}

double Molecule::Get_Ztotal() const {
    return Particle::Get_Ztotal();
}

double Molecule::Get_etrans() const {
    return Particle::Get_etrans();
}

double Molecule::Get_eint() const {
    return Particle::Get_eint();
}

double Molecule::Get_etotal() const {
    return Particle::Get_etotal();
}

double Molecule::Get_Cv_trans() const {
    return Particle::Get_Cv_trans();
}

double Molecule::Get_Cv_int() const {
    return Particle::Get_Cv_int();
}

double Molecule::Get_Cv_total() const {
    return Particle::Get_Cv_total();
}

void Molecule::Compute_Z_e_cv(const double Temp) {
    // apply computing in order;

    // Z_trans
    Compute_Ztrans(Temp);
    // Z_int
    Compute_Zint(Temp);// Z_total update;
    // etrans
    Compute_etrans(Temp);
    // eint
    Compute_eint(Temp);// etotal update;
    // Cv_trans
    Compute_Cv_trans(Temp);
    // Cv_int
    Compute_Cv_int(Temp);// Cv_total update;
}








