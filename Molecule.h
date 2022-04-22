//
// Created by alext on 4/22/2022.
//

#ifndef COMPUTE_CV_GAMMA_Z_MOLECULE_H
#define COMPUTE_CV_GAMMA_Z_MOLECULE_H

#include "Particle.h"
#include "Spectro_molecule.h"
#include "Constants.h"
#include <cmath>

class Molecule : protected Particle {
private:
    double mass;
    Spectro_molecule Molecule_data;
public:
    Molecule(double m, const Spectro_molecule &molecule_data);


    double Compute_en_ele(int n_elev);

    double Compute_en_vib(int n_elev, int i_vib);

    double Compute_en_rot(int n_elev, int i_vib, int j_rot);

    bool LessEdiss(int n_elev,int i_vib,int j_rot);

    void Compute_Ztrans(const double &Temp) override;


    void Compute_Zint(const double &Temp) override;

    void Compute_etrans(const double &Temp) override;

    void Compute_eint(const double &Temp) override;

    void Compute_Cv_trans(const double &Temp) override;

    void Compute_Cv_int(const double &Temp) override;

    double Get_Ztrans() const;

    double Get_Zint() const;

    double Get_Ztotal() const;

    double Get_etrans() const;

    double Get_eint() const;

    double Get_etotal() const;

    double Get_Cv_trans() const;

    double Get_Cv_int() const;

    double Get_Cv_total() const;

    void Compute_Z_e_cv(const double Temp);

};


#endif //COMPUTE_CV_GAMMA_Z_MOLECULE_H
