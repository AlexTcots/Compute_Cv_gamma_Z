//
// Created by alext on 4/21/2022.
//

#ifndef COMPUTE_CV_GAMMA_Z_ATOM_H
#define COMPUTE_CV_GAMMA_Z_ATOM_H

#include "Spectro_Atom.h"
#include "Particle.h"
#include "Constants.h"
#include <cmath>
class Atom : protected Particle {
private:
    double mass;
    Spectro_Atom Atom_data;
public:
    Atom(double m, const Spectro_Atom &atom_data);

    void Compute_Ztrans(const double& Temp) override;

    void Compute_Zint(const double& Temp) override;

    double Compute_en_ele(int n_elev);

    void Compute_etrans(const double& Temp) override;

    void Compute_eint(const double& Temp) override;

    void Compute_Cv_trans(const double& Temp) override;

    void Compute_Cv_int(const double& Temp) override;

    double Get_Ztrans()const;

    double Get_Zint()const;

    double Get_Ztotal()const;

    double Get_etrans()const;

    double Get_eint()const;

    double Get_etotal()const;

    double Get_Cv_trans()const;

    double Get_Cv_int()const;

    double Get_Cv_total()const;

    void Compute_Z_e_cv(const double Temp);

};


#endif //COMPUTE_CV_GAMMA_Z_ATOM_H
