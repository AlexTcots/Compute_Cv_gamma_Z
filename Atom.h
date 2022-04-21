//
// Created by alext on 4/21/2022.
//

#ifndef COMPUTE_CV_GAMMA_Z_ATOM_H
#define COMPUTE_CV_GAMMA_Z_ATOM_H

#include "Spectro_Atom.h"
#include "Particle.h"
#include "Constants.h"

class Atom : protected Particle {
private:
    double mass;
    Spectro_Atom Atom_data;
public:
    Atom(double m, const Spectro_Atom &atom_data);

    void Compute_Ztrans(double Temp) override;

    void Compute_Zint(double Temp) override;

    double Compute_en_ele(int n_elev);

    void Compute_etrans(double Temp) override;

    void Compute_eint(double Temp) override;

    void Compute_Cv_trans(double Temp);

    void Compute_Cv_int(double Temp);

    double Get_Ztrans()const;

    double Get_Zint()const;

    double Get_Ztotal()const;

};


#endif //COMPUTE_CV_GAMMA_Z_ATOM_H
