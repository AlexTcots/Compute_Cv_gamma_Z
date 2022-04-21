//
// Created by alext on 4/21/2022.
//

#include "Atom.h"

Atom::Atom(double m, const Spectro_Atom &atom_data) : Particle(),
                                                      mass(m), Atom_data(atom_data) {
    std::cout << " Atom constructed !" << '\n';

}

void Atom::Compute_Ztrans(double Temp) {
    double ztrans = pow(2 * Pi * mass * BOLTZ * Temp / (h_plc * h_plc), 3 / 2);
    Set_Z_trans(ztrans);

}

void Atom::Compute_Zint(double Temp) {
    size_t nmax = Atom_data.Get_data_size();
    double zint=0;
    for(size_t i=0;i<nmax;++i){

        zint+=Atom_data.Get_data_at(i).gn* exp(-1*Compute_en_ele(i)/(BOLTZ*Temp));
    }
    Set_Z_int(zint);
    Compute_Ztotal();
}

double Atom::Compute_en_ele(int n_elev) {
    return Atom_data.Get_data_at(n_elev).el_levl*100*h_plc*c_light;
}

void Atom::Compute_etrans(double Temp) {
    double etrans = 1.5*BOLTZ*Temp/mass;
    Set_etrans(etrans);

}

void Atom::Compute_eint(double Temp) {
    double eint =0;
    // calculate the sum part
    for(size_t n=0;n<Atom_data.Get_data_size();++n){
        double e_ele= Compute_en_ele(n);
        eint+= Atom_data.Get_data_at(n).gn*e_ele* exp(e_ele/(-BOLTZ*Temp));
    }
    Set_eint(eint/(mass*Get_Zint()));
    Compute_etotal();

}

void Atom::Compute_Cv_trans(double Temp) {


}

void Atom::Compute_Cv_int(double Temp) {

}

double Atom::Get_Ztrans() const {
    return Particle::Get_Ztrans();
}

double Atom::Get_Zint() const {
    return Particle::Get_Zint();
}

double Atom::Get_Ztotal() const {
    return Particle::Get_Ztotal();
}


