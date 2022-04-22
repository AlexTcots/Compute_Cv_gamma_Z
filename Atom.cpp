//
// Created by alext on 4/21/2022.
//

#include "Atom.h"

Atom::Atom(double m, const Spectro_Atom &atom_data) : Particle(),
                                                      mass(m), Atom_data(atom_data) {
    std::cout << " Atom constructed !" << '\n';

}

void Atom::Compute_Ztrans(const double& Temp) {
    double ztrans = pow(2 * Pi * mass * BOLTZ * Temp / (h_plc * h_plc), 3.0 / 2);
    Set_Z_trans(ztrans);

}

void Atom::Compute_Zint(const double& Temp) {
    size_t nmax = Atom_data.Get_data_size();
    double zint=0;
    for(size_t i=0;i<nmax;++i){

        zint+=Atom_data.Get_data_at(i).gn* exp(-1*Compute_en_ele(i)/(BOLTZ*Temp));
    }
    Set_Z_int(zint);
    // Z total
    Compute_Ztotal();
}

double Atom::Compute_en_ele(int n_elev) {
    return Atom_data.Get_data_at(n_elev).el_levl*100*h_plc*c_light;
}

void Atom::Compute_etrans(const double& Temp) {
    double etrans = 1.5*BOLTZ*Temp/mass;
    Set_etrans(etrans);

}

void Atom::Compute_eint(const double& Temp) {
    double eint =0;
    // calculate the sum part
    for(size_t n=0;n<Atom_data.Get_data_size();++n){
        double e_ele= Compute_en_ele(n);
        eint+= Atom_data.Get_data_at(n).gn*e_ele* exp(e_ele/(-BOLTZ*Temp));
    }
    Set_eint(eint/(mass*Get_Zint()));
    // etotal
    Compute_etotal();

}

void Atom::Compute_Cv_trans(const double& Temp) {
    Set_Cv_trans(1.5*BOLTZ/mass);


}

void Atom::Compute_Cv_int(const double& Temp) {
    double lhs,rhs;
    lhs=rhs=0;
    size_t nmax = Atom_data.Get_data_size();
    for(size_t n=0;n<nmax;++n){
        double en_ele = Compute_en_ele(n);

        lhs+= pow(en_ele/(BOLTZ*Temp),2.0)*Atom_data.Get_data_at(n).gn*exp(-en_ele/(BOLTZ*Temp));
        rhs+= (en_ele/(BOLTZ*Temp))*Atom_data.Get_data_at(n).gn* exp(-en_ele/(BOLTZ*Temp));

    }

    double cv_int = (BOLTZ/mass)*(lhs/Get_Zint()-pow(rhs/Get_Zint(),2));
    Set_Cv_int(cv_int);
    // Cv_total
    Compute_Cv_total();

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

double Atom::Get_etrans() const {
    return Particle::Get_etrans();
}

double Atom::Get_eint() const {
    return Particle::Get_eint();
}

double Atom::Get_etotal() const {
    return Particle::Get_etotal();
}

double Atom::Get_Cv_trans() const {
    return Particle::Get_Cv_trans();
}

double Atom::Get_Cv_int() const {
    return Particle::Get_Cv_int();
}

double Atom::Get_Cv_total() const {
    return Particle::Get_Cv_total();
}

void Atom::Compute_Z_e_cv(const double Temp) {
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




