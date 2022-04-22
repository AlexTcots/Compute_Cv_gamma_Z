//
// Created by alext on 4/21/2022.
//

#include "Particle.h"

Particle::Particle() : Z_trans(0), Z_int(0), Z_total(0),
                       e_trans(0), e_int(0), e_total(0),
                       Cv_trans(0), Cv_int(0), Cv_total(0) {

}

void Particle::Compute_Ztotal() {
    Z_total = Z_trans * Z_int;

}

void Particle::Compute_etotal() {
    e_total = e_trans + e_int;

}

void Particle::Compute_Cv_total() {
    Cv_total = Cv_trans + Cv_int;

}

void Particle::Set_Z_trans(double input) {
    Z_trans = input;

}

void Particle::Set_Z_int(double input) {
    Z_int = input;

}

void Particle::Set_etrans(double input) {
    e_trans = input;

}

void Particle::Set_eint(double input) {
    e_int = input;

}

void Particle::Set_Cv_trans(double input) {
    Cv_trans = input;

}

void Particle::Set_Cv_int(double input) {
    Cv_int = input;

}

const double &Particle::Get_Ztrans() const {
    return Z_trans;
}

const double &Particle::Get_Zint() const {
    return Z_int;
}

const double &Particle::Get_Ztotal() const {
    return Z_total;
}

const double &Particle::Get_etrans() const {
    return e_trans;
}

const double &Particle::Get_eint() const {
    return e_int;
}

const double &Particle::Get_etotal() const {
    return e_total;
}

const double &Particle::Get_Cv_trans() const {
    return Cv_trans;
}

const double &Particle::Get_Cv_int() const {
    return Cv_int;
}

const double &Particle::Get_Cv_total() const {
    return Cv_total;
}






