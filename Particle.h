//
// Created by alext on 4/21/2022.
//

#ifndef COMPUTE_CV_GAMMA_Z_PARTICLE_H
#define COMPUTE_CV_GAMMA_Z_PARTICLE_H


class Particle {
private:
    double Z_trans, Z_int, Z_total;
    double e_trans, e_int, e_total;
    double Cv_trans, Cv_int, Cv_total;
protected:
    void Set_Z_trans(double input);

    void Set_Z_int(double input);

    void Set_etrans(double input);

    void Set_eint(double input);

public:
    Particle();

    virtual void Compute_Ztrans(double Temp) = 0;

    virtual void Compute_Zint(double Temp) = 0;

    void Compute_Ztotal();

    virtual void Compute_etrans(double Temp) = 0;

    virtual void Compute_eint(double Temp) = 0;

    void Compute_etotal();

    virtual void Compute_Cv_trans(double Temp) = 0;

    virtual void Compute_Cv_int(double Temp) = 0;

    void Compute_Cv_total();

    const double &Get_Ztrans() const;

    const double &Get_Zint() const;

    const double &Get_Ztotal() const;

    const double &Get_etrans() const;

    const double &Get_eint() const;

    const double &Get_etotal() const;

    const double &Get_Cv_trans() const;

    const double &Get_Cv_int() const;

    const double &Get_Cv_total() const;

};


#endif //COMPUTE_CV_GAMMA_Z_PARTICLE_H
