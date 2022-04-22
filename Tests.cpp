//
// Created by alext on 4/21/2022.
//

#include "Tests.h"

void Test_spectro_data() {
    Spectro_Atom N_aotm(20, "N_atom.txt");

    for (int n = 0; n < 20; n++) {
        cout << "data of N atom:" << '\n';
        cout << N_aotm.Get_data_at(n).gn << ' ' << N_aotm.Get_data_at(n).el_levl << '\n';
    }
    Spectro_molecule O2(5, "O2.txt");
    for (int n = 0; n < 5; n++) {
        cout << "data of O2:" << '\n';
        cout << O2.Get_data_at(n).gn << ' ' << O2.Get_data_at(n).Te << '\n';
    }

}

void Test_Atom() {
    double temp = 300;
    Spectro_Atom N_atom(20, "N_atom.txt");
    Atom N(m_N, N_atom);
    N.Compute_Z_e_cv(10000);
    cout << "Z_trans = " << N.Get_Ztrans() << ' '
         << "Z_int   = " << N.Get_Zint() << ' '
         << "Z_total = " << N.Get_Ztotal() << '\n'
         << "e_trans = " << N.Get_etrans() << ' '
         << "e_int   = " << N.Get_eint() << ' '
         << "e_total = " << N.Get_etotal() << '\n'
         << "Cv_trans = " << N.Get_Cv_trans() << ' '
         << "Cv_int   = " << N.Get_Cv_int() << ' '
         << "Cv_total = " << N.Get_Cv_total() << '\n';


}
