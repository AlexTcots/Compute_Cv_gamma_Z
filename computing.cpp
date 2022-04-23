//
// Created by alext on 4/23/2022.
//

#include "computing.h"

void Computing_data(string input_filename,int size_vec,const double& temp1,const double& temp2,string output_name) {

    Spectro_molecule O2_data(5, input_filename);
    Molecule O2(m_O2, O2_data);

    ofstream output(output_name, ios::trunc | ios::out);
    output << "Z_trans  "  << ' '
           << "Z_int    "  << ' '
           << "Z_total  "  << ' '

           << "e_trans  "  << ' '
           << "e_int    "  << ' '
           << "e_total  "  << ' '

           << "Cv_trans  "  << ' '
           << "Cv_int    "  << ' '
           << "Cv_total  "  << '\n';
    double deltT = temp2 ;
    vector<double> Temps(size_vec);
    for (int i = 0; i < Temps.size(); ++i) {
        Temps[i] = temp1 + i * deltT;
    }
    int count =1;
    for (const auto &t: Temps) {
        O2.Compute_Z_e_cv(t);
        std::cout<< count <<' '<< "computing finished "<<'\n';
        ++count;

        output << t << ' '
               << O2.Get_Ztrans() << ' '
               << O2.Get_Zint() << ' '
               << O2.Get_Ztotal() << ' '

               << O2.Get_etrans() << ' '
               << O2.Get_eint() << ' '
               << O2.Get_etotal() << ' '

               << O2.Get_Cv_trans() << ' '
               << O2.Get_Cv_int() << ' '
               << O2.Get_Cv_total() << '\n';
        /*
        cout << "Z_trans = " << O2.Get_Ztrans() << ' '
             << "Z_int   = " << O2.Get_Zint() << ' '
             << "Z_total = " << O2.Get_Ztotal() << '\n'

             << "e_trans = " << O2.Get_etrans() << ' '
             << "e_int   = " << O2.Get_eint() << ' '
             << "e_total = " << O2.Get_etotal() << '\n'

             << "Cv_trans = " << O2.Get_Cv_trans() << ' '
             << "Cv_int   = " << O2.Get_Cv_int() << ' '
             << "Cv_total = " << O2.Get_Cv_total() << '\n';
             */
    }

}
