//
// Created by alext on 4/21/2022.
//

#include "Spectro_Atom.h"

Spectro_Atom::Spectro_Atom(int n_max, std::string file_name) {
    data.reserve(n_max);
    std::ifstream ifs(file_name);
    if (!ifs.is_open()) {
        throw std::invalid_argument("fail to open file " + file_name);
    } else {
        std::string line;
        for (int i = 0; i < n_max; ++i) {
            getline(ifs, line);
            std::istringstream iss(line);
            double gn, level;
            iss >> gn >> level;
            data.push_back({gn, level});
        }
    }
}

atom_data Spectro_Atom::Get_data_at(int n) const {
    return data.at(n);
}

const size_t Spectro_Atom::Get_data_size() const {
    return data.size();
}

Spectro_Atom::Spectro_Atom(const Spectro_Atom &other) {
    data=other.data;
}

