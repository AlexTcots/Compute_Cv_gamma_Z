cmake_minimum_required(VERSION 3.20)
project(Compute_Cv_gamma_Z)

set(CMAKE_CXX_STANDARD 17)

add_executable(Compute_Cv_gamma_Z main.cpp Spectro_Atom.cpp Spectro_Atom.h Constants.h Spectro_molecule.cpp Spectro_molecule.h Tests.cpp Tests.h Particle.cpp Particle.h Atom.cpp Atom.h Molecule.cpp Molecule.h computing.cpp computing.h)
