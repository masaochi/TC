// [class AtomicSpecies]
// A list of atomic species and atomic numbers

#ifndef TC_ATOMIC_SPECIES_HPP
#define TC_ATOMIC_SPECIES_HPP

class AtomicSpecies
{
private:
    std::vector<std::string> names_;

public:
    int get_number(std::string name) const; // e.g., get_number(Si) = 14
    std::string get_name(const int number) const; // e.g., get_name(14) = Si

    AtomicSpecies() :
        names_{"H ", "He",
            "Li", "Be",
            "B ", "C ", "N ", "O ", "F ", "Ne",
            "Na", "Mg",
            "Al", "Si", "P ", "S ", "Cl", "Ar",
            "K ", "Ca",
            "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Ga", "Ge", "As", "Se", "Br", "Kr",
            "Rb", "Sr",
            "Y ", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
            "In", "Sn", "Sb", "Te", "I ", "Xe",
            "Cs", "Ba",
            "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", 
            "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
            "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg",
            "Tl", "Pb", "Bi", "Po", "At", "Rn",
            "Fr", "Ra",
            "Ac", "Th", "Pa", "U ", "Np", "Pu", "Am", 
            "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No"
            "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
            "Nh"} {}; // up to 113
};

#endif // TC_ATOMIC_SPECIES_HPP
