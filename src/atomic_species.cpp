// [class AtomicSpecies]
// A list of atomic species and atomic numbers

#include "include/header.hpp"

// e.g., get_number(Si) = 14
int AtomicSpecies::get_number(std::string name) const
{
    // make a length of "name" two
    const int length = name.length();
    if (length == 1) // e.g., name = "C"
    {
        name += " "; // "C" -> "C "
    }
    else if (length != 2)
    {
        error_messages::stop("Invalid name of element in <PP_HEADER> of a upf file");
    }

    // search
    auto itr = std::find(names_.begin(), names_.end(), name);
    if (itr == names_.end()) { error_messages::stop("Search failed in atomic_species::get_number()"); }

    return std::distance(names_.begin(), itr) + 1; // std::distance starts from 0 -> should +1
}

// e.g., get_name(14) = Si
std::string AtomicSpecies::get_name(const int number) const
{
    if (number < 1) { error_messages::stop("Inappropriate (not positive) atomic number"); }
    if (number > 113) { error_messages::stop("Too large atomic number > 113: Not supported"); }

    return names_[number-1]; // names_[0] = H, names_[1] = He,...
}

