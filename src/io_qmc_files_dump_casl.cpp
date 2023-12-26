// [namespace io_qmc_files]
// Read/write parameters.casl, write pwfn.data, write jastrow.plt

#include "include/header.hpp"

void io_qmc_files::dump_casl(const FileNames &file_names,
                             const Spin &spin,
                             const Jastrow &jastrow,
                             std::ostream *ost)
{
    *ost << " Dump the default Jastrow parameters to " << file_names.tc_casl_dump() << std::endl;
    std::ofstream ofs(file_names.tc_casl_dump(), std::ios::out);
    ofs.fixed;
    ofs.precision(15);
    ofs.showpoint;
    if (ofs.fail()) { error_messages::cannot_open(file_names.tc_casl_dump()); }

    const int num_independent_spins = spin.num_independent_spins();

    ofs << "JASTROW:" << std::endl;
    ofs << "  Title: solid" << std::endl;
    // RPA term
    ofs << "  TERM 1:" << std::endl;
    ofs << "    Rank: [ 2, 0 ]" << std::endl;
    if (num_independent_spins==1)
    {
        ofs << "    Rules: [ 1-1=2-2 ]" << std::endl;
    }
    else
    {
        ofs << "    Rules: [ ]" << std::endl;
    }
    ofs << "    e-e basis:" << std::endl;
    ofs << "      Type: RPA" << std::endl;
    ofs << "      Order: 1" << std::endl;
    ofs << "      Parameters:" << std::endl;
    ofs << "        Channel 1-1:" << std::endl;
    ofs << "          F: [ " << jastrow.F_long()[0][0]<< ", fixed, limits: [ 1.000000000000000E-007, +Inf" << std::endl;
    ofs << "               ] ]" << std::endl;
    ofs << "        Channel 1-2:" << std::endl;
    ofs << "          F: [ " << jastrow.F_long()[0][1]<< ", fixed, limits: [ 1.000000000000000E-007, +Inf" << std::endl;
    ofs << "               ] ]" << std::endl;
    if (num_independent_spins==2)
    {
        ofs << "        Channel 2-2:" << std::endl;
        ofs << "          F: [ " << jastrow.F_long()[1][1]<< ", fixed, limits: [ 1.000000000000000E-007, +Inf" << std::endl;
        ofs << "               ] ]" << std::endl;
    }
    ofs << "    e-e cutoff:" << std::endl;
    ofs << "      Type: spline" << std::endl; // spline cutoff is not used in TC++, just for CASINO calculation
    ofs << "      Constants:" << std::endl;
    ofs << "         C: 2" << std::endl;
    ofs << "      Parameters:" << std::endl;
    ofs << "        Channel 1-1:" << std::endl;
    ofs << "          L: [ 500.0, fixed, limits: [ 0.100000000000000, " << std::endl;
    ofs << "               500.0 ] ]" << std::endl;
    ofs << "          x: [ 0.500000000000000, fixed, limits: [ 5.000000000000000E-002," << std::endl;
    ofs << "               0.950000000000000 ] ]" << std::endl;
    ofs << "        Channel 1-2:" << std::endl;
    ofs << "          L: [ 500.0, fixed, limits: [ 0.100000000000000, " << std::endl;
    ofs << "               500.0 ] ]" << std::endl;
    ofs << "          x: [ 0.500000000000000, fixed, limits: [ 5.000000000000000E-002," << std::endl;
    ofs << "               0.950000000000000 ] ]" << std::endl;
    if (num_independent_spins==2)
    {
        ofs << "        Channel 2-2:" << std::endl;
        ofs << "          L: [ 500.0, fixed, limits: [ 0.100000000000000, " << std::endl;
        ofs << "               500.0 ] ]" << std::endl;
        ofs << "          x: [ 0.500000000000000, fixed, limits: [ 5.000000000000000E-002," << std::endl;
        ofs << "               0.950000000000000 ] ]" << std::endl;
    }
    ofs << "    e-e cusp: T" << std::endl;
    ofs << "    Linear parameters:" << std::endl;
    ofs << "      Channel 1-1: [ ]" << std::endl;
    ofs << "      Channel 1-2: [ ]" << std::endl;
    if (num_independent_spins==2)
    {
        ofs << "      Channel 2-2: [ ]" << std::endl;
    }
    // polynomial terms
    ofs << "  TERM 2:" << std::endl;
    ofs << "    Rank: [ 2, 0 ]" << std::endl;
    if (num_independent_spins==1)
    {
        ofs << "    Rules: [ 1-1=2-2 ]" << std::endl;
    }
    else
    {
        ofs << "    Rules: [ ]" << std::endl;
    }
    ofs << "    e-e basis:" << std::endl;
    ofs << "      Type: natural power" << std::endl;
    ofs << "      Order: 0" << std::endl;
    ofs << "    e-e cutoff:" << std::endl;
    ofs << "      Type: polynomial" << std::endl;
    ofs << "      Constants:" << std::endl;
    ofs << "        C: 3" << std::endl;
    ofs << "      Parameters:" << std::endl;
    ofs << "        Channel 1-1:" << std::endl;
    ofs << "          L: [ 500.0, fixed, limits: [ 0.100000000000000, " << std::endl;
    ofs << "               500.0 ] ]" << std::endl;
    ofs << "        Channel 1-2:" << std::endl;
    ofs << "          L: [ 500.0, fixed, limits: [ 0.100000000000000, " << std::endl;
    ofs << "               500.0 ] ]" << std::endl;
    if (num_independent_spins==2)
    {
        ofs << "        Channel 2-2:" << std::endl;
        ofs << "          L: [ 500.0, fixed, limits: [ 0.100000000000000, " << std::endl;
        ofs << "               500.0 ] ]" << std::endl;
    }
    ofs << "    e-e cusp: F" << std::endl;
    ofs << "    Linear parameters:" << std::endl;
    ofs << "      Channel 1-1: [ ]" << std::endl;
    ofs << "      Channel 1-2: [ ]" << std::endl;
    if (num_independent_spins==2)
    {
        ofs << "      Channel 2-2: [ ]" << std::endl;
    }
    ofs.close();
}
