// [namespace io_qmc_files]
// Read/write parameters.casl, write pwfn.data, write jastrow.plt

#include "include/header.hpp"

// called during initialization
void io_qmc_files::dump_jastrow_plt(const FileNames &file_names,
                                    const Jastrow &jastrow,
                                    const bool dump_down_down,
                                    std::ostream *ost)
{
    const std::string &tc_jastrow_plt = file_names.tc_jastrow_plt();

    const bool is_A_zero = jastrow.is_A_zero();
    const int deg_poly = jastrow.deg_poly();

    // setup
    *ost << " Dump jastrow.plt for plotting a Jastrow function using gnuplot (" << tc_jastrow_plt << ")" << std::endl;
    std::ofstream ofs(tc_jastrow_plt, std::ios::out);
    ofs.precision(15);
    if (ofs.fail()) { error_messages::cannot_open(tc_jastrow_plt); }

    // header
    ofs << "# Please type" << std::endl;
    ofs << "#   load 'jastrow.plt'" << std::endl;
    ofs << "# in gnuplot for plotting a Jastrow function u(r). (atomic unit)" << std::endl;
    ofs << std::endl;

    // several configs
    const double xmax = deg_poly==0 ? 10.0 :
        std::max({jastrow.L_poly()[0][0], jastrow.L_poly()[0][1], jastrow.L_poly()[1][0], jastrow.L_poly()[1][1]});
    ofs << "set xr [0:" << xmax << "]" << std::endl;
    ofs << "set xlabel \"r-r'\"" << std::endl;
    ofs << "set ylabel \"u(r-r') in F = exp(-u(r-r'))\"" << std::endl;
    ofs << "set samples 1000" << std::endl;
    ofs << std::endl;

    if (!is_A_zero)
    {
        ofs << "A_uu = " << jastrow.A_long()[0][0] << std::endl;
        ofs << "A_ud = " << jastrow.A_long()[0][1] << std::endl;
        if (dump_down_down) { ofs << "A_dd = " << jastrow.A_long()[1][1] << std::endl; }
        ofs << "F_uu = " << jastrow.F_long()[0][0] << std::endl;
        ofs << "F_ud = " << jastrow.F_long()[0][1] << std::endl;
        if (dump_down_down) { ofs << "F_dd = " << jastrow.F_long()[1][1] << std::endl; }
        ofs << "f(x,A,F) = (A/x)*(1 - exp(-x/F))" << std::endl;
        ofs << std::endl;
    }

    if (deg_poly>0)
    {
        for (int is1=0; is1<2; is1++) 
        {
            for (int is2=0; is2<2; is2++)
            {
                if (is1==0 && is2==0)
                {
                    ofs << "g_uu(x) = (x>";
                }
                else if (is1==0 && is2==1)
                {
                    ofs << "g_ud(x) = (x>";
                }
                else if (is1==1 && is2==1)
                {
                    if (!dump_down_down) { continue; } // skip down_down (= up_up for non-spin-polarized calc.)
                    ofs << "g_dd(x) = (x>";
                }
                else
                {
                    continue; // skip down_up (= up_down)
                }
                ofs << jastrow.L_poly()[is1][is2] << ") ? 0 : -(";
                for (int ipoly=0; ipoly<deg_poly; ipoly++)
                {
                    ofs << jastrow.c_poly_input()[is1][is2][ipoly] << "*x**" << ipoly;
                    if (ipoly==deg_poly-1)
                    {
                        ofs << ")";
                    }
                    else
                    {
                        ofs << " + ";
                    }
                }
                ofs << "*(1-x/" << jastrow.L_poly()[is1][is2] << ")**3" << std::endl;
            } // is2
        } // is1
        ofs << std::endl;
    }

    if (is_A_zero)
    {
        if (jastrow.deg_poly()==0) // no Jastrow
        {
            ofs << "p 0 ti 'u_up_up', 0 ti 'u_up_down'";
            if (dump_down_down) { ofs << ", 0 ti 'u_down_down'"; }
            ofs << std::endl;
        }
        else // poly only
        {
            ofs << "p g_uu(x) ti 'u_up_up', g_ud(x) ti 'u_up_down'";
            if (dump_down_down) { ofs << ", g_dd(x) ti 'u_down_down'"; }
            ofs << std::endl;
        }
    }
    else
    {
        if (deg_poly==0) // RPA only
        {
            ofs << "p f(x,A_uu,F_uu) ti 'u_up_up', f(x,A_ud,F_ud) ti 'u_up_down'";
            if (dump_down_down) { ofs << ", f(x,A_dd,F_dd) ti'u_down_down'"; }
            ofs  << std::endl;
        }
        else // RPA + poly
        {
            ofs << "p f(x,A_uu,F_uu) ti 'u_up_up_RPA', f(x,A_ud,F_ud) ti 'u_up_down_RPA', ";
            if (dump_down_down) { ofs << "f(x,A_dd,F_dd) ti'u_down_down_RPA', "; }
            ofs << "g_uu(x) ti 'u_up_up_poly', g_ud(x) ti 'u_up_down_poly', ";
            if (dump_down_down) { ofs << "g_dd(x) ti 'u_down_down_poly', "; }
            ofs << "f(x,A_uu,F_uu) + g_uu(x) ti 'u_up_up', f(x,A_ud,F_ud) + g_ud(x) ti 'u_up_down'";
            if (dump_down_down) { ofs << ", f(x,A_dd,F_dd) + g_dd(x) ti 'u_down_down'"; }
            ofs << std::endl;
        }
    }
    ofs.close();
    *ost << std::endl;
}

