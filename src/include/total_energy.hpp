// [class TotalEnergy]
// total energy, ewald energy, etc.

#ifndef TC_TOTAL_ENERGY_HPP
#define TC_TOTAL_ENERGY_HPP

class TotalEnergy
{
private:
    // energy_*body[spin][kpoint][band]. used for total-energy calculation
    std::vector<std::vector<std::vector<Complex> > > energy_1body_; // <phi|H_1body|phi> (for BITC, <phi_left|H_1body|phi>)
    std::vector<std::vector<std::vector<Complex> > > energy_2body_;
    std::vector<std::vector<std::vector<Complex> > > energy_3body_;

    Complex total_energy_; 
    Complex total_energy_difference_;
    Complex total_energy_1body_;
    Complex total_energy_2body_;
    Complex total_energy_3body_;

    double ewald_energy_;

    void resize_energy(const int &num_independent_spins,
                       const int &num_irreducible_kpoints_scf,
                       const std::vector<int> &num_bands_scf,
                       const std::string &calc_method);
public:
    const std::vector<std::vector<std::vector<Complex> > > &energy_1body() const { return energy_1body_; }
    const std::vector<std::vector<std::vector<Complex> > > &energy_2body() const { return energy_2body_; }
    const std::vector<std::vector<std::vector<Complex> > > &energy_3body() const { return energy_3body_; }
    Complex total_energy() const { return total_energy_; }
    Complex total_energy_difference() const { return total_energy_difference_; }
    Complex total_energy_1body() const { return total_energy_1body_; }
    Complex total_energy_2body() const { return total_energy_2body_; }
    Complex total_energy_3body() const { return total_energy_3body_; }
    double ewald_energy() const { return ewald_energy_; }

    TotalEnergy() : total_energy_(0.0), total_energy_difference_(0.0), total_energy_1body_(0.0),
                    total_energy_2body_(0,0), total_energy_3body_(0,0), ewald_energy_(0.0) {}

    // initialization
    void set_ewald_energy(const bool is_heg, const double &ewald_energy);
    void bcast_qe_input(const int &num_independent_spins,
                        const int &num_irreducible_kpoints,
                        const std::vector<int> &num_bands_scf,
                        const std::string &calc_mode, const std::string &calc_method,
                        const bool is_heg);

    void reset_energies(const int &ispin, const int &ik,
                        const Eigen::MatrixXcd &energy_1body,
                        const Eigen::MatrixXcd &energy_2body,
                        const Eigen::MatrixXcd &energy_3body,
                        const bool uses_3body);
    void bcast_energies(const bool uses_3body);

    void calc_total_energy(const Spin &spin,
                           const Kpoints &kpoints,
                           const std::vector<std::vector<std::vector<double> > > &filling,
                           const std::vector<std::vector<std::vector<Complex> > > &eigenvalues_scf,
                           const double &fermi_energy,
                           const bool uses_3body,
                           const bool am_i_mpi_rank0,
                           std::ostream *ost);
};


#endif // TC_TOTAL_ENERGY_HPP
