// [namespace calc_hamiltonian]
// calculate the matrix elements of Hamiltonian

#ifndef TC_CALC_HAMILTONIAN_HPP
#define TC_CALC_HAMILTONIAN_HPP

namespace calc_hamiltonian
{

void all(const Parallelization &parallelization, 
         MyClock &my_clock, const Method &method,
         const CrystalStructure &crystal_structure,
         const Symmetry &symmetry, const Potentials &potentials,
         const Spin &spin, const Kpoints &kpoints,
         PlaneWaveBasis &plane_wave_basis, 
         const BlochStates &bloch_states, 
         const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
         std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H1phi,
         std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H2phi,
         std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H3phi,
         std::ostream *ost);

void kinetic(const Parallelization &parallelization, 
             const Method &method,
             const CrystalStructure &crystal_structure,
             const Spin &spin, const Kpoints &kpoints,
             PlaneWaveBasis &plane_wave_basis, 
             const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
             std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H1phi,
             std::ostream *ost);

void pseudo(const Parallelization &parallelization, 
            const Method &method,
            const CrystalStructure &crystal_structure,
            const Potentials &potentials,
            const Spin &spin, const Kpoints &kpoints,
            PlaneWaveBasis &plane_wave_basis, 
            const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
            std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H1phi,
            std::ostream *ost);

// two-body, Hartree (Hartree-Fock)
void hf2h(const Parallelization &parallelization, 
          const Method &method,
          const CrystalStructure &crystal_structure,
          const Potentials &potentials,
          const Spin &spin, const Kpoints &kpoints,
          PlaneWaveBasis &plane_wave_basis, 
          const BlochStates &bloch_states, 
          const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
          std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H2phi,
          std::ostream *ost);
// two-body, exchange (Hartree-Fock)
void hf2x(const Parallelization &parallelization, 
          const Method &method,
          const CrystalStructure &crystal_structure,
          const Potentials &potentials,
          const Spin &spin, const Kpoints &kpoints,
          PlaneWaveBasis &plane_wave_basis, 
          const BlochStates &bloch_states, 
          const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
          std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H2phi,
          std::ostream *ost);

// two-body, Hartree (TC or BITC)
void tc2h(const Parallelization &parallelization, 
          const Method &method,
          const CrystalStructure &crystal_structure,
          const Potentials &potentials,
          const Spin &spin, const Kpoints &kpoints,
          PlaneWaveBasis &plane_wave_basis, 
          const BlochStates &bloch_states, 
          const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
          std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H2phi,
          std::ostream *ost);
// two-body, exchange (TC or BITC)
void tc2x(const Parallelization &parallelization, 
          const Method &method,
          const CrystalStructure &crystal_structure,
          const Potentials &potentials,
          const Spin &spin, const Kpoints &kpoints,
          PlaneWaveBasis &plane_wave_basis, 
          const BlochStates &bloch_states, 
          const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
          std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H2phi,
          std::ostream *ost);

// three-body (TC or BITC)
void tc3a1(const Parallelization &parallelization, 
           const Method &method,
           const CrystalStructure &crystal_structure,
           const Potentials &potentials,
           const Spin &spin, const Kpoints &kpoints,
           PlaneWaveBasis &plane_wave_basis, 
           const BlochStates &bloch_states, 
           const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
           std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H3phi,
           std::ostream *ost);
void tc3b1(const Parallelization &parallelization, 
           const Method &method,
           const CrystalStructure &crystal_structure,
           const Potentials &potentials,
           const Spin &spin, const Kpoints &kpoints,
           PlaneWaveBasis &plane_wave_basis, 
           const BlochStates &bloch_states, 
           const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
           std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H3phi,
           std::ostream *ost);

// merged calculation (for reducing computational cost)
void tc3a3b3b4b6(const Parallelization &parallelization, 
                 const Method &method,
                 const CrystalStructure &crystal_structure,
                 const Potentials &potentials,
                 const Spin &spin, const Kpoints &kpoints,
                 PlaneWaveBasis &plane_wave_basis, 
                 const BlochStates &bloch_states, 
                 const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
                 std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H3phi,
                 std::ostream *ost);
void tc3a2a4b2b5(const Parallelization &parallelization, 
                 const Method &method,
                 const CrystalStructure &crystal_structure,
                 const Symmetry &symmetry, const Potentials &potentials,
                 const Spin &spin, const Kpoints &kpoints,
                 PlaneWaveBasis &plane_wave_basis, 
                 const BlochStates &bloch_states, 
                 const std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &phi,
                 std::vector<std::vector<std::vector<std::vector<Eigen::VectorXcd> > > > &H3phi,
                 std::ostream *ost);
    
std::vector<Eigen::Vector3cd> force(const Parallelization &parallelization, 
                                    const Method &method,
                                    const CrystalStructure &crystal_structure,
                                    const Potentials &potentials,
                                    const Spin &spin, const Kpoints &kpoints,
                                    PlaneWaveBasis &plane_wave_basis,
                                    const BlochStates &bloch_states,
                                    const TotalEnergy &total_energy, 
                                    std::ostream *ost);
} // namespace calc_hamiltonian

#endif // TC_CALC_HAMILTONIAN_HPP
