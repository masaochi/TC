Input parameters in input.in
============================

Input parameters that should be specified in ``input.in`` are listed below in the following format:

+-------------------+----------------------+
| **Keyword**       |   TYPE               |
+-------------------+----------------------+
|   Default value   |   Available values   |
+-------------------+----------------------+
|   Detailed description                   |
+------------------------------------------+

All inputs are case-sensitive. Keywords listed in the alphabetical order:

[**Mandatory keywords**] 
calc_method_ , calc_mode_ , pseudo_dir_ , qe_save_dir_

[**Optional keywords**]
A_dn_dn_ , A_up_dn_ , A_up_up_ , charge_tolerance_ , energy_tolerance_ , force_tolerance_ , includes_div_correction_ ,
is_heg_ , max_num_blocks_david_ , max_num_ionic_steps_ , max_num_iterations_ , mixes_density_matrix_ , mixing_beta_ , num_bands_tc_ ,
num_refresh_david_ , reads_crystal_structure_ , restarts_ , smearing_mode_ , smearing_width_ 

Mandatory Keywords
------------------

.. _calc_method:

+------------------------------------+------------------------------------------------------+
| **calc_method**                    | STRING                                               |
+------------------------------------+------------------------------------------------------+
| **Mandatory** (*no default value*) | FREE, HF, TC, BITC                                   |
+------------------------------------+------------------------------------------------------+
| | Calculation method: free-electron mode (FREE), Hartree-Fock (HF), transcorrelated (TC), |
| | bi-orthogonal transcorrelated (BITC) methods. No electron-electron interaction is       |
| | considered for FREE, i.e., the kinetic energy and pseudopotentials are only considered. |
+-------------------------------------------------------------------------------------------+

.. _calc_mode:

+------------------------------------+-----------------------------------------------------+
| **calc_mode**                      | STRING                                              |
+------------------------------------+-----------------------------------------------------+
| **Mandatory** (*no default value*) | SCF, BAND                                           |
+------------------------------------+-----------------------------------------------------+
| Calculation mode. BAND calculation should be performed after SCF calculation.            |
+------------------------------------------------------------------------------------------+

.. _pseudo_dir:

+------------------------------------+-------------------------------------------------------------+
| **pseudo_dir**                     | STRING                                                      |
+------------------------------------+-------------------------------------------------------------+
| **Mandatory** (*no default value*) | *any*                                                       |
+------------------------------------+-------------------------------------------------------------+
| A directory where pseudopotential files are placed, e.g., /home/user/where_pseudopot_are_placed  |
+--------------------------------------------------------------------------------------------------+

.. _qe_save_dir:

+------------------------------------+---------------------------------------------------------+
| **qe_save_dir**                    | STRING                                                  |
+------------------------------------+---------------------------------------------------------+
| **Mandatory** (*no default value*) | *any*                                                   |
+------------------------------------+---------------------------------------------------------+
| A ``save`` directory created by QE, e.g., /home/user/where_QEcalc_was_performed/prefix.save  |
+----------------------------------------------------------------------------------------------+

Optional keywords
-----------------

.. _A_dn_dn:

+------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| **A_dn_dn**                        | REAL                                                                                                                                                                                 |
+------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| *Default:* 1.0                     | *any*                                                                                                                                                                                |
+------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| | :math:`A_{\downarrow,\downarrow}` in the RPA-type Jastrow function, :math:`\frac{A_{\sigma, \sigma'}}{|{\bf r}-{\bf r'}|}(1-e^{-|{\bf r}-{\bf r'}|/C_{\sigma,\sigma'}})`, normalized by :math:`A_0 = 1/\sqrt{4 \pi n}`, |
| | the value for the homogeneous electron gas having the same electron density :math:`n` of the system.                                                                                                                    |
| | For example, **A_dn_dn** = 0.8 means :math:`A_{\downarrow,\downarrow}=0.8A_0`.                                                                                                                                          | 
| | :math:`C_{\downarrow,\downarrow}` is determined so as to satisfy the cusp condition. Not used for calc_method_ = FREE or HF.                                                                                            |
+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

.. _A_up_dn:

+------------------------------------+----------------------------------------------------------------------------------------+
| **A_up_dn**                        | REAL                                                                                   |
+------------------------------------+----------------------------------------------------------------------------------------+
| *Default:* 1.0                     | *any*                                                                                  |
+------------------------------------+----------------------------------------------------------------------------------------+
| :math:`A_{\uparrow,\downarrow}`, same as A_dn_dn_.    :math:`A_{\uparrow,\downarrow}=A_{\downarrow,\uparrow}` is assumed.   |
+-----------------------------------------------------------------------------------------------------------------------------+

.. _A_up_up:

+------------------------------------+-----------------------------------------------------------------------------+
| **A_up_up**                        | REAL                                                                        |
+------------------------------------+-----------------------------------------------------------------------------+
| *Default:* 1.0                     | *any*                                                                       |
+------------------------------------+-----------------------------------------------------------------------------+
| | :math:`A_{\uparrow,\uparrow}`, same as A_dn_dn_.    Users cannot specify different values for **A_up_up** and  |
| | A_dn_dn_ in non-spin-polarized calculation.                                                                    |
+------------------------------------------------------------------------------------------------------------------+

.. _charge_tolerance:

+------------------------------------------------+-----------------------------------------------------------------+
| **charge_tolerance**                           | REAL                                                            |
+------------------------------------------------+-----------------------------------------------------------------+
| *Default:* 1e-5 (1e-4 for ver.1.2 and before)  | :math:`\geq 0`                                                  |
+------------------------------------------------+-----------------------------------------------------------------+
| In :math:`e^-`. Convergence criteria for the charge density, used only for calc_mode_ = SCF.                     |
+------------------------------------------------------------------------------------------------------------------+

.. _energy_tolerance:

+-----------------------------------------------+--------------------------------------------------------+
| **energy_tolerance**                          | REAL                                                   |
+-----------------------------------------------+--------------------------------------------------------+
| *Default:* 1e-6 (1e-5 for ver.1.2 and before) | :math:`\geq 0`                                         |
+-----------------------------------------------+--------------------------------------------------------+
| In Ht. Convergence criteria for the total energy (calc_mode_ = SCF) or a sum of eigenvalues (BAND).    |
+--------------------------------------------------------------------------------------------------------+

.. _force_tolerance:

+------------------------------------+-------------------------------------------------------------------+
| **force_tolerance** (from ver.1.3) | REAL                                                              |
+------------------------------------+-------------------------------------------------------------------+
| *Default:* 1e-2                    | :math:`\geq 0`                                                    |
+------------------------------------+-------------------------------------------------------------------+
| In eV/ang. Convergence criteria for the force in structural optimization (see max_num_ionic_steps_)    |
+--------------------------------------------------------------------------------------------------------+

.. _includes_div_correction:

+------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| **includes_div_correction**        | BOOLEAN                                                                                                                |
+------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| *Default:* true                    | true, false                                                                                                            |
+------------------------------------+------------------------------------------------------------------------------------------------------------------------+
| | Whether the divergence correction in the reciprocal space is included. See our paper:                                                                     |
| | M. Ochi, `Comput. Phys. Commun. 287, 108687 (2023) <https://doi.org/10.1016/j.cpc.2023.108687>`_ . [`arXiv <http://arxiv.org/abs/2302.07420>`_]           |
+-------------------------------------------------------------------------------------------------------------------------------------------------------------+

.. _is_heg:

+------------------------------------+---------------------------------------------------------------------+
| **is_heg**                         | BOOLEAN                                                             |
+------------------------------------+---------------------------------------------------------------------+
| *Default:* false                   | true, false                                                         |
+------------------------------------+---------------------------------------------------------------------+
| | Switches on the homogeneous-electron-gas mode where pseudopotentials and the Ewald energy              |
| | are ignored. (i.e., energy contribution from a lattice is ignored.)                                    |
+----------------------------------------------------------------------------------------------------------+

.. _max_num_blocks_david:

+------------------------------------+---------------------------------------------------------------------+
| **max_num_blocks_david**           | INTEGER                                                             |
+------------------------------------+---------------------------------------------------------------------+
| *Default:* 2                       | :math:`\geq 2`                                                      |
+------------------------------------+---------------------------------------------------------------------+
| | This keywords determines a size of subspace dimension for diagonalization.                             |
| | subspace dimension = **max_num_blocks_david** :math:`\times` num_bands_tc_.                            |
| | Increasing this value can improve convergence while computational time is proportional to it.          |
+----------------------------------------------------------------------------------------------------------+

.. _max_num_iterations:

+---------------------------------------------------------------+------------------------------------------+
| **max_num_iterations**                                        | INTEGER                                  |
+---------------------------------------------------------------+------------------------------------------+
| *Default:* 30 for calc_mode_ = SCF, 15 for calc_mode_ = BAND  | :math:`\geq 0`                           |
+---------------------------------------------------------------+------------------------------------------+
| Maximum number of iterations for the self-consistent-field loop (Also needed for calc_mode_ = BAND).     |
+----------------------------------------------------------------------------------------------------------+

.. _max_num_ionic_steps:

+---------------------------------------------------------------+------------------------------------------+
| **max_num_ionic_steps** (from ver.1.3)                        | INTEGER                                  |
+---------------------------------------------------------------+------------------------------------------+
| *Default:* 0                                                  | :math:`\geq 0`                           |
+---------------------------------------------------------------+------------------------------------------+
| | Maximum number of structural-optimization steps. Structural optimization is switched on when           |
| | calc_mode_ = SCF and is_heg_ = false and max_num_ionic_steps_ > 0 and calc_method_ = HF or BITC.       |
+----------------------------------------------------------------------------------------------------------+

.. _mixes_density_matrix:

+------------------------------------+--------------------------------------------------------------------------+
| **mixes_density_matrix**           | BOOLEAN                                                                  |
+------------------------------------+--------------------------------------------------------------------------+
| *Default:* false                   | true, false                                                              |
+------------------------------------+--------------------------------------------------------------------------+
| | The density matrix (true) or the density (false) is used for mixing. The former with a small mixing_beta_   |
| | can improve the convergence of calculation while computational time will be longer (by about a              |
| | factor of two). Used only for calc_mode_ = SCF. From ver.1.1.                                               |
+---------------------------------------------------------------------------------------------------------------+

.. _mixing_beta:

+------------------------------------+--------------------------------------------------------------------------------------+
| **mixing_beta**                    | REAL                                                                                 |
+------------------------------------+--------------------------------------------------------------------------------------+
| *Default:* 0.7                     | :math:`> 0`                                                                          |
+------------------------------------+--------------------------------------------------------------------------------------+
| | Mixing ratio for linear mixing: new density (or density matrix) = **mixing_beta** :math:`\times` new density :math:`+`  |
| | :math:`(1-` **mixing_beta** :math:`)\times` old density. Used only for calc_mode_ = SCF.                                |
+---------------------------------------------------------------------------------------------------------------------------+

.. _num_bands_tc:

+------------------------------------+----------------------------------------------------------------------------------------+
| **num_bands_tc**                   | INTEGER                                                                                |
+------------------------------------+----------------------------------------------------------------------------------------+
| *Default:* **nbnd** in QE          | :math:`1 \leq` **num_bands_tc** :math:`\leq` **nbnd** in QE                            |
+------------------------------------+----------------------------------------------------------------------------------------+
|   The number of bands for calculating eigenvalues. Also see max_num_blocks_david_.                                          |
+-----------------------------------------------------------------------------------------------------------------------------+

.. _num_refresh_david:

+------------------------------------+--------------------------------------------------------------------------+
| **num_refresh_david**              | INTEGER                                                                  |
+------------------------------------+--------------------------------------------------------------------------+
| *Default:* 1                       | :math:`\geq 1`                                                           |
+------------------------------------+--------------------------------------------------------------------------+
| Trial vectors are updated by **num_refresh_david** times for each update of the Fock operator.                |
+---------------------------------------------------------------------------------------------------------------+

.. _reads_crystal_structure:

+---------------------------------------------+------------------------------------------------------------+
| **reads_crystal_structure** (from ver.1.3)  | BOOLEAN                                                    |
+---------------------------------------------+------------------------------------------------------------+
| *Default:* false                            | true, false                                                |
+---------------------------------------------+------------------------------------------------------------+
| | When **reads_crystal_structure** = true, TC++ reads ``tc_crystal_structure.dat``,                      |
| | which is dumped in a previous structural-optimization run.                                             |
| | Note that this option is available only for the structural optimization (i.e. max_num_ionic_steps_ >0).|
| | Please also see :doc:`how_to_use` describing ``tc_crystsal_structure.dat``.                            |
+----------------------------------------------------------------------------------------------------------+

.. _restarts:

+------------------------------------+---------------------------------------------------------------------+
| **restarts**                       | BOOLEAN                                                             |
+------------------------------------+---------------------------------------------------------------------+
| *Default:* false                   | true, false                                                         |
+------------------------------------+---------------------------------------------------------------------+
| When **restarts** = true, TC++ restarts calculation from a previous run.                                 |
+----------------------------------------------------------------------------------------------------------+

.. _smearing_mode:

+------------------------------------+---------------------------------------------------------------------+
| **smearing_mode**                  | STRING                                                              |
+------------------------------------+---------------------------------------------------------------------+
| *Default:* gaussian                | fixed, gaussian                                                     |
+------------------------------------+---------------------------------------------------------------------+
| | *fixed*: fixed occupation for each k-point, *gaussian*: Gaussian smearing with smearing_width_.        |
| | Recommended values are *fixed* for insulators and *gaussian* for metals.                               |
+----------------------------------------------------------------------------------------------------------+


.. _smearing_width:

+------------------------------------+--------------------------------------------------------------------------+
| **smearing_width**                 | REAL                                                                     |
+------------------------------------+--------------------------------------------------------------------------+
| *Default:* 0.01                    | :math:`\geq 0` (A negative value will be ignored.)                       |
+------------------------------------+--------------------------------------------------------------------------+
| In Ht. Not used for smearing_mode_ = fixed.                                                                   |
+---------------------------------------------------------------------------------------------------------------+

Examples of input.in
--------------------

Example 1 (Minimum ``input.in`` for insulators)

::
   
   calc_method    TC
   calc_mode      SCF
   pseudo_dir     /home/user/pseudopot
   qe_save_dir    /home/user/QE/Si/prefix.save
   smearing_mode  fixed

Example 2 (Minimum ``input.in`` for metals)

::
   
   calc_method     TC
   calc_mode       SCF
   pseudo_dir      /home/user/pseudopot
   qe_save_dir     /home/user/QE/Al/prefix.save
   smearing_mode   gaussian
   smearing_width  0.01

Example 3 (Restart calculation after Example 1)

::
   
   calc_method     TC
   calc_mode       SCF
   pseudo_dir      /home/user/pseudopot
   qe_save_dir     /home/user/QE/Si/prefix.save
   smearing_mode   fixed
   restarts        true

Example 4

::
   
   calc_method     BITC
   calc_mode       SCF
   pseudo_dir      /home/user/pseudopot
   qe_save_dir     /home/user/QE/something/prefix.save
   smearing_mode   fixed
   A_up_up         0.2
   A_up_dn         0.2
   A_dn_dn         0.2
   max_num_iterations    15
   max_num_blocks_david  5
   mixes_density_matrix  true
   mixing_beta     0.2

Example 5 (Structural optimization)

::
   
   calc_method     BITC
   calc_mode       SCF
   pseudo_dir      /home/user/pseudopot
   qe_save_dir     /home/user/QE/something/prefix.save
   smearing_mode   fixed
   max_num_iterations    20
   max_num_ionic_steps   10   

Example 6 (Restarted structural optimization reading ``tc_crystal_structure.dat``)

::
   
   calc_method     BITC
   calc_mode       SCF
   pseudo_dir      /home/user/pseudopot
   qe_save_dir     /home/user/QE/something/prefix.save
   smearing_mode   fixed
   max_num_iterations       20
   max_num_ionic_steps      10
   reads_crystal_structure  true
