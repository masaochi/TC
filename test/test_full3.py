import subprocess


def copy_file(filename):
    try:
        subprocess.run('cp ' + filename + ' .', shell=True, check=True)
    except subprocess.CalledProcessError:
        print('\033[31m'+'Error'+'\033[0m'+': Failed in preparation')
        print('Make sure that ' + filename + ' exists.')
        exit()


def copy_scfoutput(foldername):
    subprocess.run('cp tc_energy_scf.dat ' + foldername + '/.', shell=True, check=True)
    subprocess.run('cp tc_wfc_scf.dat ' + foldername + '/.', shell=True, check=True)
    subprocess.run('cp tc_scfinfo.dat ' + foldername + '/.', shell=True, check=True)


def copy_bandoutput(foldername):
    subprocess.run('cp tc_energy_scf.dat ' + foldername + '/.', shell=True, check=True)
    subprocess.run('cp tc_wfc_scf.dat ' + foldername + '/.', shell=True, check=True)
    subprocess.run('cp tc_energy_band.dat ' + foldername + '/.', shell=True, check=True)
    subprocess.run('cp tc_wfc_band.dat ' + foldername + '/.', shell=True, check=True)
    subprocess.run('cp tc_scfinfo.dat ' + foldername + '/.', shell=True, check=True)


def check_tc(output, output_ref, checks_totE):
    try:
        subprocess.run('./tc++', shell=True, check=True)
    except subprocess.CalledProcessError:
        print('Error: Failed to run ./tc++')
        print('See output.out (if exists), make sure that tc++ is present in the current directory (i.e., test/tc++).')
        exit()

    with open(output) as f:
        lines = f.readlines()

    with open(output_ref) as f_ref:
        lines_ref = f_ref.readlines()

    # Total energy check
    if checks_totE:
        lines_ref_totE = [line for line in lines_ref if 'Total energy =' in line]
        lines_ref_totE_split = lines_ref_totE[len(lines_ref_totE)-1].split()
        totE_ref = [lines_ref_totE_split[4]]
        lines_ref_totE_split = lines_ref_totE[len(lines_ref_totE)-2].split()
        totE_ref.append(lines_ref_totE_split[4])
        
        print('  Reference total energies (two SCF loops)', totE_ref)
        
        lines_totE = [line for line in lines if 'Total energy =' in line]
        lines_totE_split = lines_totE[len(lines_totE)-1].split()
        totE = [lines_totE_split[4]]
        lines_totE_split = lines_totE[len(lines_totE)-2].split()
        totE.append(lines_totE_split[4])

        print('  Calculated total energies (two SCF loops)', totE)
        
        if (abs(float(totE[0]) - float(totE_ref[0])) < 1e-4 and abs(float(totE[1]) - float(totE_ref[1])) < 1e-4):
            print('\033[32m'+'  Total-energy check passed.'+'\033[0m')
        else:
            print('\033[31m'+'Error'+'\033[0m'+': Total-energy difference is larger than 1e-4')
            exit()

    # Orbital energy check
    lines_ref_orbE = [line for line in lines_ref if '   6 ' in line]
    lines_ref_orbE_split = lines_ref_orbE[len(lines_ref_orbE)-1].split() # final line
    lines_ref_orbE_split2 = lines_ref_orbE[len(lines_ref_orbE)-2].split()
    orbE_ref = [lines_ref_orbE_split[1], lines_ref_orbE_split2[1]]
    print('  Reference orbital energies', orbE_ref)

    lines_orbE = [line for line in lines if '   6 ' in line]
    lines_orbE_split = lines_orbE[len(lines_orbE)-1].split() # final line
    lines_orbE_split2 = lines_orbE[len(lines_orbE)-2].split()
    orbE = [lines_orbE_split[1], lines_orbE_split2[1]]
    print('  Calculated orbital energies', orbE)

    if (abs(float(orbE[0]) - float(orbE_ref[0])) < 1e-4 and abs(float(orbE[1]) - float(orbE_ref[1])) < 1e-4):
        print('\033[32m'+'  Orbital-energy check passed.'+'\033[0m')
    else:
        print('\033[31m'+'Error'+'\033[0m'+': Orbital-energy difference is larger than 1e-4')
        exit()


def remove_temporary_files():
    print('Remove temporary files')
    subprocess.run('rm ./tc_energy_scf.dat', shell=True)
    subprocess.run('rm example*/tc_energy_scf.dat', shell=True)
    subprocess.run('rm ./tc_wfc_scf.dat', shell=True)
    subprocess.run('rm example*/tc_wfc_scf.dat', shell=True)
    subprocess.run('rm ./tc_energy_band.dat', shell=True)
    subprocess.run('rm ./tc_wfc_band.dat', shell=True)
    subprocess.run('rm ./tc_scfinfo.dat', shell=True)
    subprocess.run('rm example*/tc_scfinfo.dat', shell=True)
    subprocess.run('rm ./tc_bandplot.dat', shell=True)
    subprocess.run('rm ./tc_bandplot_up.dat', shell=True)
    subprocess.run('rm ./tc_bandplot_dn.dat', shell=True)
    subprocess.run('rm ./input.in', shell=True)
    subprocess.run('rm ./output.out', shell=True)
    subprocess.run('rm ./pwfn.data', shell=True)
    subprocess.run('rm ./parameters.casl', shell=True)
    subprocess.run('rm ./parameters.casl.dump', shell=True)
    subprocess.run('rm ./jastrow.plt', shell=True)
    print('\n')


print('\033[35m'+'Start test calculations.'+'\033[0m'+'\n')

print('\033[31m'+'Note! Full test will take > 30 minutes in total, rather for developers.')
print('If you would like to perform a minimal check of installation, stop this test and instead type'+'\033[0m')
print('  python3 test.py'+'\n')

print('\033[31m'+'Note! DO NOT use the input files (including pseudopot.) provided here for your research.')
print('These were made only for test calculation...'+'\033[0m'+'\n')

#remove_temporary_files()

# common in the following examples
output = './output.out'
num_test = str(33)

print('Test 1-25 were already done in test_full1.py & test_full2.py...'+'\n')

## Example 26
print('\033[35m'+'Test 26/'+num_test+'\033[0m'+': HF SCF calculation for bulk Si (gamma-point calc.)')
copy_file('./example26/input.in')
check_tc(output, './example26/output.out', True)
print('\033[35m'+'Test 26/'+num_test+'\033[32m'+' passed.'+'\033[0m'+'\n')

## Example 27
print('\033[35m'+'Test 27/'+num_test+'\033[0m'+': HF SCF calculation for bulk Si (gamma-point spin-pol. calc.)')
copy_file('./example27/input.in')
check_tc(output, './example27/output.out', True)
print('\033[35m'+'Test 27/'+num_test+'\033[32m'+' passed.'+'\033[0m'+'\n')

## Example 28
print('\033[35m'+'Test 28/'+num_test+'\033[0m'+': TC SCF calculation for bulk Si (poly Jastrow)')
copy_file('./example28/input.in')
copy_file('./example28/parameters.casl')
check_tc(output, './example28/output.out', True)
print('\033[35m'+'Test 28/'+num_test+'\033[32m'+' passed.'+'\033[0m'+'\n')

# copy output files for the subsequent tests...
copy_scfoutput('./example29')

## Example 29
print('\033[35m'+'Test 29/'+num_test+'\033[0m'+': TC BAND calculation for bulk Si (poly Jastrow)')
copy_file('./example29/input.in')
copy_file('./example29/tc_energy_scf.dat')
copy_file('./example29/tc_wfc_scf.dat')
copy_file('./example29/tc_scfinfo.dat')
check_tc(output, './example29/output.out', False)
print('\033[35m'+'Test 29/'+num_test+'\033[32m'+' passed.'+'\033[0m'+'\n')

## Example 30
print('\033[35m'+'Test 30/'+num_test+'\033[0m'+': TC SCF calculation for bulk Si (RPA+poly Jastrow)')
copy_file('./example30/input.in')
copy_file('./example30/parameters.casl')
check_tc(output, './example30/output.out', True)
print('\033[35m'+'Test 30/'+num_test+'\033[32m'+' passed.'+'\033[0m'+'\n')

# copy output files for the subsequent tests...
copy_scfoutput('./example31')

## Example 31
print('\033[35m'+'Test 31/'+num_test+'\033[0m'+': TC BAND calculation for bulk Si (RPA+poly Jastrow)')
copy_file('./example31/input.in')
copy_file('./example31/tc_energy_scf.dat')
copy_file('./example31/tc_wfc_scf.dat')
copy_file('./example31/tc_scfinfo.dat')
check_tc(output, './example31/output.out', False)
print('\033[35m'+'Test 31/'+num_test+'\033[32m'+' passed.'+'\033[0m'+'\n')

## Example 32
print('\033[35m'+'Test 32/'+num_test+'\033[0m'+': TC SCF calculation for bulk Si (spin-pol. calc., RPA+poly Jastrow)')
copy_file('./example32/input.in')
copy_file('./example32/parameters.casl')
check_tc(output, './example32/output.out', True)
print('\033[35m'+'Test 32/'+num_test+'\033[32m'+' passed.'+'\033[0m'+'\n')

# copy output files for the subsequent tests...
copy_scfoutput('./example33')

## Example 33
print('\033[35m'+'Test 33/'+num_test+'\033[0m'+': TC BAND calculation for bulk Si (spin-pol. calc., RPA+poly Jastrow)')
copy_file('./example33/input.in')
copy_file('./example33/tc_energy_scf.dat')
copy_file('./example33/tc_wfc_scf.dat')
copy_file('./example33/tc_scfinfo.dat')
check_tc(output, './example33/output.out', False)
print('\033[35m'+'Test 33/'+num_test+'\033[32m'+' passed.'+'\033[0m'+'\n')

remove_temporary_files()

print('\033[35m'+'Since these tests are done in serial calculation, it is recommended to check MPI parallelization does not change the results by yourself...'+'\033[0m')

