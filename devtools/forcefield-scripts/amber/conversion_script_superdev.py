from __future__ import print_function
import parmed
from parmed.utils.six.moves import StringIO, zip
import simtk.openmm.app as app
from simtk.unit import *
import os
import re
from numpy import testing
import tempfile
import yaml
from distutils.spawn import find_executable
import hashlib
from collections import OrderedDict
import warnings
from parmed.exceptions import ParameterWarning
warnings.filterwarnings('error', category=ParameterWarning)

_loadoffre = re.compile(r'loadoff (\S*)', re.I)

def convert(filename, ignore=None, provenance=None):
    basename = os.path.basename(filename)
    if not os.path.exists('ffxml/'):
        os.mkdir('ffxml')
    ffxml_name = 'ffxml/' + basename + '.xml'
    print('Preparing %s for conversion...' % basename)
    with open(filename) as f:
        lines = map(lambda line:
                line if '#' not in line else line[:line.index('#')], f)
    if ignore is not None:
        new_lines = []
        for line in lines:
            if _loadoffre.findall(line) and _loadoffre.findall(line)[0] in ignore:
                continue
            else:
                new_lines.append(line)
    else:
        new_lines = lines
    leaprc = StringIO(''.join(new_lines))
    print('Converting to ffxml...')
    params = parmed.amber.AmberParameterSet.from_leaprc(leaprc)
    params = parmed.openmm.OpenMMParameterSet.from_parameterset(params)
    # If there are no residue templates (GAFF) - write_unused must be True
    if params.residues:
        params.write(ffxml_name, provenance=provenance, write_unused=False)
    else:
        params.write(ffxml_name, provenance=provenance, write_unused=True)
    print('Ffxml successfully written!')
    return ffxml_name

def validate_protein(ffxml_name, leaprc_name, united_atom=False):
    if not find_executable('tleap'):
        raise RuntimeError('tleap not available from PATH')
    print('Preparing temporary files for validation...')
    ala3_top = tempfile.mkstemp()
    ala3_crd = tempfile.mkstemp()
    villin_top = tempfile.mkstemp()
    villin_crd = tempfile.mkstemp()
    leap_script_ala3_file = tempfile.mkstemp()
    leap_script_villin_file = tempfile.mkstemp()

    print('Preparing LeaP scripts...')
    if not united_atom:
        leap_script_ala3_string = """source %s
x = loadPdb files/ala3.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, ala3_top[1], ala3_crd[1])
        leap_script_villin_string = """source %s
x = loadPdb files/villin.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, villin_top[1], villin_crd[1])
    else:
        leap_script_ala3_string = """source %s
x = loadPdb files/ala3_ua.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, ala3_top[1], ala3_crd[1])
        leap_script_villin_string = """source %s
x = loadPdb files/villin_ua.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, villin_top[1], villin_crd[1])

    os.write(leap_script_ala3_file[0], leap_script_ala3_string)
    os.write(leap_script_villin_file[0], leap_script_villin_string)

    print('Running LEaP...')
    os.system('tleap -f %s' % leap_script_ala3_file[1])
    if os.path.getsize(ala3_top[1]) == 0 or os.path.getsize(ala3_crd[1]) == 0:
        raise RuntimeError('Ala_ala_ala LEaP fail for %s' % leaprc_name)
    os.system('tleap -f %s' % leap_script_villin_file[1])
    if os.path.getsize(villin_top[1]) == 0 or os.path.getsize(villin_crd[1]) == 0:
        raise RuntimeError('Villin headpiece LEaP fail for %s' % leaprc_name)

    print('Calculating ala_ala_ala energies...')
    # AMBER
    parm_amber = parmed.load_file(ala3_top[1], ala3_crd[1])
    system_amber = parm_amber.createSystem(splitDihedrals=True)
    ala3_amber_energies = parmed.openmm.energy_decomposition_system(parm_amber, system_amber, nrg=kilojoules_per_mole)
    # OpenMM
    ff = app.ForceField(ffxml_name)
    system_omm = ff.createSystem(parm_amber.topology)
    parm_omm = parmed.openmm.load_topology(parm_amber.topology, system_omm, xyz=parm_amber.positions)
    system_omm = parm_omm.createSystem(splitDihedrals=True)
    ala3_omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=kilojoules_per_mole)

    print('Calculating villin headpiece energies...')
    # AMBER
    parm_amber = parmed.load_file(villin_top[1], villin_crd[1])
    system_amber = parm_amber.createSystem(splitDihedrals=True)
    villin_amber_energies = parmed.openmm.energy_decomposition_system(parm_amber, system_amber, nrg=kilojoules_per_mole)
    # OpenMM
    ff = app.ForceField(ffxml_name)
    system_omm = ff.createSystem(parm_amber.topology)
    parm_omm = parmed.openmm.load_topology(parm_amber.topology, system_omm, xyz=parm_amber.positions)
    system_omm = parm_omm.createSystem(splitDihedrals=True)
    villin_omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=kilojoules_per_mole)

    print('Deleting temp files...')
    for f in (ala3_top, ala3_crd, villin_top, villin_crd, leap_script_ala3_file,
             leap_script_villin_file):
        os.close(f[0])
        os.unlink(f[1])

    # temporary for testing
    #return (ala3_amber_energies, ala3_omm_energies, villin_amber_energies,
    #        villin_omm_energies)
    # debugging mod for turning dihedrals off
    e = ala3_omm_energies
    if len(e) == 4:
        ala3_omm_energies = [e[0], e[1], ('PeriodicTorsionForce', 0.0), ('PeriodicTorsionForce', 0.0), e[2], e[3]]
    elif len(e) == 5:
        ala3_omm_energies = [e[0], e[1], ('PeriodicTorsionForce', 0.0), e[2], e[3], e[4]]

    e = villin_omm_energies
    if len(e) == 4:
        villin_omm_energies = [e[0], e[1], ('PeriodicTorsionForce', 0.0), ('PeriodicTorsionForce', 0.0), e[2], e[3]]
    elif len(e) == 5:
        villin_omm_energies = [e[0], e[1], ('PeriodicTorsionForce', 0.0), e[2], e[3], e[4]]

    # dev if you don't want to assert
    return (ala3_amber_energies, ala3_omm_energies, villin_amber_energies,
           villin_omm_energies)

    print('Asserting ala_ala_ala energies...')
    counter = 0
    for i, j in zip(ala3_amber_energies, ala3_omm_energies):
        if counter != 3: # Not Impropers
            testing.assert_allclose(j[1], i[1], rtol=1e-5,
            err_msg=('Ala_ala_ala energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
        else: # Impropers - higher tolerance
            try:
                testing.assert_allclose(j[1], i[1], rtol=1e-2)
            except AssertionError:
                testing.assert_allclose(j[1], i[1], rtol=2e-2,
                err_msg=('Ala_ala_ala energies outside of allowed tolerance for %s' % ffxml_name))
                warnings.warn('Ala_ala_ala impropers failed assertion at 1% tolerance, but '
                              'have been asserted at the higher 2% tolerance')
            finally:
                counter += 1
    print('Ala_ala_ala energy validation successful!')

    print('Asserting villin headpiece energies...')
    counter = 0
    for i, j in zip(villin_amber_energies, villin_omm_energies):
        if counter != 3: # Not Impropers
            testing.assert_allclose(j[1], i[1], rtol=1e-5,
            err_msg=('Villin energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
        else: # Impropers - higher tolerance
            try:
                testing.assert_allclose(j[1], i[1], rtol=1e-2)
            except AssertionError:
                testing.assert_allclose(j[1], i[1], rtol=2e-2,
                err_msg=('Villin energies outside of allowed tolerance for %s' % ffxml_name))
                warnings.warn('Villin impropers failed assertion at 1% tolerance, but '
                              'have been asserted at the higher 2% tolerance')
            finally:
                counter += 1
    print('Villin headpiece energy validation successful!')
    # dev
    #return (ala3_amber_energies, ala3_omm_energies, villin_amber_energies,
    #       villin_omm_energies)
    print('Done!')

def validate_nucleic(ffxml_name, leaprc_name):
    if not find_executable('tleap'):
        raise RuntimeError('tleap not available from PATH')
    print('Preparing temporary files for validation...')
    dna_top = tempfile.mkstemp()
    dna_crd = tempfile.mkstemp()
    leap_script_dna_file = tempfile.mkstemp()
    rna_top = tempfile.mkstemp()
    rna_crd = tempfile.mkstemp()
    leap_script_rna_file = tempfile.mkstemp()
    leap_script_rna_file_alt = tempfile.mkstemp()

    print('Preparing LeaP scripts...')
    leap_script_dna_string = """addPdbAtomMap {
{ "H1'" "H1*" }
{ "H2'" "H2'1" }
{ "H2''" "H2'2" }
{ "H3'" "H3*" }
{ "H4'" "H4*" }
{ "H5'" "H5'1" }
{ "H5''" "H5'2" }
{ "HO2'" "HO'2" }
{ "HO5'" "H5T"  }
{ "HO3'" "H3T" }
{ "OP1" "O1P" }
{ "OP2" "O2P" }
}
source %s
addPdbResMap {
{ 0 "DG" "DG5"  } { 1 "DG" "DG3"  }
{ 0 "DA" "DA5"  } { 1 "DA" "DA3"  }
{ 0 "DC" "DC5"  } { 1 "DC" "DC3"  }
{ 0 "DT" "DT5"  } { 1 "DT" "DT3"  }
}
x = loadPdb files/4rzn_dna.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, dna_top[1], dna_crd[1])

    leap_script_rna_string = """
addPdbAtomMap {
{ "H1'" "H1*" }
{ "H2'" "H2'1" }
{ "H2''" "H2'2" }
{ "H3'" "H3*" }
{ "H4'" "H4*" }
{ "H5'" "H5'1" }
{ "H5''" "H5'2" }
{ "HO2'" "HO'2" }
{ "HO5'" "H5T"  }
{ "HO3'" "H3T" }
{ "OP1" "O1P" }
{ "OP2" "O2P" }
}
source %s
addPdbResMap {
{ 0 "G" "G5"  } { 1 "G" "G3"  } { "G" "G" }
{ 0 "A" "A5"  } { 1 "A" "A3"  } { "A" "A" }
{ 0 "C" "C5"  } { 1 "C" "C3"  } { "C" "C" }
{ 0 "U" "U5"  } { 1 "U" "U3"  } { "U" "U" }
}
x = loadPdb files/5c5w_rna.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, rna_top[1], rna_crd[1])

    leap_script_rna_string_alt = """
addPdbAtomMap {
{ "H1'" "H1*" }
{ "H2'" "H2'1" }
{ "H2''" "H2'2" }
{ "H3'" "H3*" }
{ "H4'" "H4*" }
{ "H5'" "H5'1" }
{ "H5''" "H5'2" }
{ "HO2'" "HO'2" }
{ "HO5'" "H5T"  }
{ "HO3'" "H3T" }
{ "OP1" "O1P" }
{ "OP2" "O2P" }
}
source %s
addPdbResMap {
{ 0 "G" "RG5"  } { 1 "G" "RG3"  } { "G" "RG" }
{ 0 "A" "RA5"  } { 1 "A" "RA3"  } { "A" "RA" }
{ 0 "C" "RC5"  } { 1 "C" "RC3"  } { "C" "RC" }
{ 0 "U" "RU5"  } { 1 "U" "RU3"  } { "U" "RU" }
}
x = loadPdb files/5c5w_rna.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, rna_top[1], rna_crd[1])

    os.write(leap_script_dna_file[0], leap_script_dna_string)
    os.write(leap_script_rna_file[0], leap_script_rna_string)
    os.write(leap_script_rna_file_alt[0], leap_script_rna_string_alt)

    print('Running LEaP...')
    os.system('tleap -f %s' % leap_script_dna_file[1])
    if os.path.getsize(dna_top[1]) == 0 or os.path.getsize(dna_crd[1]) == 0:
        raise RuntimeError('DNA LEaP fail for %s' % leaprc_name)
    os.system('tleap -f %s' % leap_script_rna_file[1])
    if os.path.getsize(rna_top[1]) == 0 or os.path.getsize(rna_crd[1]) == 0:
        # try alternative name mappings
        os.system('tleap -f %s' % leap_script_rna_file_alt[1])
    if os.path.getsize(rna_top[1]) == 0 or os.path.getsize(rna_crd[1]) == 0:
        raise RuntimeError('RNA LEaP fail for %s' % leaprc_name)

    print('Calculating DNA energies...')
    # AMBER
    parm_amber = parmed.load_file(dna_top[1], dna_crd[1])
    system_amber = parm_amber.createSystem(splitDihedrals=True)
    dna_amber_energies = parmed.openmm.energy_decomposition_system(parm_amber, system_amber, nrg=kilojoules_per_mole)
    # OpenMM
    ff = app.ForceField(ffxml_name)
    system_omm = ff.createSystem(parm_amber.topology)
    parm_omm = parmed.openmm.load_topology(parm_amber.topology, system_omm, xyz=parm_amber.positions)
    system_omm = parm_omm.createSystem(splitDihedrals=True)
    dna_omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=kilojoules_per_mole)

    print('Calculating RNA energies...')
    # AMBER
    parm_amber = parmed.load_file(rna_top[1], rna_crd[1])
    system_amber = parm_amber.createSystem(splitDihedrals=True)
    rna_amber_energies = parmed.openmm.energy_decomposition_system(parm_amber, system_amber, nrg=kilojoules_per_mole)
    # OpenMM
    ff = app.ForceField(ffxml_name)
    system_omm = ff.createSystem(parm_amber.topology)
    parm_omm = parmed.openmm.load_topology(parm_amber.topology, system_omm, xyz=parm_amber.positions)
    system_omm = parm_omm.createSystem(splitDihedrals=True)
    rna_omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=kilojoules_per_mole)

    print('Deleting temp files...')
    for f in (dna_top, dna_crd, leap_script_dna_file, rna_top, rna_crd,
              leap_script_rna_file, leap_script_rna_file_alt):
        os.close(f[0])
        os.unlink(f[1])

    # temp for testing
    #return (dna_amber_energies, dna_omm_energies, rna_amber_energies,
    #        rna_omm_energies)
    # temp for debugging fixes
    e = dna_omm_energies
    if len(e) == 4:
        dna_omm_energies = [e[0], e[1], ('PeriodicTorsionForce', 0.0), ('PeriodicTorsionForce', 0.0), e[2], e[3]]
    elif len(e) == 5:
        dna_omm_energies = [e[0], e[1], ('PeriodicTorsionForce', 0.0), e[2], e[3], e[4]]

    e = rna_omm_energies
    if len(e) == 4:
        rna_omm_energies = [e[0], e[1], ('PeriodicTorsionForce', 0.0), ('PeriodicTorsionForce', 0.0), e[2], e[3]]
    elif len(e) == 5:
        rna_omm_energies = [e[0], e[1], ('PeriodicTorsionForce', 0.0), e[2], e[3], e[4]]

    # dev
    return (dna_amber_energies, dna_omm_energies, rna_amber_energies,
           rna_omm_energies)

    print('Asserting DNA energies...')
    counter = 0
    for i, j in zip(dna_amber_energies, dna_omm_energies):
        if counter != 3: # Not Impropers
            testing.assert_allclose(j[1], i[1], rtol=1e-5,
            err_msg=('DNA energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
        else: # Impropers - higher tolerance
            testing.assert_allclose(j[1], i[1], rtol=1e-2,
            err_msg=('DNA energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
    print('DNA energy validation successful!')

    print('Asserting RNA energies...')
    counter = 0
    for i, j in zip(rna_amber_energies, rna_omm_energies):
        if counter != 3: # Not Impropers
            testing.assert_allclose(j[1], i[1], rtol=1e-5,
            err_msg=('RNA energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
        else: # Impropers - higher tolerance
            testing.assert_allclose(j[1], i[1], rtol=1e-2,
            err_msg=('RNA energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
    print('RNA energy validation successful!')
    # dev
    #return (dna_amber_energies, dna_omm_energies, rna_amber_energies,
    #       rna_omm_energies)

def validate_gaff(ffxml_name, leaprc_name):
    if not find_executable('tleap'):
        raise RuntimeError('tleap not available from PATH')
    print('Preparing temporary files for validation...')
    imatinib_top = tempfile.mkstemp()
    imatinib_crd = tempfile.mkstemp()
    leap_script_imatinib_file = tempfile.mkstemp()

    print('Preparing LeaP scripts...')
    leap_script_imatinib_string = """source %s
x = loadMol2 files/imatinib.mol2
saveAmberParm x %s %s
quit""" % (leaprc_name, imatinib_top[1], imatinib_crd[1])
    os.write(leap_script_imatinib_file[0], leap_script_imatinib_string)

    print('Running LEaP...')
    os.system('tleap -f %s' % leap_script_imatinib_file[1])
    if os.path.getsize(imatinib_top[1]) == 0 or os.path.getsize(imatinib_crd[1]) == 0:
        raise RuntimeError('imatinib LEaP fail for %s' % leaprc_name)

    print('Calculating imatinib energies...')
    # AMBER
    parm_amber = parmed.load_file(imatinib_top[1], imatinib_crd[1])
    system_amber = parm_amber.createSystem()
    imatinib_amber_energies = parmed.openmm.energy_decomposition_system(parm_amber, system_amber, nrg=kilojoules_per_mole)
    # OpenMM
    ff = app.ForceField(ffxml_name, 'files/imatinib_frcmod.xml', 'files/imatinib.xml')
    system_omm = ff.createSystem(parm_amber.topology)
    parm_omm = parmed.openmm.load_topology(parm_amber.topology, system_omm, xyz=parm_amber.positions)
    system_omm = parm_omm.createSystem()
    imatinib_omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=kilojoules_per_mole)

    print('Deleting temp files...')
    for f in (imatinib_top, imatinib_crd, leap_script_imatinib_file):
        os.close(f[0])
        os.unlink(f[1])

    print('Asserting imatinib energies...')
    for i, j in zip(imatinib_amber_energies, imatinib_omm_energies):
        testing.assert_allclose(j[1], i[1], rtol=1e-5,
        err_msg=('Imatinib energies outside of allowed tolerance for %s' % ffxml_name))
    print('Imatinib energy validation successful!')
    print('Done!')

def diagnostics_modify_leaprc():
    if not find_executable('tleap'):
        raise RuntimeError('tleap not available from PATH')
    if os.getenv('AMBERHOME'):
        AMBERHOME = os.getenv('AMBERHOME')
    else:
        tleap_path = find_executable('tleap')
        AMBERHOME = os.path.split(tleap_path)[0]
        AMBERHOME = os.path.join(AMBERHOME, '../')
        parmed.amber.AMBERHOME = AMBERHOME
    import copy
    from parmed.utils.six import iteritems
    import shutil
    # open log
    log = open('conv.log', 'w')
    log.write('Start log\n')
    ignore = {'solvents.lib', 'atomic_ions.lib', 'ions94.lib', 'ions91.lib',
          'phosphoaa10.lib'}
    leaprc_path = os.path.join(AMBERHOME, 'dat/leap/cmd/leaprc.ff14SB')
    # pre-process leaprc through the ignore mechanism
    with open(leaprc_path) as q:
        lines = map(lambda line:
                line if '#' not in line else line[:line.index('#')], q)
        new_lines = []
        for line in lines:
            if _loadoffre.findall(line) and _loadoffre.findall(line)[0] in ignore:
                continue
            else:
                new_lines.append(line)
    leaprc_new = StringIO(''.join(new_lines))

    parm = parmed.amber.AmberParameterSet.from_leaprc(leaprc_new)
    # turn off dihedrals
    for dih in parm.dihedral_types:
        for typ in parm.dihedral_types[dih]:
            typ.phi_k = 0 # set all k's to zero to turn off

    # remove impropers with unused types
    keep_types = set()
    for name, residue in iteritems(parm.residues):
        for atom in residue.atoms:
            keep_types.add(atom.type)
    skip_types = {typ for typ in parm.atom_types if typ not in keep_types}
    improper = OrderedDict()
    for types, imp in iteritems(parm.improper_periodic_types):
        if not (types[0] in skip_types or types[1] in skip_types or
           types[2] in skip_types or types[3] in skip_types):
           improper[types] = imp
    for typ in parm.atom_types:
        if typ in skip_types:
            del parm.atom_types[typ]
    # prep sums of energies
    a_sum=b_sum=c_sum=d_sum=e_sum=f_sum=g_sum=h_sum = 0
    # now turn off all impropers BUT one
    for imp in improper:
        log = open('conv.log', 'a')
        improper_new = copy.deepcopy(improper)
        leaprc_mod_file = tempfile.mkstemp()
        dat_mod_file = tempfile.mkstemp()
        for imp2 in improper_new:
            if imp == imp2: continue
            improper_new[imp2].phi_k = 0
        parm.improper_periodic_types = improper_new
        parm.write(dat_mod_file[1], style='parm')
        leaprc_mod_string = """
logFile leap.log
#
# ----- leaprc for loading the ff14SB force field
# ----- NOTE: this is designed for PDB format 3!
#    Uses frcmod.ff14SB for proteins; ff99bsc0 for DNA; ff99bsc0_chiOL3 for RNA
#
#	load atom type hybridizations
#
addAtomTypes {
{ "H"   "H" "sp3" }
{ "HO"  "H" "sp3" }
{ "HS"  "H" "sp3" }
{ "H1"  "H" "sp3" }
{ "H2"  "H" "sp3" }
{ "H3"  "H" "sp3" }
{ "H4"  "H" "sp3" }
{ "H5"  "H" "sp3" }
{ "HW"  "H" "sp3" }
{ "HC"  "H" "sp3" }
{ "HA"  "H" "sp3" }
{ "HP"  "H" "sp3" }
{ "HZ"  "H" "sp3" }
{ "OH"  "O" "sp3" }
{ "OS"  "O" "sp3" }
{ "O"   "O" "sp2" }
{ "O2"  "O" "sp2" }
{ "OP"  "O" "sp2" }
{ "OW"  "O" "sp3" }
{ "CT"  "C" "sp3" }
{ "CX"  "C" "sp3" }
{ "C8"  "C" "sp3" }
{ "2C"  "C" "sp3" }
{ "3C"  "C" "sp3" }
{ "CH"  "C" "sp3" }
{ "CS"  "C" "sp2" }
{ "C"   "C" "sp2" }
{ "CO"   "C" "sp2" }
{ "C*"  "C" "sp2" }
{ "CA"  "C" "sp2" }
{ "CB"  "C" "sp2" }
{ "CC"  "C" "sp2" }
{ "CN"  "C" "sp2" }
{ "CM"  "C" "sp2" }
{ "CK"  "C" "sp2" }
{ "CQ"  "C" "sp2" }
{ "CD"  "C" "sp2" }
{ "C5"  "C" "sp2" }
{ "C4"  "C" "sp2" }
{ "CP"  "C" "sp2" }
{ "CI"  "C" "sp3" }
{ "CJ"  "C" "sp2" }
{ "CW"  "C" "sp2" }
{ "CV"  "C" "sp2" }
{ "CR"  "C" "sp2" }
{ "CA"  "C" "sp2" }
{ "CY"  "C" "sp2" }
{ "C0"  "Ca" "sp3" }
{ "MG"  "Mg" "sp3" }
{ "N"   "N" "sp2" }
{ "NA"  "N" "sp2" }
{ "N2"  "N" "sp2" }
{ "N*"  "N" "sp2" }
{ "NP"  "N" "sp2" }
{ "NQ"  "N" "sp2" }
{ "NB"  "N" "sp2" }
{ "NC"  "N" "sp2" }
{ "NT"  "N" "sp3" }
{ "NY"  "N" "sp2" }
{ "N3"  "N" "sp3" }
{ "S"   "S" "sp3" }
{ "SH"  "S" "sp3" }
{ "P"   "P" "sp3" }
{ "LP"  ""  "sp3" }
{ "EP"  ""  "sp3" }
{ "F"   "F" "sp3" }
{ "Cl"  "Cl" "sp3" }
{ "Br"  "Br" "sp3" }
{ "I"   "I"  "sp3" }
{ "F-"   "F" "sp3" }
{ "Cl-"  "Cl" "sp3" }
{ "Br-"  "Br" "sp3" }
{ "I-"   "I"  "sp3" }
{ "Li+"  "Li"  "sp3" }
{ "Na+"  "Na"  "sp3" }
{ "K+"  "K"  "sp3" }
{ "Rb+"  "Rb"  "sp3" }
{ "Cs+"  "Cs"  "sp3" }
{ "Mg+"  "Mg"  "sp3" }
# glycam
{ "OG"  "O" "sp3" }
{ "OL"  "O" "sp3" }
{ "AC"  "C" "sp3" }
{ "EC"  "C" "sp3" }
}
#
#	Load the main parameter set.
#
parm10 = loadamberparams %s
#frcmod14SB = loadamberparams frcmod.ff14SB
#
#	Load main chain and terminating amino acid libraries, nucleic acids
#
loadOff amino12.lib
loadOff aminoct12.lib
loadOff aminont12.lib
loadOff nucleic12.lib
#
#       Load water and ions
#
loadOff atomic_ions.lib
loadOff solvents.lib
HOH = TP3
WAT = TP3

#
#	Define the PDB name map for the amino acids and nucleic acids
#
addPdbResMap {
{ 0 "HYP" "NHYP" } { 1 "HYP" "CHYP" }
{ 0 "ALA" "NALA" } { 1 "ALA" "CALA" }
{ 0 "ARG" "NARG" } { 1 "ARG" "CARG" }
{ 0 "ASN" "NASN" } { 1 "ASN" "CASN" }
{ 0 "ASP" "NASP" } { 1 "ASP" "CASP" }
{ 0 "CYS" "NCYS" } { 1 "CYS" "CCYS" }
{ 0 "CYX" "NCYX" } { 1 "CYX" "CCYX" }
{ 0 "GLN" "NGLN" } { 1 "GLN" "CGLN" }
{ 0 "GLU" "NGLU" } { 1 "GLU" "CGLU" }
{ 0 "GLY" "NGLY" } { 1 "GLY" "CGLY" }
{ 0 "HID" "NHID" } { 1 "HID" "CHID" }
{ 0 "HIE" "NHIE" } { 1 "HIE" "CHIE" }
{ 0 "HIP" "NHIP" } { 1 "HIP" "CHIP" }
{ 0 "ILE" "NILE" } { 1 "ILE" "CILE" }
{ 0 "LEU" "NLEU" } { 1 "LEU" "CLEU" }
{ 0 "LYS" "NLYS" } { 1 "LYS" "CLYS" }
{ 0 "MET" "NMET" } { 1 "MET" "CMET" }
{ 0 "PHE" "NPHE" } { 1 "PHE" "CPHE" }
{ 0 "PRO" "NPRO" } { 1 "PRO" "CPRO" }
{ 0 "SER" "NSER" } { 1 "SER" "CSER" }
{ 0 "THR" "NTHR" } { 1 "THR" "CTHR" }
{ 0 "TRP" "NTRP" } { 1 "TRP" "CTRP" }
{ 0 "TYR" "NTYR" } { 1 "TYR" "CTYR" }
{ 0 "VAL" "NVAL" } { 1 "VAL" "CVAL" }
{ 0 "HIS" "NHIS" } { 1 "HIS" "CHIS" }
{ 0 "G" "G5"  } { 1 "G" "G3"  }
{ 0 "A" "A5"  } { 1 "A" "A3"  }
{ 0 "C" "C5"  } { 1 "C" "C3"  }
{ 0 "U" "U5"  } { 1 "U" "U3"  }
{ 0 "DG" "DG5"  } { 1 "DG" "DG3"  }
{ 0 "DA" "DA5"  } { 1 "DA" "DA3"  }
{ 0 "DC" "DC5"  } { 1 "DC" "DC3"  }
{ 0 "DT" "DT5"  } { 1 "DT" "DT3"  }
#  some old Amber residue names for RNA:
{ 0  "RA5" "A5" } { 1 "RA3" "A3"} {"RA" "A" }
{ 0  "RC5" "C5" } { 1 "RC3" "C3"} {"RC" "C" }
{ 0  "RG5" "G5" } { 1 "RG3" "G3"} {"RG" "G" }
{ 0  "RU5" "U5" } { 1 "RU3" "U3"} {"RU" "U" }
#  some really old Amber residue names, assuming DNA:
{ 0 "GUA" "DG5"  } { 1 "GUA" "DG3"  } { "GUA" "DG" }
{ 0 "ADE" "DA5"  } { 1 "ADE" "DA3"  } { "ADE" "DA" }
{ 0 "CYT" "DC5"  } { 1 "CYT" "DC3"  } { "CYT" "DC" }
{ 0 "THY" "DT5"  } { 1 "THY" "DT3"  } { "THY" "DT" }
#  uncomment out the following if you have this old style RNA files:
# { 0 "GUA" "G5"  } { 1 "GUA" "G3"  } { "GUA" "G" }
# { 0 "ADE" "A5"  } { 1 "ADE" "A3"  } { "ADE" "A" }
# { 0 "CYT" "C5"  } { 1 "CYT" "C3"  } { "CYT" "C" }
# { 0 "URA" "R5"  } { 1 "URA" "R3"  } { "URA" "R" }

}

#  try to be good about reading in really old atom names as well:
addPdbAtomMap {
{ "O5*" "O5'" }
{ "C5*" "C5'" }
{ "C4*" "C4'" }
{ "O4*" "O4'" }
{ "C3*" "C3'" }
{ "O3*" "O3'" }
{ "C2*" "C2'" }
{ "O2*" "O2'" }
{ "C1*" "C1'" }
{ "C5M" "C7"  }
{ "H1*" "H1'" }
{ "H2*1" "H2'" }
{ "H2*2" "H2''" }
{ "H2'1" "H2'" }
{ "H2'2" "H2''" }
{ "H3*" "H3'" }
{ "H4*" "H4'" }
{ "H5*1" "H5'" }
{ "H5*2" "H5''" }
{ "H5'1" "H5'" }
{ "H5'2" "H5''" }
{ "HO'2" "HO2'" }
{ "H5T"  "HO5'" }
{ "H3T"  "HO3'" }
{ "O1'" "O4'" }
{ "OA"  "OP1" }
{ "OB"  "OP2" }
{ "O1P" "OP1" }
{ "O2P" "OP2" }
}

#
# assume that most often proteins use HIE
#
NHIS = NHIE
HIS = HIE
CHIS = CHIE""" % dat_mod_file[1]

        os.write(leaprc_mod_file[0], leaprc_mod_string)

        #run stuff
        dat_closed = False
        try:
            ffxml_name = convert(leaprc_mod_file[1], ignore=ignore)
            #validate_protein(ffxml_name, leaprc_mod_file[1])
            #validate_nucleic(ffxml_name, leaprc_mod_file[1])
            (a, b, c, d) = validate_protein(ffxml_name, leaprc_mod_file[1])
            (e, f, g, h)= validate_nucleic(ffxml_name, leaprc_mod_file[1])
            # add to sums
            a_sum += a[3][1]
            b_sum += b[3][1]
            c_sum += c[3][1]
            d_sum += d[3][1]
            e_sum += e[3][1]
            f_sum += f[3][1]
            g_sum += g[3][1]
            h_sum += h[3][1]
            # dev - save dat file too
            dat_save_name = ffxml_name.split('.')[0] + '.dat'
            os.close(dat_mod_file[0])
            dat_closed = True
            shutil.copyfile(dat_mod_file[1], dat_save_name)
        finally:
            os.close(leaprc_mod_file[0])
            if not dat_closed:
                os.close(dat_mod_file[0])
            os.unlink(leaprc_mod_file[1])
            os.unlink(dat_mod_file[1])
        #(a, b, c, d) = validate_protein(ffxml_name, leaprc_mod_file[1])
        #(e, f, g, h)= validate_nucleic(ffxml_name, leaprc_mod_file[1])
        #for i in (a,b,c,d,e,f,g,h):
        #    print(i)
        log.write(str(imp))
        log.write('\n')
        log.write(ffxml_name)
        log.write('\n')
        log.write(str(a[3]))
        log.write('\n')
        log.write(str(b[3]))
        log.write('\n\n')
        log.write(str(c[3]))
        log.write('\n')
        log.write(str(d[3]))
        log.write('\n\n')
        log.write(str(e[3]))
        log.write('\n')
        log.write(str(f[3]))
        log.write('\n\n')
        log.write(str(g[3]))
        log.write('\n')
        log.write(str(h[3]))
        log.write('\n\n')
        #log.write(str(improper_new))
        #log.write('\n\n')
        log.close()

    #write sums
    log = open('conv.log', 'a')
    log.write('Totals\n')
    log.write(str(a_sum))
    log.write('\n')
    log.write(str(b_sum))
    log.write('\n\n')
    log.write(str(c_sum))
    log.write('\n')
    log.write(str(d_sum))
    log.write('\n\n')
    log.write(str(e_sum))
    log.write('\n')
    log.write(str(f_sum))
    log.write('\n\n')
    log.write(str(g_sum))
    log.write('\n')
    log.write(str(h_sum))
    log.write('\n\n')
    log.close()

if __name__ == '__main__':
    diagnostics_modify_leaprc()
