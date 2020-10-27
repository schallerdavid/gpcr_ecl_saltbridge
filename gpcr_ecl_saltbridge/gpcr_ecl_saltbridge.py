"""
gpcr_ecl_saltbridge.py
Identify GPCR structures with saltbridges between extracellular loops 2 and 3.

Handles all functions.
"""
# python standard library
from collections import OrderedDict
import json
import logging
import pathlib
import pickle
from urllib.error import HTTPError
from urllib.request import urlopen, urlretrieve

# external libraries
from Bio import pairwise2


DATA = pathlib.Path(__file__).parent.absolute() / 'data'


def distance(x, y):
    """ This function returns the euclidean distance between two point in three dimensional space. """
    return ((x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2 + (x[2] - y[2]) ** 2) ** 0.5


def get_gpcr_structures():
    url = 'https://gpcrdb.org/services/structure/'
    response = urlopen(url)
    structures = json.loads(response.read().decode('utf-8'))
    return structures


def get_gpcr_sequence_dict(protein_name):
    url = 'https://gpcrdb.org/services/residues/{}'.format(protein_name)
    response = urlopen(url)
    protein = json.loads(response.read().decode('utf-8'))
    sequence_dict = OrderedDict()
    for residue in protein:
        sequence_dict[residue['sequence_number']] = [residue['amino_acid'], residue['display_generic_number']]
    return sequence_dict


def generate_pdb_code_dict(directory=DATA):
    """ Retrieves data from the GPCRDB to build or update a dictionary with pdb code as key and a list
    with gene name and preferred chain as value. It can be saved as pickle-file in a data directory in the provided
    directory. """
    directory.mkdir(parents=True, exist_ok=True)
    pdb_code_dict_path = directory / 'pdb_code_dict.pkl'

    if pdb_code_dict_path.is_file():
        pdb_code_dict = pickle.load(open(pdb_code_dict_path, 'rb'))
    else:
        pdb_code_dict = {}

    structures = get_gpcr_structures()
    for structure in structures:
        if structure['pdb_code'] not in pdb_code_dict.keys():
            pdb_code_dict[structure['pdb_code']] = [structure['protein'], structure['preferred_chain'].split(',')[0]]

    with open(pdb_code_dict_path, 'wb') as pickle_file:
        pickle.dump(pdb_code_dict, pickle_file, pickle.HIGHEST_PROTOCOL)

    return pdb_code_dict


def generate_sequence_dict(directory=DATA):
    """ This functions retrieves data from the GPCRDB to build or update a dictionary with gene name as key and a
    dictionary as value with residue id as key and a list with amino acid and GPCRDB number as value. It can be saved
    as pickle-file in a data directory in the provided directory. """
    directory.mkdir(parents=True, exist_ok=True)
    sequence_dict_path = directory / 'sequence_dict.pkl'

    if sequence_dict_path.is_file():
        sequence_dict = pickle.load(open(sequence_dict_path, 'rb'))
    else:
        sequence_dict = {}

    structures = get_gpcr_structures()
    for structure in structures:
        protein_name = structure['protein']
        if protein_name not in sequence_dict.keys():
            sequence_dict[protein_name] = get_gpcr_sequence_dict(protein_name)

    with open(sequence_dict_path, 'wb') as pickle_file:
        pickle.dump(sequence_dict, pickle_file, pickle.HIGHEST_PROTOCOL)

    return sequence_dict


def download_pdb_files(pdb_codes, directory=DATA):
    """ Download pdb files from the PDB specified in the provided pdb code list and saves the files in
    the provided directory. """
    directory.mkdir(parents=True, exist_ok=True)

    for pdb_code in pdb_codes:
        try:
            file_path = DATA / f'{pdb_code}.pdb'
            if not file_path.is_file():
                urlretrieve(f'https://files.rcsb.org/download/{pdb_code}.pdb', file_path)
                logging.debug(f"Downloaded pdb file for {pdb_code} ...")
        except HTTPError:
            file_path = DATA / f'{pdb_code}.cif'
            if not file_path.is_file():
                urlretrieve(f'https://files.rcsb.org/download/{pdb_code}.cif', file_path)
                logging.debug(f"Downloaded mmcif file for {pdb_code} ...")

    return


def update_data(directory=DATA):
    """ This function updates data retrieved from GPCRDB and PDB. """

    pdb_code_dict = generate_pdb_code_dict(directory=directory)
    sequences_dict = generate_sequence_dict(directory=directory)
    download_pdb_files(pdb_code_dict.keys(), directory=directory)

    return pdb_code_dict, sequences_dict


def read_pdb_structure(pdb_code, directory=DATA):
    from Bio.PDB import PDBParser, MMCIFParser

    file_path = directory / f'{pdb_code}.pdb'
    if file_path.is_file():
        parser = PDBParser()
    else:
        parser = MMCIFParser()
        file_path = directory / f'{pdb_code}.cif'
    structure = parser.get_structure(pdb_code, file_path)
    return structure


def generate_pdb_sequence_dict(pdb_code, preferred_chain, directory=DATA):

    from Bio.SeqUtils import seq1

    structure = read_pdb_structure(pdb_code, directory=directory)
    pdb_sequence_dict = OrderedDict()

    for residue in structure[0][preferred_chain].get_residues():
        if residue.get_id()[0] == ' ':  # no hetero atoms
            pdb_sequence_dict[residue.get_id()[1]] = seq1(residue.get_resname())

    return pdb_sequence_dict


def assign_generic_numbers_to_pdb(sequence_dict, pdb_sequence_dict):

    sequence = ''.join([sequence_dict[residue][0] for residue in sequence_dict.keys()])
    pdb_sequence = ''.join([pdb_sequence_dict[residue] for residue in pdb_sequence_dict.keys()])
    alignment = pairwise2.align.globalxs(sequence, pdb_sequence, -10, 0)[0]
    pdb_generic_numbers_dict = {}
    pdb_sequence_ids = list(pdb_sequence_dict.keys())
    counter = 1
    pdb_counter = 0
    for residue, pdb_residue in zip(alignment[0], alignment[1]):
        if residue != '-' and pdb_residue != '-':
            pdb_generic_numbers_dict[pdb_sequence_ids[pdb_counter]] = sequence_dict[counter][1]
            counter += 1
            pdb_counter += 1
        else:
            if residue != '-':
                counter += 1
            if pdb_residue != '-':
                pdb_counter += 1
    return pdb_generic_numbers_dict


def salt_bridges(directory=DATA):
    """ This function analyzes pdb files to contain salt bridges between ECL2 and ECL3 and returns the pdb codes. """
    ni_residues = {'ASP': ['OD1', 'OD2'], 'GLU': ['OE1', 'OE2']}
    pi_residues = {'ARG': ['NE', 'NH1', 'NH2'], 'HIS': ['ND1', 'NE2'], 'LYS': ['NZ']}
    pdb_code_dict, sequences_dict = update_data(directory=directory)
    salt_bridge_dict = {}
    for pdb_code in pdb_code_dict.keys():
        print(f'Analyzing {pdb_code} ...')
        distances = []
        ecl2_ni, ecl2_pi, ecl3_ni, ecl3_pi = [], [], [], []
        protein_name, preferred_chain = pdb_code_dict[pdb_code]
        try:
            pdb_sequence_dict = generate_pdb_sequence_dict(pdb_code, preferred_chain)
        except KeyError:
            print(f'Error for {pdb_code} ...')
            continue
        pdb_generic_numbers_dict = assign_generic_numbers_to_pdb(sequences_dict[protein_name], pdb_sequence_dict)
        structure = read_pdb_structure(pdb_code, directory=directory)
        h4, h5, h6, h7 = False, False, False, False
        h7_counter = 0
        for residue in structure[0][preferred_chain].get_residues():
            resid = residue.get_id()[1]
            if resid in pdb_generic_numbers_dict.keys():
                generic_number = pdb_generic_numbers_dict[resid]
                if generic_number is not None:
                    if generic_number.split('.')[0] == '4':
                        h4 = True
                    elif generic_number.split('.')[0] == '5':
                        h5 = True
                    elif generic_number.split('.')[0] == '6':
                        h6 = True
                    elif generic_number.split('.')[0] == '7':
                        h7 = True
                else:
                    generic_number = 'x.x'
                residue_name = residue.get_resname()
                for atom in residue.get_atoms():
                    atom_name = atom.get_name()
                    position = [x for x in atom.get_vector()]
                    if h4 and not h5:
                        if generic_number.split('.')[0] != '4':
                            if residue_name in ni_residues.keys():
                                if atom_name in ni_residues[residue_name]:
                                    ecl2_ni.append(position)
                            if residue_name in pi_residues.keys():
                                if atom_name in pi_residues[residue_name]:
                                    ecl2_pi.append(position)
                    if h6 and not h7:
                        if generic_number.split('.')[0] != '6':
                            if residue_name in ni_residues.keys():
                                if atom_name in ni_residues[residue_name]:
                                    ecl3_ni.append(position)
                            if residue_name in pi_residues.keys():
                                if atom_name in pi_residues[residue_name]:
                                    ecl3_pi.append(position)
                    if h7:
                        if h7_counter <= 3:
                            if atom_name == 'N':
                                h7_counter += 1
                            if residue_name in ni_residues.keys():
                                if atom_name in ni_residues[residue_name]:
                                    ecl3_ni.append(position)
                            if residue_name in pi_residues.keys():
                                if atom_name in pi_residues[residue_name]:
                                    ecl3_pi.append(position)
        for ni in ecl2_ni:
            for pi in ecl3_pi:
                distances.append(distance(ni, pi))
        for ni in ecl3_ni:
            for pi in ecl2_pi:
                distances.append(distance(ni, pi))
        if len(distances) > 0:
            if min(distances) < 5:
                print('Found salt bridge!')
                if protein_name in salt_bridge_dict.keys():
                    salt_bridge_dict[protein_name].append(pdb_code)
                else:
                    salt_bridge_dict[protein_name] = [pdb_code]
    for protein_name in salt_bridge_dict.keys():
        print(protein_name, ':', *salt_bridge_dict[protein_name])
    return salt_bridge_dict
