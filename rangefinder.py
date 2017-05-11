#! /usr/bin/python3
# -*- coding: utf-8 -*-
"""Predict dynamic range of sensors from crystal structures.

Josh Mitchell 2016"""
from __future__ import print_function

### Startup message ###

print("\
                             RANGEFINDER \n\
\n \
    Please cite:\n\
  Mitchell, J. A.; Whitfield, J. H.; Zhang, W. H.; O'Mara, M. L.; Jackson,\n\
  C. J. 2016. Rangefinder: A semi-synthetic FRET sensor design algorithm. \n\
  Under revision.\n\
")

import sys

# Crash if we're not on Python 3
try:
    if sys.version_info[0] < 3:
        sys.exit()
except SystemExit:
    file = sys.stdout
    python_version_number = '%i.%i.%i' % (sys.version_info[0], sys.version_info[1], sys.version_info[2])
    file.write('Rangefinder requires Python 3+, but was run with version %s. Exiting.\n' % python_version_number)
    sys.exit(2)

# Try to import Biopython, and give a nice error if it fails
try:
    import Bio
    from Bio import SeqIO, PDB, AlignIO
except ImportError:
    file = sys.stdout
    file.write('Rangefinder requires Biopython (http://biopython.org/wiki/Biopython), which can not be found. Exiting.\n')
    sys.exit(2)

# Import everything else
import argparse
import os.path
from math import sqrt
from random import shuffle,sample
from io import StringIO
import copy as _copy




### PARSE ARGUMENTS ###


# First some classes needed for argument parsing
def parse_range(astr):
    """Parser for comma-seperated series of ranges."""
    if astr=='all':
        return 'all'
    else:
        result = set()
        for part in astr.split(','):
            x = part.split('-')
            result.update(range(int(x[0]), int(x[-1]) + 1))
        return sorted(result)

class DyesArgumentError(Exception):
    """Error class for AppendDyesAction"""
    pass

class AppendDyesAction(argparse.Action):
    """argparse action to allow several dyes to be defined simply 

    Takes a single argument with one mandatory and 2 optional settings for each dye

    Adapted from argparse.py"""
    def _ensure_value(self, namespace, name, value):
        """Return namespace.name, or value if namespace.name is None """
        if getattr(namespace, name, None) is None:
            setattr(namespace, name, value)
        return getattr(namespace, name)

    def __init__(self,
                 option_strings,
                 dest,
                 nargs='+',
                 const=None,
                 default=None,
                 type=None,
                 choices=None,
                 required=False,
                 help=None,
                 metavar=None):
        if nargs != '+':
            raise ValueError('nargs must be * for AppendDyesAction to work as intended')
        super().__init__(
            option_strings=option_strings,
            dest=dest,
            nargs=nargs,
            const=const,
            default=default,
            type=type,
            choices=choices,
            required=required,
            help=help,
            metavar=metavar)


    def __call__(self, parser, namespace, values, option_string=None):
        items = _copy.copy(self._ensure_value(namespace, self.dest, []))
        # How many arguments have been passed to this instance?
        if len(values) == 1:
            # if it's just one, make up a title and fudge factor
            title = "R0={}/S={}".format(values[0], default_fudge_factor)
            out_tuple = (title, float(values[0]), default_fudge_factor) 
        elif len(values) == 2:
            # if it's two, use the default fudge factor
            out_tuple = (values[1], float(values[0]), default_fudge_factor) 
        elif len(values) == 3:
            # if it's three, use them all
            out_tuple = (values[1], float(values[0]), float(values[2])) 
        else:
            raise DyesArgumentError("argument " + option_string + ": expected one to three arguments")
        print("New dye \"{}\" defined with Förster distance {} and fudge factor {}.".format(*out_tuple))
        items.append(out_tuple)
        setattr(namespace, self.dest, items)

# Now any settings that are needed in the arguments
default_fudge_factor = 0.75

# Define the ArgumentParser
parser = argparse.ArgumentParser(description="Predict dynamic range of sensors from crystal structures", add_help = False)

required_named = parser.add_argument_group('required named arguments')
required_named.add_argument('-a', '--open', '--apo', type=str, nargs='+', required=True,
                            help="Structure(s) in the open/apo/ligand-unbound conformation.")
required_named.add_argument('-b', '--closed', '--holo', type=str, nargs='+', required=True,
                            help="Structure(s) in the closed/holo/ligand-bound conformation.")
required_named.add_argument('-f', '--fit', type=parse_range, required=True,
                            help="A comma-seperated list of ranges defining the residues by which to fit all structures.")

opt_args = parser.add_argument_group('optional arguments')
opt_args.add_argument('-h', '--help', action='help',
                      help='Show this help message and exit.')
opt_args.add_argument('-o', '--output', type=str, default='sensor_predict.csv',
                      help="The comma-seperated plain text file to output to.")
opt_args.add_argument('--fluorophore', type=int, default='0',
                      help="The residue number of the donor fluorophore - 0: use Rangefinder model (default).")
opt_args.add_argument('-r', '--residues', type=parse_range, default='all',
                      help="A comma-seperated list of ranges defining the residues for which dynamic ranges will be predicted.")
opt_args.add_argument('-n', '--name', type=str, default='CA', help="Atomname to use")
opt_args.add_argument('-l', '--linker', type=float, default=20.0, help="Distance from first residue CA to fluorophore (Angstroms)")
opt_args.add_argument('-d', '--dye', type=str, action=AppendDyesAction, metavar="PROPERTY",
                      help="Define a different pair of fluorophores for your experiment. Takes up to 3 properties: the Förster \
                      distance of the interaction, a name for the fluorophore pair, and a fudge factor describing cross-talk. \
                      Omitted properties are assumed. May be specified multiple times for multiple dyes.")
opt_args.add_argument('--advanced', action='store_true', help="Return predictions (rather than 0) for aligned residues.")
opt_args.add_argument('-v', '--verbose', action='store_true', help="Print additional output to terminal.")

# Parse arguments, and if the new AppendDyesAction throws an error, produce a pretty error to match argparse
try:
    args = parser.parse_args()
except DyesArgumentError:
    err = sys.exc_info()[1]
    file = sys.stdout
    file.write(parser.format_usage())
    args = {'prog': sys.argv[0], 'message': err}
    file.write('%(prog)s: error: %(message)s\n' % args)
    sys.exit(2)

# Set global variables used below to output from arguments
structures_to_load_open = args.open
structures_to_load_closed = args.closed

residues_to_check = args.residues
residues_to_align = args.fit

fp_resnum = str(args.fluorophore)

atom_name = args.name

linker_length = args.linker

advanced = args.advanced

quiet = not args.verbose

# Define dyes to predict
# List of (name, R_nought,fudge)
# If no dyes are defined in arguments, use the default (CFP-AF532)
if args.dye == None:
    global_dyes = []
    global_dyes.append(("CFP-AF532", 4.80, 0.82))
    # global_dyes.append(("CFP-OG488", 4.80, 0.70)) # More default dyes can be added like this
# Otherwise, convert the input dyes into tuples and stick 'em in a list
else:
    global_dyes = args.dye


calc_data = dict()

total_coords = 0
total_time = 0

def save_structure(structure, filename, quiet = False):
    """Write a PDB file from structure"""
    if not quiet:
        print_line = "Writing " + filename + "..."
        print(print_line, end="\r")
        longest_line_len = max(longest_line_len,len(print_line))
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(filename)
    if not quiet:
        final_print_str = "Wrote " + filename + "."
        num_spaces = max(0, longest_line_len - len(final_print_str))
        print(final_print_str + " "*num_spaces)

def centre_of_geometry(atoms):
    """ Return the centre of geometry of the input """
    coords = [a.get_coord() for a in atoms]
    x = mean([float(x) for x,y,z in coords])
    y = mean([float(y) for x,y,z in coords])
    z = mean([float(z) for x,y,z in coords])
    return (x, y, z)
    
def load_structures(files_to_load, quiet = False):
    """Load PDB files from a list and return a list of the structures"""
    parser = PDB.PDBParser(QUIET=True, PERMISSIVE=True)
    structures = []
    longest_line_len = 0
    for file in files_to_load:
        name = os.path.splitext(file)[0]
        if not quiet:
            print_line = "Loading " + name + "..."
            print(print_line, end="\r")
            longest_line_len = max(longest_line_len,len(print_line))
        new_structure = parser.get_structure(name, file)
        # Remove residue 0 to dedicate it to the donor fluorophore
        for new_model in new_structure:
            for new_chain in new_model:
                for residue in new_chain:
                    if residue.id[1] == 0:
                        new_chain.detach_child(residue.id)
        structures.append(new_structure)
        # save_structure(new_structure, name + ".no0.pdb")
    if not quiet:
        final_print_str = "Loaded " + str(len(files_to_load)) + " structures."
        num_spaces = max(0, longest_line_len - len(final_print_str))
        print(final_print_str + " "*num_spaces)
    return structures


def mean(data):
    """Return the sample arithmetic mean of data."""
    n = len(data)
    if n < 1:
        raise ValueError('mean requires at least one data point')
    return sum(data)/n # in Python 2 use sum(data)/float(n)

def select_subset(entity, structure = None, model = None, chain = None, residue = None, atom = None):
    """ Returns a list of atoms that are in the given entity, restricted by the parameters

    For example, if residue is given as 47, a list of atoms in residue 47 across all chains, models and structures will be returned"""
    # Associate each Biopython type with the appropriate parameter
    entity_restrictions = {
        Bio.PDB.Structure.Structure: structure, 
        Bio.PDB.Model.Model: model, 
        Bio.PDB.Chain.Chain: chain, 
        Bio.PDB.Residue.Residue: residue, 
        Bio.PDB.Atom.Atom: atom }

    # Associate each Biopython type with its child
    child_entity_type = {
        Bio.PDB.Structure.Structure: Bio.PDB.Model.Model, 
        Bio.PDB.Model.Model: Bio.PDB.Chain.Chain, 
        Bio.PDB.Chain.Chain: Bio.PDB.Residue.Residue, 
        Bio.PDB.Residue.Residue: Bio.PDB.Atom.Atom, 
        Bio.PDB.Atom.Atom: None }

    # Associate each Biopython type with the type its child must be selected by
    entity_selection_type = {
        Bio.PDB.Structure.Structure: int, 
        Bio.PDB.Model.Model: str, 
        Bio.PDB.Chain.Chain: int, 
        Bio.PDB.Residue.Residue: str, 
        Bio.PDB.Atom.Atom: None }

    # Get the type of the given entity's children
    type_of_entitys_children = child_entity_type[type(entity)]

    # Base case: if we're at the end of the entity chain (ie, atoms), return the current entity
    if type_of_entitys_children == None:
        return [entity]

    # Apply the appropriate restriction to the entity's children
    selected_children = []    
    restriction = entity_restrictions[type_of_entitys_children] 
    if restriction == None:
        selected_children = [child for child in entity]
    else:
        restriction_correct_type = entity_selection_type[type(entity)](restriction)
        try: # Try and select the desired child entity
            selected_children.append(entity[restriction_correct_type])
        except KeyError: # If the desired child isn't there, continue iterating over the parent entities
            pass

    subset = []
    for element in selected_children:
        subset += select_subset(element, structure, model, chain, residue, atom)
    return subset

def beads_from_model(residues, model, bead_name=atom_name):
    """Return a list of all beads with the specified name in the specified residue number range in model"""
    bead_dict = dict()
    for residue in model.get_residues():
        residue_seq_id = get_res_seq_id(residue)
        try:
            bead = residue[bead_name]
        except KeyError:
            continue
        if residue_seq_id in residues: 
            if residue_seq_id in bead_dict:
                raise Exception("Bead with residue sequence ID " + str(residue_seq_id) + " already in bead_dict: " + str(residue))
            bead_dict[residue_seq_id] = bead

    missing_residues = [key for key in residues if key not in bead_dict]
    if len(missing_residues) > 0:
        raise Exception("Residues " + str(missing_residues) + " are missing from model " + model)

    return bead_dict

def get_res_seq_id(residue):
    return residue.get_id()[1]


def align_structures(residues, structures, bead_name=atom_name, quiet = False):
    """Superimpose all models in structures by the desired beads in the nominated residues and return the target 

    The target is simply the first structure in the list"""
    # Take the first model in the first structure in structures, and pick out the 
    # BB beads in the appropriate residues so we have a target for superposition
    target_model = structures[0][0]

    target_beads_dict = beads_from_model(residues, target_model, bead_name)
    target_beads = list(target_beads_dict.values())

    longest_line_len = 0
    n=0
    conf = "open"
    for structure in structures:
        if not quiet:
            print_line = "Aligning " + str(structure) + "..."
            print(print_line, end="\r")
            longest_line_len = max(longest_line_len,len(print_line))
        for model in structure:
            # Get list of beads to align by
            mobile_beads_dict = beads_from_model(residues, model, bead_name)
            mobile_beads = list(mobile_beads_dict.values())
            # Superimpose the model on the target
            sup = PDB.Superimposer()
            sup.set_atoms(target_beads, mobile_beads) 
            sup.apply(model.get_atoms())
        n += 1
        if n == (len(structures)/2.0) + 1:
            conf = "closed"
            n=1
        structure_filename = conf + str(n) + ".pdb"
        # if not os.path.isfile(structure_filename):
        #     save_structure(structure, structure_filename)

    if not quiet:
        final_print_str = "Performed " + str(len(structures)) + " structure alignments."
        num_spaces = max(0, longest_line_len - len(final_print_str) + 12)
        print(final_print_str + " "*num_spaces)

    return target_model

def get_all_interfluorophore_distances(fp_resnum, 
                                       residues_to_check, 
                                       residues_to_align, 
                                       first_structures, 
                                       second_structures, 
                                       bead_name = atom_name, 
                                       quiet = False):
    """ Gets distances from all input structures based on alignment to all input structures in first_structures

    For consistency, I should align open to closed and closed to open and use fixed beads across everything
    So I should align everything to a single open structure, get all 60 open distances, then redo that with a closed structure. 
    Cool.
    That means this can be a single, simple function that I call twice with input lists concatenated in different orders.
    I will lose all the variation from movement in the sensing domain, but the whole point of this is that that's nonphysical
     OK. Sounds like a PLAN!
    To avoid bias based on the single structure I align to, I repeat the above with every structure as the alignment target. 
    Takes a bit longer but it works."""
    all_distances = {}
    longest_line_len = 0
    for structure in first_structures:
        # Make a list with all the structures except the current iteration
        other_first_structures = [x for x in first_structures if x != structure]
        # Now we can put the current iteration at the front
        all_structures = [structure] + other_first_structures + second_structures

        # Align all models in structures (Side effect)
        aligned_target_model = align_structures(residues_to_align, all_structures, bead_name=atom_name, quiet=quiet) 
 
        # Get all the fixed mutation candidate beads from the target structure
        fixed_mutation_beads = {}
        for mutation in residues_to_check:
            selected_beads = select_subset(aligned_target_model, residue = mutation, atom = bead_name)
            beads_found = len(selected_beads)
            if beads_found == 1:
                fixed_mutation_beads[mutation] = selected_beads[0]
            elif beads_found == 0:
                raise Exception("No bead found for mutation " + str(mutation))
            else:
                raise Exception(str(beads_found) + " beads found for mutation " + str(mutation))

        # Now grab the distances
        new_distances = {}
        for mutation,mutation_bead in fixed_mutation_beads.items():
            ### Get a list of pairs of fluorophore and mutation for each model
            beads_to_measure = []
            # Iterate over the chains in structures, remembering where we are at each step
            for structure in all_structures:
                # Throw up an error if we haven't been given a list of structures
                if type(structure) != Bio.PDB.Structure.Structure:
                    raise TypeError("all_structures is not an iterable of structures!")
                num_expected_beads = len(structure) # Should have 1 bead per model
                fluorophore_beads = select_subset(structure, residue = fp_resnum, atom = bead_name)
                # Make sure the numbers of everything match up, just to catch wierd errors
                if len(fluorophore_beads) == num_expected_beads:
                    beads_to_measure = beads_to_measure + [(fluor_bead, mutation_bead) for fluor_bead in fluorophore_beads]
                else:
                    raise Exception("Bead count mismatch in structure " + str(structure))
            
            # Now we have a list of pairs of beads between which distances must be measured, and we can just iterate over it!
            new_distances[mutation] = []
            for bead_pair in beads_to_measure:
                distance = bead_pair[0] - bead_pair[1]
                # distance = bead_pair["Fluorophore"] - bead_pair["Mutation"]
                new_distances[mutation].append(distance)
        #Stick them in the all_distances dict
        for mutation,distances in new_distances.items():
            try: 
                all_distances[mutation] +=  distances
            except KeyError:
                all_distances[mutation] =  distances

    return all_distances

def pick_single_FPpoint(atoms, distance = 20, bead_name=atom_name):
    """ Return the coordinates of a point distance units away from the first 
    atom in atoms with name bead_name along a line through the COG """
    cog = PDB.Vector(*centre_of_geometry(atoms))
    first_bead = PDB.Vector(0,0,0)
    # Pick the first bead_name atom in atoms and read coords
    for atom in atoms:
        if atom.get_name() == bead_name:
            first_bead = atom.get_vector()
            break
    # Find a unit vector pointing along the line from cog to first_bead
    unit_vector = first_bead - cog
    unit_vector.normalize()
    # Scale the unit vector to the correct distance
    scaled_vector = unit_vector ** distance
    # Add the scaled vector to the first bead to get coordinates for a fluorophore
    new_vector = scaled_vector + first_bead
    # Convert vector to tuple
    coords_out = (new_vector[0], new_vector[1], new_vector[2])
    return coords_out

def place_fluorophore(struct, dist, bead):
    """ Place bead located by pick_single_FPpoint in all models in structure """
    atoms = [x for x in struct.get_atoms()]
    single_fp_coord = pick_single_FPpoint(atoms, distance = dist, bead_name=bead)
    fp_atom = PDB.Atom.Atom(bead, single_fp_coord, 0, 1,' ', bead.ljust(3), 0, "C")
    fp_residue = PDB.Residue.Residue((' ', 0, ' '), 'SFP', 0)
    fp_residue.add(fp_atom)
    for model in struct:
        for chain in model:
            chain.add(fp_residue)


def load_distances_main(files_to_load_closed, 
                        files_to_load_open, 
                        fp_resnum, 
                        residues_to_check, 
                        residues_to_align, 
                        extra_distance, 
                        bead_name = atom_name, 
                        quiet = False):
    """This is where the magic happens."""
    # First up, load the structures from the list given as an argument
    open_structures = load_structures(files_to_load_open, quiet = quiet)
    closed_structures = load_structures(files_to_load_closed, quiet = quiet)

    # Now figure out which residues to check, if it's 'all'
    if residues_to_check == 'all':
        result = set()
        # Loop over loaded structures
        for struct in open_structures + closed_structures:
            for model in struct:
                # Store all the residues with alpha carbons from this model in a set
                model_resis = set()
                for residue in model.get_residues():
                    if "CA" in residue:
                        model_resis.add(residue.get_id()[1])
                # If we don't have a result yet, save this model's residues
                if result == set():
                    result = model_resis
                # Otherwise save the intersection of all residues
                else:
                    result.intersection_update(model_resis)
        residues_to_check = result


    # Now calculate the Rangefinder fluorophore location
    distances_dict = {}

    for n, structure in enumerate(open_structures):
        place_fluorophore(structure, extra_distance, bead_name)
    #     save_structure(structure, "open_fp" + str(n) + ".pdb")

    for n, structure in enumerate(closed_structures):
        place_fluorophore(structure, extra_distance, bead_name)
    #     save_structure(structure, "closed_fp" + str(n) + ".pdb")


    # # Now, measure all the interfluorophore distances for the models we just read
    open_distances = get_all_interfluorophore_distances(fp_resnum, 
                                                        residues_to_check, 
                                                        residues_to_align, 
                                                        open_structures, 
                                                        closed_structures, 
                                                        bead_name = bead_name, 
                                                        quiet = quiet)
    closed_distances = get_all_interfluorophore_distances(fp_resnum, 
                                                          residues_to_check, 
                                                          residues_to_align, 
                                                          closed_structures, 
                                                          open_structures, 
                                                          bead_name = bead_name, 
                                                          quiet = quiet)
    
    distances_dict["Rangefinder predictions"] = {"open": open_distances, "closed": closed_distances}

    return distances_dict

def calculate_fret_efficiency(r, r_nought):
    """Calculate the FRET efficiency at r for a fluorophore pair with Forster distance r_nought."""
    r_nm = float(r / 10.0)
    efficiency = 1 / (1 + (r_nm / r_nought)**6)
    return efficiency

def dynamic_range(efficiency_bound, efficiency_apo, fudge):
    """Predict the dynamic range from two efficiencies and the overlap fudge factor."""
    try:
        dynamic_range = (efficiency_bound - efficiency_apo) / ((1 - efficiency_bound) * (1 + (fudge * (efficiency_apo - 1))))
        # fudge ~ 0.70 for CFP-Venus or CFP-OG488
        # fudge ~ 0.82 for CFP-AF532
        return dynamic_range
    except ZeroDivisionError:
        return float('nan')

def distances_to_efficiencies(distances_dict, dyes):
    """Converts a dict of distances to one of efficiencies given the dict of dyes

    distances_dict should be indexed by method, conformation and then mutation, with the final level being a 
        list of distances, eg as produced by load_distances_main()
    dyes should have r_nought as the second element (index 1)
    Output is indexed by dyes, method, mutation and then conformation to aid in dynamic range generation
    Output is given as the mean, as efficiency is an ensemble average property"""

    efficiencies = {}
    for dye in dyes:
        r_nought = dye[1]
        efficiencies[dye] = {}
        for method,conformations in distances_dict.items():
            efficiencies[dye][method] = {}
            for conformation,mutations in conformations.items():
                for mutation,distances in mutations.items():
                    # We want to select by mutation before conformation in the next step, which means 
                    # changing things up here, but we don't want to lose objects if we've been here before
                    if mutation not in efficiencies[dye][method]:
                        efficiencies[dye][method][mutation] = {}
                    # Use list comprehension to calculate FRET efficiencies
                    efficiencies_list = [calculate_fret_efficiency(r, r_nought) for r in distances]
                    efficiencies[dye][method][mutation][conformation] = mean(efficiencies_list)

    return efficiencies

def efficiencies_to_dynamic_ranges(efficiencies, aligned_residues, advanced_mode):
    """Converts a dict of efficiencies to one of dynamic ranges, averaged over the entire set of structures

    distances_dict should be indexed by dye (third element should be fudge factor), mutation and then conformation, with 
        the final level being a list of efficiencies, eg as produced by distances_to_efficiencies()
    Output is indexed by conditions, which includes both dye and method, and then mutation, with each mutation's dynamic
        range"""


    dynamic_ranges = {}
    for dye,methods in efficiencies.items():
        fudge_factor = dye[2]
        dye_name = str(dye[0])
        for method,mutations in methods.items():
            condition = dye_name + " " + str(method)
            dynamic_ranges[condition] = {}
            dynamic_ranges[condition + " closed eff"] = {}
            dynamic_ranges[condition + " open eff"] = {}
            dynamic_ranges[condition + " closed dist"] = {}
            dynamic_ranges[condition + " open dist"] = {}
            for mutation,conformations in mutations.items():
                # print(method, mutation, conformations)
                # We need efficiencies from both conformations for this, so zip them together
                e1 = conformations["closed"]
                e2 = conformations["open"]
                if mutation in aligned_residues and not advanced_mode:
                    this_dynamic_range = 0.0
                else:
                    this_dynamic_range = dynamic_range(e1, e2, fudge_factor)
                dynamic_ranges[condition][mutation] = this_dynamic_range
    return dynamic_ranges

def dynamic_ranges_to_csv(in_dict, aligned_residues, filename):
    """Traverses the dynamic_ranges dict and generates a table in csv format"""
   
    # Reorganise the dict so that it's in the same order as the CSV
    organised_dict = {}
    for condition,mutations in in_dict.items():
        for mutation,dynamic_range in mutations.items():
            # Initialise the mutation dict iff we haven't been here before
            if mutation not in organised_dict:
                organised_dict[mutation] = {}
            organised_dict[mutation][condition] = dynamic_range
            # We want a column that indicates the residues to which we aligned structures
            # aligned_title = "aligned?"
            # if aligned_title not in organised_dict[mutation]: # only check aligned_residues once per mutation
            #     if mutation in aligned_residues:
            #         organised_dict[mutation][aligned_title] = 1
            #     else:
            #         organised_dict[mutation][aligned_title] = 0

    dict_to_csv("residue", organised_dict, filename)

def dict_to_csv(key_title, in_dict, filename):
    """Takes a dictionary and prints it to csv

    in_dict should be in format {key:{y_title[y]}}, with the title for the x values given as a parameter. 
    Each dictionary in in_dict is thus a data line, printed as key,y1,y2,y3...
    Header line is taken from each y_dict and sorted, supports missing elements, and is
    presented as key_title,y_title1,y_title2,y_title3..."""

    with open(filename, 'w') as output_file:
        # Get a sorted tuple of all y_titles
        y_titles = set()
        for key,y_dicts in in_dict.items():
            for y_title in y_dicts:
                y_titles.add(y_title)
        y_titles = tuple(sorted(y_titles))

        # Write the header line
        # key_title,y_title1,y_title2,y_title3 etc.
        header_line = str(key_title) + ","
        for title in y_titles:
            header_line += str(title) + ","
        header_line += "\n" # Finish with newline
        output_file.write(header_line) # Write to disk

        # Write the actual csv
        for key,y_dict in in_dict.items():
            this_line = str(key) + ","
            for y_title in y_titles:
                try:
                    y_value = y_dict[y_title]
                except KeyError:
                    pass                
                else:
                    this_line += str(y_value)
                this_line += ","
            this_line += "\n"
            output_file.write(this_line)
        
        print("Wrote predictions to file " + output_file.name)
    


# Read everything from PDB using the above functions
global_distances_dict = load_distances_main(structures_to_load_closed, 
                                            structures_to_load_open, 
                                            fp_resnum, 
                                            residues_to_check, 
                                            residues_to_align, 
                                            linker_length, 
                                            quiet = quiet)
    
# Iterate over distances_dict and dyes to generate an analogous dictionary of efficiencies
global_efficiencies = distances_to_efficiencies(global_distances_dict, global_dyes)  

# Now we have efficiencies in a reasonable format we can get dynamic ranges
global_dynamic_ranges = efficiencies_to_dynamic_ranges(global_efficiencies, residues_to_align, advanced)

# Finally, print out a csv with the results
dynamic_ranges_to_csv(global_dynamic_ranges, residues_to_align, args.output)
