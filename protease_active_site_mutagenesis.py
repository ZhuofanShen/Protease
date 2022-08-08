import argparse
import contextlib
from os.path import basename, isfile
import math
from pymol import cmd
from pyrosetta import *
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import ExtraRotamers, \
    IncludeCurrent, OperateOnResidueSubset, RestrictAbsentCanonicalAASRLT, \
    RestrictToRepackingRLT, PreventRepackingRLT
from pyrosetta.rosetta.core.pack import *
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.scoring.constraints import add_fa_constraints_from_cmdline, \
    AmbiguousConstraint, AngleConstraint, AtomPairConstraint
from pyrosetta.rosetta.core.scoring.func import CircularHarmonicFunc, FlatHarmonicFunc
from pyrosetta.rosetta.core.scoring.symmetry import SymmetricScoreFunction
from pyrosetta.rosetta.core.select.residue_selector import \
    AndResidueSelector, NotResidueSelector, OrResidueSelector, \
    InterGroupInterfaceByVectorSelector, \
    ChainSelector, ResidueIndexSelector, ResidueNameSelector
from pyrosetta.rosetta.core.simple_metrics.metrics import \
    InteractionEnergyMetric, RMSDMetric, TotalEnergyMetric
from pyrosetta.rosetta.protocols.constraint_generator import \
    AddConstraints, CoordinateConstraintGenerator
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, \
    AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.symmetry import SetupForSymmetryMover


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb', type=str)
    parser.add_argument('-xtal', type=str, required=True)
    parser.add_argument('-params', type=str, nargs='*')
    parser.add_argument('-cst', type=str)
    parser.add_argument('-enzdescst', type=str)
    parser.add_argument('-site', type=int, required=True, help='pdb numbering of the mutated site')
    parser.add_argument('-aa', type=str, required=True)
    parser.add_argument('-symm', '--symmetry', type=str)
    parser.add_argument('-no_coord', '--no_coordinate_csts', type=int, nargs='*', 
        help='pdb numbering of unconstrained sites')
    parser.add_argument('-cart', '--cartesian', action='store_true')
    parser.add_argument('-n', '--n_decoy', type=int, default=10)
    return parser.parse_args()

def initialize_pyrosetta(args):
    opts = '-ex1 -ex2 -use_input_sc -flip_HNQ -no_optH false'
    if args.params:
        opts += ' -extra_res_fa ' + ' '.join(args.params)
    if args.cst:
        opts += ' -constraints:cst_fa_file {}'.format(args.cst)
    if args.enzdescst:
        opts += ' -enzdes:cstfile {} -run:preserve_header'.format(args.enzdescst)
    init(opts)

def set_score_function(symmetry, cartesian):
    if symmetry:
        score_function = SymmetricScoreFunction()
        if cartesian:
            score_function.add_weights_from_file('ref2015_cart_cst')
        else:
            score_function.add_weights_from_file('ref2015_cst')
    else:
        if cartesian:
            score_function = create_score_function('ref2015_cart_cst')
        else:
            score_function = create_score_function('ref2015_cst')
    score_function.set_weight(ScoreType.fa_intra_rep_nonprotein, 0.545)
    score_function.set_weight(ScoreType.fa_intra_atr_nonprotein, 1)
    return score_function

def create_pose(args):
    pose = pose_from_pdb(args.pdb)
    xtal_pose = pose_from_pdb(args.xtal)
    return pose, xtal_pose

def set_fold_tree(pose, is_symmetric):
    protease_vector = ChainSelector('A').apply(pose)
    protease_length = protease_vector.count(True)
    monomer_length = protease_vector.u()
    if not is_symmetric:
        monomer_length = int(monomer_length / 2)
    not_protease_length = monomer_length - protease_length
    chain_B_length = ChainSelector('B').apply(pose).count(True)
    if chain_B_length != not_protease_length:
        substrate_length = chain_B_length
    else:
        substrate_length = 0
    fold_tree = FoldTree()
    fold_tree.add_edge(1, protease_length, -1)
    jumps = 0
    if substrate_length > 0:
        fold_tree.add_edge(1, protease_length + 1, 1)
        jumps = 1
        if substrate_length > 1: # and is_protein
            fold_tree.add_edge(protease_length + 1, protease_length + substrate_length, -1)
    for edge, wat in enumerate(range(protease_length + substrate_length + 1, monomer_length + 1)):
        if substrate_length > 0:
            fold_tree.add_edge(1, wat, edge + 2)
        else:
            fold_tree.add_edge(1, wat, edge + 1)
        jumps += 1
    if not is_symmetric:
        current_edge = edge + 3
        fold_tree.add_edge(1, monomer_length + 1, current_edge)
        jumps += 1
        fold_tree.add_edge(monomer_length + 1, monomer_length + protease_length, -1)
        if substrate_length > 0:
            current_edge += 1
            fold_tree.add_edge(monomer_length + 1, monomer_length + protease_length + 1, current_edge)
            jumps += 1
            if substrate_length > 1: # and is_protein
                fold_tree.add_edge(monomer_length + protease_length + 1, monomer_length + protease_length + substrate_length, -1)
        for edge2, wat in enumerate(range(monomer_length + protease_length + substrate_length + 1, 2 * monomer_length + 1)):
            if substrate_length > 0:
                fold_tree.add_edge(monomer_length + 1, wat, current_edge + edge2 + 1)
            else:
                fold_tree.add_edge(monomer_length + 1, wat, current_edge + edge2)
            jumps += 1
    pose.fold_tree(fold_tree)
    protease_selector = ResidueIndexSelector('1-' + str(protease_length) + ',' + str(monomer_length + 1) + '-' + str(monomer_length + protease_length))
    return protease_selector, jumps, substrate_length, monomer_length

def apply_enzyme_design_constraints(pose, enzdescst):
    enzdescst = AddOrRemoveMatchCsts()
    enzdescst.set_cst_action(ADD_NEW)
    enzdescst.apply(pose)

def apply_constraints(pose, xtal_ref_pdb, substrate_length, is_first_monomer, site):
    if is_first_monomer:
        pro_chain = 'A'
        if substrate_length > 0:
            subs_chain = 'B'
            wat_chain = 'C'
            wat_chain2 = 'F'
        else:
            subs_chain = '0'
            wat_chain = 'B'
            wat_chain2 = 'D'
    else:
        if substrate_length > 0:
            pro_chain = 'D'
            subs_chain = 'E'
            wat_chain = 'F'
            wat_chain2 = 'C'
        else:
            pro_chain = 'C'
            subs_chain = '0'
            wat_chain = 'D'
            wat_chain2 = 'B'
    # Define chain A atoms
    pose_idx_A1S = pose.pdb_info().pdb2pose(pro_chain, 1)
    atom_A1S_N = AtomID(pose.residue(pose_idx_A1S).atom_index("N"), pose_idx_A1S)
    atom_A1S_OG = AtomID(pose.residue(pose_idx_A1S).atom_index("OG"), pose_idx_A1S)
    pose_idx_A26T = pose.pdb_info().pdb2pose(pro_chain, 26)
    atom_A26T_N = AtomID(pose.residue(pose_idx_A26T).atom_index("N"), pose_idx_A26T)
    atom_A26T_CA = AtomID(pose.residue(pose_idx_A26T).atom_index("CA"), pose_idx_A26T)
    atom_A26T_C = AtomID(pose.residue(pose_idx_A26T).atom_index("C"), pose_idx_A26T)
    atom_A26T_O = AtomID(pose.residue(pose_idx_A26T).atom_index("O"), pose_idx_A26T)
    pose_idx_A41H = pose.pdb_info().pdb2pose(pro_chain, 41)
    atom_A41H_N = AtomID(pose.residue(pose_idx_A41H).atom_index("N"), pose_idx_A41H)
    atom_A41H_ND1 = AtomID(pose.residue(pose_idx_A41H).atom_index("ND1"), pose_idx_A41H)
    pose_idx_A119N = pose.pdb_info().pdb2pose(pro_chain, 119)
    atom_A119N_N = AtomID(pose.residue(pose_idx_A119N).atom_index("N"), pose_idx_A119N)
    atom_A119N_ND2 = AtomID(pose.residue(pose_idx_A119N).atom_index("ND2"), pose_idx_A119N)
    atom_A119N_OD1 = AtomID(pose.residue(pose_idx_A119N).atom_index("OD1"), pose_idx_A119N)
    pose_idx_A140F = pose.pdb_info().pdb2pose(pro_chain, 140)
    atom_A140F_C = AtomID(pose.residue(pose_idx_A140F).atom_index("C"), pose_idx_A140F)
    atom_A140F_O = AtomID(pose.residue(pose_idx_A140F).atom_index("O"), pose_idx_A140F)
    pose_idx_A142N = pose.pdb_info().pdb2pose(pro_chain, 142)
    atom_A142N_N = AtomID(pose.residue(pose_idx_A142N).atom_index("N"), pose_idx_A142N)
    atom_A142N_OD1 = AtomID(pose.residue(pose_idx_A142N).atom_index("OD1"), pose_idx_A142N)
    pose_idx_A143G = pose.pdb_info().pdb2pose(pro_chain, 143)
    atom_A143G_N = AtomID(pose.residue(pose_idx_A143G).atom_index("N"), pose_idx_A143G)
    atom_A143G_CA = AtomID(pose.residue(pose_idx_A143G).atom_index("CA"), pose_idx_A143G)
    atom_A143G_O = AtomID(pose.residue(pose_idx_A143G).atom_index("O"), pose_idx_A143G)
    pose_idx_A145C = pose.pdb_info().pdb2pose(pro_chain, 145)
    atom_A145C_N = AtomID(pose.residue(pose_idx_A145C).atom_index("N"), pose_idx_A145C)
    atom_A145C_CA = AtomID(pose.residue(pose_idx_A145C).atom_index("CA"), pose_idx_A145C)
    pose_idx_A164H = pose.pdb_info().pdb2pose(pro_chain, 164)
    atom_A164H_O = AtomID(pose.residue(pose_idx_A164H).atom_index("O"), pose_idx_A164H)    
    atom_A164H_ND1 = AtomID(pose.residue(pose_idx_A164H).atom_index("ND1"), pose_idx_A164H)
    pose_idx_A166E = pose.pdb_info().pdb2pose(pro_chain, 166)
    atom_A166E_N = AtomID(pose.residue(pose_idx_A166E).atom_index("N"), pose_idx_A166E)
    atom_A166E_CA = AtomID(pose.residue(pose_idx_A166E).atom_index("CA"), pose_idx_A166E)
    atom_A166E_C = AtomID(pose.residue(pose_idx_A166E).atom_index("C"), pose_idx_A166E)
    atom_A166E_O = AtomID(pose.residue(pose_idx_A166E).atom_index("O"), pose_idx_A166E)
    pose_idx_A187D = pose.pdb_info().pdb2pose(pro_chain, 187)
    atom_A187D_OD2 = AtomID(pose.residue(pose_idx_A187D).atom_index("OD2"), pose_idx_A187D)
    # Define substrate atoms
    pose_idx_B401 = pose.pdb_info().pdb2pose(subs_chain, 401)
    pose_idx_B505 = pose.pdb_info().pdb2pose(subs_chain, 505)
    atom_B505_N = None
    atom_B505_CA = None
    atom_B505_C = None
    atom_B505_O = None
    if pose_idx_B505 != 0:
        atom_B505_N = AtomID(pose.residue(pose_idx_B505).atom_index("N"), pose_idx_B505)
        atom_B505_N_str = '/505/N'
        atom_B505_CA = AtomID(pose.residue(pose_idx_B505).atom_index("CA"), pose_idx_B505)
        atom_B505_CA_str = '/505/CA'
        atom_B505_C = AtomID(pose.residue(pose_idx_B505).atom_index("C"), pose_idx_B505)
        atom_B505_C_str = '/505/C'
        atom_B505_O = AtomID(pose.residue(pose_idx_B505).atom_index("O"), pose_idx_B505)
        atom_B505_O_str = '/505/O'
    elif pose_idx_B401 != 0:
        with contextlib.suppress(RuntimeError):
            atom_B505_N = AtomID(pose.residue(pose_idx_B401).atom_index("NB5"), pose_idx_B401)
            atom_B505_N_str = '/401/NB5'
        with contextlib.suppress(RuntimeError):
            atom_B505_CA = AtomID(pose.residue(pose_idx_B401).atom_index("CAB5"), pose_idx_B401)
            atom_B505_CA_str = '/401/CAB5'
        with contextlib.suppress(RuntimeError):
            atom_B505_C = AtomID(pose.residue(pose_idx_B401).atom_index("CB5"), pose_idx_B401)
            atom_B505_C_str = '/401/CB5'
        with contextlib.suppress(RuntimeError):
            atom_B505_O = AtomID(pose.residue(pose_idx_B401).atom_index("OB5"), pose_idx_B401)
            atom_B505_O_str = '/401/OB5'
    pose_idx_B507 = pose.pdb_info().pdb2pose(subs_chain, 507)
    atom_B507_N = None
    atom_B507_CA = None
    atom_B507_C = None
    atom_B507_O = None
    atom_B507_CD = None
    atom_B507_NE2 = None
    if pose_idx_B507 != 0:
        atom_B507_N = AtomID(pose.residue(pose_idx_B507).atom_index("N"), pose_idx_B507)
        atom_B507_N_str = '/507/N'
        atom_B507_CA = AtomID(pose.residue(pose_idx_B507).atom_index("CA"), pose_idx_B507)
        atom_B507_CA_str = '/507/CA'
        atom_B507_C = AtomID(pose.residue(pose_idx_B507).atom_index("C"), pose_idx_B507)
        atom_B507_C_str = '/507/C'
        atom_B507_O = AtomID(pose.residue(pose_idx_B507).atom_index("O"), pose_idx_B507)
        atom_B507_O_str = '/507/O'
        atom_B507_CD = AtomID(pose.residue(pose_idx_B507).atom_index("CD"), pose_idx_B507)
        atom_B507_CD_str = '/507/CD'
        atom_B507_NE2 = AtomID(pose.residue(pose_idx_B507).atom_index("NE2"), pose_idx_B507)
        atom_B507_NE2_str = '/507/NE2'
    elif pose_idx_B401 != 0:
        with contextlib.suppress(RuntimeError):
            atom_B507_N = AtomID(pose.residue(pose_idx_B401).atom_index("NB7"), pose_idx_B401)
            atom_B507_N_str = '/401/NB7'
        with contextlib.suppress(RuntimeError):
            atom_B507_CA = AtomID(pose.residue(pose_idx_B401).atom_index("CAB7"), pose_idx_B401)
            atom_B507_CA_str = '/401/CAB7'
        with contextlib.suppress(RuntimeError):
            atom_B507_C = AtomID(pose.residue(pose_idx_B401).atom_index("CB7"), pose_idx_B401)
            atom_B507_C_str = '/401/CB7'
        with contextlib.suppress(RuntimeError):
            atom_B507_O = AtomID(pose.residue(pose_idx_B401).atom_index("OB7"), pose_idx_B401)
            atom_B507_O_str = '/401/OB7'
        with contextlib.suppress(RuntimeError):
            atom_B507_CD = AtomID(pose.residue(pose_idx_B401).atom_index("CDB7"), pose_idx_B401)
            atom_B507_CD_str = '/401/CDB7'
        with contextlib.suppress(RuntimeError):
            atom_B507_NE2 = AtomID(pose.residue(pose_idx_B401).atom_index("NEB7"), pose_idx_B401)
            atom_B507_NE2_str = '/401/NEB7'
    pose_idx_B509 = pose.pdb_info().pdb2pose(subs_chain, 509)
    atom_B509_N = None
    atom_B509_CA = None
    atom_B509_C = None
    atom_B509_O = None
    if pose_idx_B509 != 0:
        atom_B509_N = AtomID(pose.residue(pose_idx_B509).atom_index("N"), pose_idx_B509)
        atom_B509_N_str = '/509/N'
        atom_B509_CA = AtomID(pose.residue(pose_idx_B509).atom_index("CA"), pose_idx_B509)
        atom_B509_CA_str = '/509/CA'
        atom_B509_C = AtomID(pose.residue(pose_idx_B509).atom_index("C"), pose_idx_B509)
        atom_B509_C_str = '/509/C'
        atom_B509_O = AtomID(pose.residue(pose_idx_B509).atom_index("O"), pose_idx_B509)
        atom_B509_O_str = '/509/O'
    elif pose_idx_B401 != 0:
        with contextlib.suppress(RuntimeError):
            atom_B509_N = AtomID(pose.residue(pose_idx_B401).atom_index("NB9"), pose_idx_B401)
            atom_B509_N_str = '/401/NB9'
        with contextlib.suppress(RuntimeError):
            atom_B509_CA = AtomID(pose.residue(pose_idx_B401).atom_index("CAB9"), pose_idx_B401)
            atom_B509_CA_str = '/401/CAB9'
        with contextlib.suppress(RuntimeError):
            atom_B509_C = AtomID(pose.residue(pose_idx_B401).atom_index("CB9"), pose_idx_B401)
            atom_B509_C_str = '/401/CB9'
        with contextlib.suppress(RuntimeError):
            atom_B509_O = AtomID(pose.residue(pose_idx_B401).atom_index("OB9"), pose_idx_B401)
            atom_B509_O_str = '/401/OB9'
    # Define water oxygen atoms
    pose_idx_C601WAT = pose.pdb_info().pdb2pose(wat_chain, 601)
    if pose_idx_C601WAT != 0:
        atom_C601WAT_O = AtomID(pose.residue(pose_idx_C601WAT).atom_index("O"), pose_idx_C601WAT)
    pose_idx_C602WAT = pose.pdb_info().pdb2pose(wat_chain, 602)
    if pose_idx_C602WAT != 0:
        atom_C602WAT_O = AtomID(pose.residue(pose_idx_C602WAT).atom_index("O"), pose_idx_C602WAT)
    pose_idx_C603WAT = pose.pdb_info().pdb2pose(wat_chain, 603)
    if pose_idx_C603WAT != 0:
        atom_C603WAT_O = AtomID(pose.residue(pose_idx_C603WAT).atom_index("O"), pose_idx_C603WAT)
    pose_idx_C604WAT = pose.pdb_info().pdb2pose(wat_chain, 604)
    pose_idx_F604WAT = pose.pdb_info().pdb2pose(wat_chain2, 604)
    if pose_idx_C604WAT != 0:
        atom_C604WAT_O = AtomID(pose.residue(pose_idx_C604WAT).atom_index("O"), pose_idx_C604WAT)
        atom_F604WAT_O = AtomID(pose.residue(pose_idx_F604WAT).atom_index("O"), pose_idx_F604WAT)
    pose_idx_C605WAT = pose.pdb_info().pdb2pose(wat_chain, 605)
    pose_idx_F605WAT = pose.pdb_info().pdb2pose(wat_chain2, 605)
    if pose_idx_C605WAT != 0:
        atom_C605WAT_O = AtomID(pose.residue(pose_idx_C605WAT).atom_index("O"), pose_idx_C605WAT)
        atom_F605WAT_O = AtomID(pose.residue(pose_idx_F605WAT).atom_index("O"), pose_idx_F605WAT)
    
    # Load xtal symm pdb
    cmd.load(xtal_ref_pdb, 'xtal')
    # Add C601WAT constraints
    if pose_idx_C601WAT != 0:
        if site != 164:
            # harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + wat_chain + '/601/O','xtal//' + pro_chain + '/164/ND1'), 0.35, 0.3)
            harmonic_fc = FlatHarmonicFunc(3.0, 0.35, 0.3)
            pose.add_constraint(AtomPairConstraint(atom_C601WAT_O, atom_A164H_ND1, harmonic_fc))
        if site != 187:
            harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + wat_chain + '/601/O','xtal//' + pro_chain + '/187/OD2'), 0.35, 0.3)
            pose.add_constraint(AtomPairConstraint(atom_C601WAT_O, atom_A187D_OD2, harmonic_fc))
        harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + wat_chain + '/601/O','xtal//' + pro_chain + '/41/N'), 0.35, 0.3)
        pose.add_constraint(AtomPairConstraint(atom_C601WAT_O, atom_A41H_N, harmonic_fc))
        if site != 41:
            # harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + wat_chain + '/601/O','xtal//' + pro_chain + '/41/ND1'), 0.35, 0.3)
            harmonic_fc = FlatHarmonicFunc(3.0, 0.35, 0.3)
            pose.add_constraint(AtomPairConstraint(atom_C601WAT_O, atom_A41H_ND1, harmonic_fc))
    # Add C602WAT constraints
    if pose_idx_C602WAT != 0:
        harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + wat_chain + '/602/O','xtal//' + pro_chain + '/119/N'), 0.35, 0.3)
        pose.add_constraint(AtomPairConstraint(atom_C602WAT_O, atom_A119N_N, harmonic_fc))
        C602WAT_cst_2 = AmbiguousConstraint()
        harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + wat_chain + '/602/O','xtal//' + pro_chain + '/119/ND2'), 0.35, 0.3)
        C602WAT_cst_2.add_individual_constraint(AtomPairConstraint(atom_C602WAT_O, atom_A119N_ND2, harmonic_fc))
        harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + wat_chain + '/602/O','xtal//' + pro_chain + '/119/OD1'), 0.35, 0.3)
        C602WAT_cst_2.add_individual_constraint(AtomPairConstraint(atom_C602WAT_O, atom_A119N_OD1, harmonic_fc))
        pose.add_constraint(C602WAT_cst_2)
        harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + wat_chain + '/602/O','xtal//' + pro_chain + '/143/O'), 0.35, 0.3)
        pose.add_constraint(AtomPairConstraint(atom_C602WAT_O, atom_A143G_O, harmonic_fc))
        harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + wat_chain + '/602/O','xtal//' + pro_chain + '/26/O'), 0.35, 0.3)
        pose.add_constraint(AtomPairConstraint(atom_C602WAT_O, atom_A26T_O, harmonic_fc))
    # Add C603WAT constraints
    if pose_idx_C603WAT != 0:
        if site != 142:
            harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + wat_chain + '/603/O','xtal//' + pro_chain + '/142/OD1'), 0.35, 0.3)
            pose.add_constraint(AtomPairConstraint(atom_C603WAT_O, atom_A142N_OD1, harmonic_fc))
            if atom_B507_NE2:
                harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + wat_chain + '/603/O','xtal//' + subs_chain + atom_B507_NE2_str), 0.35, 0.3)
                pose.add_constraint(AtomPairConstraint(atom_C603WAT_O, atom_B507_NE2, harmonic_fc))
            if pose_idx_C604WAT != 0:
                harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + wat_chain + '/603/O','xtal//' + wat_chain + '/604/O'), 0.35, 0.3)
                pose.add_constraint(AtomPairConstraint(atom_C603WAT_O, atom_C604WAT_O, harmonic_fc))
            if pose_idx_C605WAT != 0:
                harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + wat_chain + '/603/O','xtal//' + wat_chain2 + '/605/O'), 0.35, 0.3)
                pose.add_constraint(AtomPairConstraint(atom_C603WAT_O, atom_F605WAT_O, harmonic_fc))
    # Add C604WAT constraints
    if pose_idx_C604WAT != 0:
        harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + wat_chain + '/604/O','xtal//' + pro_chain + '/142/N'), 0.35, 0.3)
        pose.add_constraint(AtomPairConstraint(atom_C604WAT_O, atom_A142N_N, harmonic_fc))
        if pose_idx_C605WAT != 0:
            harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + wat_chain + '/604/O','xtal//' + wat_chain2 + '/605/O'), 0.35, 0.3)
            pose.add_constraint(AtomPairConstraint(atom_C604WAT_O, atom_F605WAT_O, harmonic_fc))
            circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * cmd.angle('tmp','xtal//' + pro_chain + '/142/N','xtal//' + wat_chain + '/604/O','xtal//' + wat_chain2 + '/605/O'), 0.2)
            pose.add_constraint(AngleConstraint(atom_A142N_N, atom_C604WAT_O, atom_F605WAT_O, circular_harmonic_fc))
    # Add C605WAT constraints
    if pose_idx_C605WAT != 0:
        harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + wat_chain + '/605/O','xtal//' + pro_chain + '/1/N'), 0.35, 0.3)
        pose.add_constraint(AtomPairConstraint(atom_C605WAT_O, atom_A1S_N, harmonic_fc))
        harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + wat_chain + '/605/O','xtal//' + pro_chain + '/1/OG'), 0.35, 0.3)
        pose.add_constraint(AtomPairConstraint(atom_C605WAT_O, atom_A1S_OG, harmonic_fc))
        if pose_idx_C604WAT != 0:
            circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * cmd.angle('tmp','xtal//' + pro_chain + '/1/N','xtal//' + wat_chain + '/605/O','xtal//' + wat_chain2 + '/604/O'), 0.2)
            pose.add_constraint(AngleConstraint(atom_A1S_N, atom_C605WAT_O, atom_F604WAT_O, circular_harmonic_fc))
    # Add beta sheet constraints between A166E and B505
    if atom_B505_O:
        harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + pro_chain + '/166/N','xtal//' + subs_chain + atom_B505_O_str), 0.15, 0.3)
        pose.add_constraint(AtomPairConstraint(atom_A166E_N, atom_B505_O, harmonic_fc))
        circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * cmd.angle('tmp','xtal//' + pro_chain + '/166/CA','xtal//' + pro_chain + '/166/N','xtal//' + subs_chain + atom_B505_O_str), 0.2)
        pose.add_constraint(AngleConstraint(atom_A166E_CA, atom_A166E_N, atom_B505_O, circular_harmonic_fc))
        if atom_B505_C:
            circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * cmd.angle('tmp','xtal//' + pro_chain + '/166/N','xtal//' + subs_chain + atom_B505_O_str,'xtal//' + subs_chain + atom_B505_C_str), 0.2)
            pose.add_constraint(AngleConstraint(atom_A166E_N, atom_B505_O, atom_B505_C, circular_harmonic_fc))
    if atom_B505_N:
        harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + pro_chain + '/166/O','xtal//' + subs_chain + atom_B505_N_str), 0.15, 0.3)
        pose.add_constraint(AtomPairConstraint(atom_A166E_O, atom_B505_N, harmonic_fc))
        circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * cmd.angle('tmp','xtal//' + pro_chain + '/166/C','xtal//' + pro_chain + '/166/O','xtal//' + subs_chain + atom_B505_N_str), 0.2)
        pose.add_constraint(AngleConstraint(atom_A166E_C, atom_A166E_O, atom_B505_N, circular_harmonic_fc))
        if atom_B505_CA:
            circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * cmd.angle('tmp','xtal//' + pro_chain + '/166/O','xtal//' + subs_chain + atom_B505_N_str,'xtal//' + subs_chain + atom_B505_CA_str), 0.2)
            pose.add_constraint(AngleConstraint(atom_A166E_O, atom_B505_N, atom_B505_CA, circular_harmonic_fc))
    # Add H-bonding constraints between A164H and B507
    if atom_B507_N:
        harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + pro_chain + '/164/O','xtal//' + subs_chain + atom_B507_C_str), 0.15, 0.3)
        pose.add_constraint(AtomPairConstraint(atom_A164H_O, atom_B507_N, harmonic_fc))
        if atom_B507_CA:
            circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * cmd.angle('tmp','xtal//' + pro_chain + '/164/O','xtal//' + subs_chain + atom_B507_N_str,'xtal//' + subs_chain + atom_B507_CA_str), 0.2)
            pose.add_constraint(AngleConstraint(atom_A164H_O, atom_B507_N, atom_B507_CA, circular_harmonic_fc))
    # Add H-bonding constraints between A143G and B507
    if atom_B507_O:
        harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + pro_chain + '/143/N','xtal//' + subs_chain + atom_B507_O_str), 0.15, 0.3)
        pose.add_constraint(AtomPairConstraint(atom_A143G_N, atom_B507_O, harmonic_fc))
        circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * cmd.angle('tmp','xtal//' + pro_chain + '/143/CA','xtal//' + pro_chain + '/143/N','xtal//' + subs_chain + atom_B507_O_str), 0.2)
        pose.add_constraint(AngleConstraint(atom_A143G_CA, atom_A143G_N, atom_B507_O, circular_harmonic_fc))
        if atom_B507_C:
            circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * cmd.angle('tmp','xtal//' + pro_chain + '/143/N','xtal//' + subs_chain + atom_B507_O_str,'xtal//' + subs_chain + atom_B507_C_str), 0.2)
            pose.add_constraint(AngleConstraint(atom_A143G_N, atom_B507_O, atom_B507_C, circular_harmonic_fc))
    # Add H-bonding constraints between A145C and B507
    if atom_B507_O:
        harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + pro_chain + '/145/N','xtal//' + subs_chain + atom_B507_O_str), 0.15, 0.3)
        pose.add_constraint(AtomPairConstraint(atom_A145C_N, atom_B507_O, harmonic_fc))
        circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * cmd.angle('tmp','xtal//' + pro_chain + '/145/CA','xtal//' + pro_chain + '/145/N','xtal//' + subs_chain + atom_B507_O_str), 0.2)
        pose.add_constraint(AngleConstraint(atom_A145C_CA, atom_A145C_N, atom_B507_O, circular_harmonic_fc))
        if atom_B507_C:
            circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * cmd.angle('tmp','xtal//' + pro_chain + '/145/N','xtal//' + subs_chain + atom_B507_O_str,'xtal//' + subs_chain + atom_B507_C_str), 0.2)
            pose.add_constraint(AngleConstraint(atom_A145C_N, atom_B507_O, atom_B507_C, circular_harmonic_fc))
    # Add H-bonding constraints between A140F and B507
    if atom_B507_NE2:
        harmonic_fc = FlatHarmonicFunc(cmd.distance('tmp','xtal//' + pro_chain + '/140/O','xtal//' + subs_chain + atom_B507_NE2_str), 0.15, 0.3)
        pose.add_constraint(AtomPairConstraint(atom_A140F_O, atom_B507_NE2, harmonic_fc))
        circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * cmd.angle('tmp','xtal//' + pro_chain + '/140/C','xtal//' + pro_chain + '/140/O','xtal//' + subs_chain + atom_B507_NE2_str), 0.2)
        pose.add_constraint(AngleConstraint(atom_A140F_C, atom_A140F_O, atom_B507_NE2, circular_harmonic_fc))
        if atom_B507_CD:
            circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * cmd.angle('tmp','xtal//' + pro_chain + '/140/O','xtal//' + subs_chain + atom_B507_NE2_str,'xtal//' + subs_chain + atom_B507_CD_str), 0.2)
            pose.add_constraint(AngleConstraint(atom_A140F_O, atom_B507_NE2, atom_B507_CD, circular_harmonic_fc))
    # Add beta sheet constraints between A26T and B509
    if atom_B509_O:
        dst = cmd.distance('tmp','xtal//' + pro_chain + '/26/N','xtal//' + subs_chain + atom_B509_O_str)
        if dst < 4:
            harmonic_fc = FlatHarmonicFunc(dst, 0.15, 0.3)
            pose.add_constraint(AtomPairConstraint(atom_A26T_N, atom_B509_O, harmonic_fc))
            circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * cmd.angle('tmp','xtal//' + pro_chain + '/26/CA','xtal//' + pro_chain + '/26/N','xtal//' + subs_chain + atom_B509_O_str), 0.2)
            pose.add_constraint(AngleConstraint(atom_A26T_CA, atom_A26T_N, atom_B509_O, circular_harmonic_fc))
            if atom_B509_C:
                circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * cmd.angle('tmp','xtal//' + pro_chain + '/26/N','xtal//' + subs_chain + atom_B509_O_str,'xtal//' + subs_chain + atom_B509_C_str), 0.2)
                pose.add_constraint(AngleConstraint(atom_A26T_N, atom_B509_O, atom_B509_C, circular_harmonic_fc))
    if atom_B509_N:
        dst = cmd.distance('tmp','xtal//' + pro_chain + '/26/O','xtal//' + subs_chain + atom_B509_N_str)
        if dst < 4:
            harmonic_fc = FlatHarmonicFunc(dst, 0.15, 0.3)
            pose.add_constraint(AtomPairConstraint(atom_A26T_O, atom_B509_N, harmonic_fc))
            circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * cmd.angle('tmp','xtal//' + pro_chain + '/26/C','xtal//' + pro_chain + '/26/O','xtal//' + subs_chain + atom_B509_N_str), 0.2)
            pose.add_constraint(AngleConstraint(atom_A26T_C, atom_A26T_O, atom_B509_N, circular_harmonic_fc))
            if atom_B509_CA:
                circular_harmonic_fc = CircularHarmonicFunc(math.pi / 180 * cmd.angle('tmp','xtal//' + pro_chain + '/26/O','xtal//' + subs_chain + atom_B509_N_str,'xtal//' + subs_chain + atom_B509_CA_str), 0.2)
                pose.add_constraint(AngleConstraint(atom_A26T_O, atom_B509_N, atom_B509_CA, circular_harmonic_fc))

def add_coordinate_constraint(pose, xtal_pose, monomer_length, no_coordinate_csts=None):
    coord_cst_gen = CoordinateConstraintGenerator()
    coord_cst_gen.set_reference_pose(xtal_pose)
    if no_coordinate_csts:
        no_coord_cst_res_pose_indexes = list()
        for no_coord_cst_res in no_coordinate_csts:
            if no_coord_cst_res < 400:
                no_coord_cst_res_pose_idx = pose.pdb_info().pdb2pose('A', no_coord_cst_res)
                no_coord_cst_res_pose_indexes.append(str(no_coord_cst_res_pose_idx))
                no_coord_cst_res_pose_indexes.append(str(monomer_length + no_coord_cst_res_pose_idx))
            elif no_coord_cst_res > 400 and no_coord_cst_res < 600:
                no_coord_cst_res_pose_idx = pose.pdb_info().pdb2pose('B', no_coord_cst_res)
                no_coord_cst_res_pose_indexes.append(str(no_coord_cst_res_pose_idx))
                no_coord_cst_res_pose_indexes.append(str(monomer_length + no_coord_cst_res_pose_idx))
        coord_cst_gen.set_residue_selector(NotResidueSelector(ResidueIndexSelector(','.join(no_coord_cst_res_pose_indexes))))
    add_coord_cst = AddConstraints()
    add_coord_cst.add_generator(coord_cst_gen)
    add_coord_cst.apply(pose)

def create_residue_selector(mutation_selector, protease_selector):
    interface_selector = InterGroupInterfaceByVectorSelector()
    interface_selector.group1_selector(mutation_selector)
    interface_selector.group2_selector(NotResidueSelector(mutation_selector))
    movable_selector = OrResidueSelector(interface_selector, NotResidueSelector(protease_selector))
    repacking_selector = AndResidueSelector(movable_selector, NotResidueSelector(mutation_selector))
    static_selector = NotResidueSelector(movable_selector)

    interface_selector_2 = InterGroupInterfaceByVectorSelector()
    interface_selector_2.group1_selector(movable_selector)
    interface_selector_2.group2_selector(static_selector)
    minimization_selector = OrResidueSelector(movable_selector, interface_selector_2)
    minimization_vector = minimization_selector.apply(pose)
    return repacking_selector, static_selector, minimization_vector
    
def create_fast_relax_mover(score_function, cartesian, aa, mutation_selector, repacking_selector, \
        static_selector, minimization_vector, jumps):
    # task factory for the mutated pose
    task_factory = TaskFactory()
    task_factory.push_back(IncludeCurrent())
    # Mutated residue
    aa_force = RestrictAbsentCanonicalAASRLT()
    aa_force.aas_to_keep(aa)
    task_factory.push_back(OperateOnResidueSubset(aa_force, mutation_selector))
    # Repacking residues
    repack = RestrictToRepackingRLT()
    task_factory.push_back(OperateOnResidueSubset(repack, repacking_selector))
    # Immobile side chains
    prevent = PreventRepackingRLT()
    task_factory.push_back(OperateOnResidueSubset(prevent, static_selector))

    # move map for both the mutated pose and the reference pose
    move_map = MoveMap()
    move_map.set_bb(True)
    move_map.set_chi(minimization_vector)
    for jump in range(1, jumps + 1):
        move_map.set_jump(jump, True) # allow rigid body movement

    # fast relax mover for the mutated pose
    fast_relax = FastRelax()
    fast_relax.set_scorefxn(score_function)
    fast_relax.set_task_factory(task_factory)
    fast_relax.set_movemap(move_map)
    if cartesian:
        fast_relax.cartesian(True)

    return fast_relax

def apply_fast_relax(fast_relax, pose, n_decoy):
    for decoy in range(n_decoy):
        pose_copy = Pose(pose)
        fast_relax.apply(pose_copy)
        current_energy = calculate_energy(score_function, pose_copy)
        if decoy == 0 or current_energy < lowest_energy:
            lowest_energy = current_energy
            point_mutated_pose = pose_copy
    return point_mutated_pose

def calculate_energy(score_function, pose, selector=None, score_type=None):
    # Create the metric
    metric = TotalEnergyMetric()
    metric.set_scorefunction(score_function)
    if score_type:
        exec("metric.set_scoretype(ScoreType.{})".format(score_type))
    # Add the selector
    if selector:
        metric.set_residue_selector(selector)
    return metric.calculate(pose)

def calculate_interaction_energy(score_function, pose, selector_1, selector_2):
    # Create the metric
    metric = InteractionEnergyMetric()
    metric.set_scorefunction(score_function)
    # Add selectors
    metric.set_residue_selectors(selector_1, selector_2)
    return metric.calculate(pose)

def calculate_rmsd(pose, xtal_pose):
    # Create the RMSDMetric, setting 'xtal_pose' as the reference
    rmsd_metric = RMSDMetric()
    rmsd_metric.set_comparison_pose(xtal_pose)
    # rmsd_metric.set_residue_selector(protease_selector)
    # Use the RMSDMetirc to calculate the RMSD of 'pose'
    rmsd = rmsd_metric.calculate(pose)
    return rmsd

def mutagenesis(score_function, cartesian, site, monomer_length, aa, pose, xtal_pose, protease_selector, jumps, \
        n_decoy, prefix):
    site_pose_idx = pose.pdb_info().pdb2pose('A', site)
    mutation_selector = ResidueIndexSelector(str(site_pose_idx) + ',' + str(monomer_length + site_pose_idx))
    # create mutation, repacking, static, minimization residues selector
    repacking_selector, static_selector, minimization_vector \
            = create_residue_selector(mutation_selector, protease_selector)
    # create relax task factories
    fast_relax = create_fast_relax_mover(score_function, cartesian, aa, mutation_selector, \
            repacking_selector, static_selector, minimization_vector, jumps)
    # if the mutated pdb file exists, load this file
    if isfile(prefix + '_' + str(site) + aa + '.pdb'):
        point_mutated_pose = pose_from_pdb(prefix + '_' + str(site) + aa + '.pdb')
    # if not, apply the FastRelax mover to a copy of the pose
    else:
        point_mutated_pose = apply_fast_relax(fast_relax, pose, n_decoy)
        calculate_rmsd(point_mutated_pose, xtal_pose)
        point_mutated_pose.dump_pdb(prefix + '_' + str(site) + aa + '.pdb')
    return point_mutated_pose

def calculate_delta_G(score_function, point_mutated_pose, site, aa, protease_selector, monomer_length):
    # calculate delta total energy
    d_total_energy = calculate_energy(score_function, point_mutated_pose) \
            - calculate_energy(score_function, point_mutated_pose, score_type='coordinate_constraint')
    # calculate delta shell energy
    shell_res_pose_indexes = list()
    shell_res_list = [168,191,50,189,190,165,167,173,192,142,166,170,44,49,52,54,164,181,40,187,140,144,163,172,25,27,42]
    for shell_res in shell_res_list:
        shell_res_pose_idx = pose.pdb_info().pdb2pose('A', shell_res)
        shell_res_pose_indexes.append(str(shell_res_pose_idx))
        shell_res_pose_indexes.append(str(monomer_length + shell_res_pose_idx))
    shell_selector = ResidueIndexSelector(','.join(shell_res_pose_indexes))
    d_shell_energy = calculate_energy(score_function, point_mutated_pose, selector=shell_selector) \
            - calculate_energy(score_function, point_mutated_pose, selector=shell_selector, score_type='coordinate_constraint')
    # calculate delta substrate energy
    substrate_selector = NotResidueSelector(OrResidueSelector(protease_selector, ResidueNameSelector('TP3')))
    d_substrate_energy = calculate_energy(score_function, point_mutated_pose, selector=substrate_selector) \
            - calculate_energy(score_function, point_mutated_pose, selector=substrate_selector, score_type='coordinate_constraint')
    # calculate delta residue energy
    site_pose_idx = pose.pdb_info().pdb2pose('A', site)
    mutation_selector = ResidueIndexSelector(str(site_pose_idx) + ',' + str(monomer_length + site_pose_idx))
    d_residue_energy = calculate_energy(score_function, point_mutated_pose, selector=mutation_selector) \
            - calculate_energy(score_function, point_mutated_pose, selector=mutation_selector, score_type='coordinate_constraint')
    # calculate delta interaction energy
    delta_interaction_energy = calculate_interaction_energy(score_function, point_mutated_pose, protease_selector, substrate_selector)
    # calculate delta constraint energy
    delta_constraint_energy = calculate_energy(score_function, point_mutated_pose, score_type='atom_pair_constraint') + \
            calculate_energy(score_function, point_mutated_pose, score_type='angle_constraint') + \
            calculate_energy(score_function, point_mutated_pose, score_type='dihedral_constraint')
    with open(str(site) + '_' + aa + '.dat', 'w') as pf:
        pf.write(str(round(d_total_energy, 3)) + ',' + str(round(d_shell_energy, 3)) + ',' + \
                str(round(d_substrate_energy, 3)) + ',' + str(round(d_residue_energy, 3)) + ',' + \
                str(round(delta_interaction_energy, 3)) + ',' + str(round(delta_constraint_energy, 3)) + '\n')

def calculate_delta_G_2(score_function, point_mutated_pose, site, aa, protease_selector, monomer_length):
    # calculate delta total energy
    first_monomer_selector = ResidueIndexSelector("1-" + str(monomer_length))
    second_monomer_selector = ResidueIndexSelector(str(monomer_length + 1) + "-" + str(2 * monomer_length))
    d_total_energy_1 = calculate_energy(score_function, point_mutated_pose, selector=first_monomer_selector) \
            - calculate_energy(score_function, point_mutated_pose, selector=first_monomer_selector, score_type='coordinate_constraint')
    d_total_energy_2 = calculate_energy(score_function, point_mutated_pose, selector=second_monomer_selector) \
            - calculate_energy(score_function, point_mutated_pose, selector=second_monomer_selector, score_type='coordinate_constraint')
    # calculate delta shell energy
    shell_res_pose_indexes_1 = list()
    shell_res_pose_indexes_2 = list()
    shell_res_list = [168,191,50,189,190,165,167,173,192,142,166,170,44,49,52,54,164,181,40,187,140,144,163,172,25,27,42]
    for shell_res in shell_res_list:
        shell_res_pose_idx = pose.pdb_info().pdb2pose('A', shell_res)
        shell_res_pose_indexes_1.append(str(shell_res_pose_idx))
        shell_res_pose_indexes_2.append(str(monomer_length + shell_res_pose_idx))
    shell_selector_1 = ResidueIndexSelector(','.join(shell_res_pose_indexes_1))
    shell_selector_2 = ResidueIndexSelector(','.join(shell_res_pose_indexes_2))
    d_shell_energy_1 = calculate_energy(score_function, point_mutated_pose, selector=shell_selector_1) \
            - calculate_energy(score_function, point_mutated_pose, selector=shell_selector_1, score_type='coordinate_constraint')
    d_shell_energy_2 = calculate_energy(score_function, point_mutated_pose, selector=shell_selector_2) \
            - calculate_energy(score_function, point_mutated_pose, selector=shell_selector_2, score_type='coordinate_constraint')
    # calculate delta substrate energy
    substrate_selector = NotResidueSelector(OrResidueSelector(protease_selector, ResidueNameSelector('TP3')))
    substrate_selector_1 = AndResidueSelector(substrate_selector, first_monomer_selector)
    substrate_selector_2 = AndResidueSelector(substrate_selector, second_monomer_selector)
    d_substrate_energy_1 = calculate_energy(score_function, point_mutated_pose, selector=substrate_selector_1) \
            - calculate_energy(score_function, point_mutated_pose, selector=substrate_selector_1, score_type='coordinate_constraint')
    d_substrate_energy_2 = calculate_energy(score_function, point_mutated_pose, selector=substrate_selector_2) \
            - calculate_energy(score_function, point_mutated_pose, selector=substrate_selector_2, score_type='coordinate_constraint')
    # calculate delta residue energy
    site_pose_idx = pose.pdb_info().pdb2pose('A', site)
    mutation_selector_1 = ResidueIndexSelector(str(site_pose_idx))
    mutation_selector_2 = ResidueIndexSelector(str(monomer_length + site_pose_idx))
    d_residue_energy_1 = calculate_energy(score_function, point_mutated_pose, selector=mutation_selector_1) \
            - calculate_energy(score_function, point_mutated_pose, selector=mutation_selector_1, score_type='coordinate_constraint')
    d_residue_energy_2 = calculate_energy(score_function, point_mutated_pose, selector=mutation_selector_2) \
            - calculate_energy(score_function, point_mutated_pose, selector=mutation_selector_2, score_type='coordinate_constraint')
    # calculate delta interaction energy
    protease_selector_1 = AndResidueSelector(protease_selector, first_monomer_selector)
    protease_selector_2 = AndResidueSelector(protease_selector, second_monomer_selector)
    delta_interaction_energy_1 = calculate_interaction_energy(score_function, point_mutated_pose, protease_selector_1, substrate_selector_1)
    delta_interaction_energy_2 = calculate_interaction_energy(score_function, point_mutated_pose, protease_selector_2, substrate_selector_2)
    # calculate delta constraint energy
    delta_constraint_energy_1 = calculate_energy(score_function, point_mutated_pose, selector=first_monomer_selector, score_type='atom_pair_constraint') \
            + calculate_energy(score_function, point_mutated_pose, selector=first_monomer_selector, score_type='angle_constraint') \
            + calculate_energy(score_function, point_mutated_pose, selector=first_monomer_selector, score_type='dihedral_constraint')
    delta_constraint_energy_2 = calculate_energy(score_function, point_mutated_pose, selector=second_monomer_selector, score_type='atom_pair_constraint') \
            + calculate_energy(score_function, point_mutated_pose, selector=second_monomer_selector, score_type='angle_constraint') \
            + calculate_energy(score_function, point_mutated_pose, selector=second_monomer_selector, score_type='dihedral_constraint')
    with open(str(site) + '_' + aa + '.dat', 'w') as pf:
        pf.write(str(round(d_total_energy_1, 3)) + ',' + str(round(d_shell_energy_1, 3)) + ',' + \
                str(round(d_substrate_energy_1, 3)) + ',' + str(round(d_residue_energy_1, 3)) + ',' + \
                str(round(delta_interaction_energy_1, 3)) + ',' + str(round(delta_constraint_energy_1, 3)) + '\n')
        pf.write(str(round(d_total_energy_2, 3)) + ',' + str(round(d_shell_energy_2, 3)) + ',' + \
                str(round(d_substrate_energy_2, 3)) + ',' + str(round(d_residue_energy_2, 3)) + ',' + \
                str(round(delta_interaction_energy_2, 3)) + ',' + str(round(delta_constraint_energy_2, 3)) + '\n')


if __name__ == '__main__':
    args = parse_arguments()
    initialize_pyrosetta(args)
    score_function = set_score_function(args.symmetry, args.cartesian)
    pose, xtal_pose = create_pose(args)
    protease_selector, jumps, substrate_length, monomer_length = set_fold_tree(pose, args.symmetry)
    if args.symmetry:
        symmetry_mover = SetupForSymmetryMover(args.symmetry)
        symmetry_mover.apply(pose)
        symmetry_mover.apply(xtal_pose)
        xtal_ref_pdb = args.xtal[:-3] + 'symm.pdb'
        if not isfile(xtal_ref_pdb):
            score_function(xtal_pose)
            xtal_pose.dump_pdb(xtal_ref_pdb)
    else:
        xtal_ref_pdb = args.xtal
    if args.cst:
        add_fa_constraints_from_cmdline(pose, score_function)
    else:
        apply_constraints(pose, xtal_ref_pdb, substrate_length, True, args.site)
        if not args.symmetry:
                apply_constraints(pose, xtal_ref_pdb, substrate_length, False, args.site)
    if args.enzdescst:
        apply_enzyme_design_constraints(pose, args.enzdescst)
    add_coordinate_constraint(pose, xtal_pose, monomer_length, args.no_coordinate_csts)
    prefix = basename(args.pdb).replace('.pdb', '').replace('.gz', '')
    point_mutated_pose = mutagenesis(score_function, args.cartesian, args.site, monomer_length, args.aa, \
            pose, xtal_pose, protease_selector, jumps, args.n_decoy, prefix)
    if args.symmetry:
        calculate_delta_G(score_function, point_mutated_pose, args.site, args.aa, protease_selector, monomer_length)
    else:
        calculate_delta_G_2(score_function, point_mutated_pose, args.site, args.aa, protease_selector, monomer_length)
