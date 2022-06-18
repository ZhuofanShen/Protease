import argparse
from os.path import basename, isfile
from pyrosetta import *
from pyrosetta.rosetta.core.chemical import ResidueProperty
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import ExtraRotamers, \
    IncludeCurrent, OperateOnResidueSubset, RestrictAbsentCanonicalAASRLT, \
    RestrictToRepackingRLT, PreventRepackingRLT
from pyrosetta.rosetta.core.pack import *
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.scoring.constraints import *
from pyrosetta.rosetta.core.scoring.symmetry import SymmetricScoreFunction
from pyrosetta.rosetta.core.select.residue_selector import \
    AndResidueSelector, ChainSelector, ResidueIndexSelector, \
    NotResidueSelector, OrResidueSelector, \
    InterGroupInterfaceByVectorSelector, NeighborhoodResidueSelector
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
    parser.add_argument('-xtal', type=str, require=True)
    parser.add_argument('-ft', '--fold_tree', type=str, nargs='*', require=True)
    parser.add_argument('-params', type=str, nargs='*')
    parser.add_argument('-cst', type=str)
    parser.add_argument('-enzdescst', type=str, require=True)
    parser.add_argument('-site', type=int, require=True)
    parser.add_argument('-aa', type=str, require=True)
    parser.add_argument('-symm', '--symmetry', type=str, require=True)
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

def set_symmetric_score_function():
    # score_function = create_score_function('ref2015_cst')
    score_function = SymmetricScoreFunction()
    score_function.add_weights_from_file('ref2015_cst')
    score_function.set_weight(ScoreType.fa_intra_rep_nonprotein, 0.545)
    score_function.set_weight(ScoreType.fa_intra_atr_nonprotein, 1)
    return score_function

def create_pose(args):
    pose = pose_from_pdb(args.pdb)
    xtal_pose = pose_from_pdb(args.xtal)
    return pose, xtal_pose

def set_fold_tree(edge_list):
    # 'FOLD_TREE  EDGE 1 306 -1  EDGE 1 307 1  EDGE 307 316 -1 '
    # 'FOLD_TREE  EDGE 1 306 -1  EDGE 1 307 1 '
    fold_tree = FoldTree()
    for edge_str in edge_list:
        edge_num = edge_str.split(',')
        if len(edge_num) is 3:
            fold_tree.add_edge(int(edge_num[0]), int(edge_num[1]), int(edge_num[2]))
        elif len(edge_num) is 4:
            fold_tree.add_edge(int(edge_num[0]), int(edge_num[1]), edge_num[2], edge_num[3])
        else:
            raise Exception("The number of arguments in add_edge function should be 3 or 4")
    return fold_tree

def apply_symmetry(pose, symmetry):
    symmetry_mover = SetupForSymmetryMover(symmetry)
    symmetry_mover.set_keep_pdb_info_labels(True)
    symmetry_mover.apply(pose)

def get_protease_selector(pose):
    # fold_tree = pose.fold_tree().to_string().split('EDGE')
    # protease_index = fold_tree[1].strip(' ').split(' ')
    # protease_selector = ResidueIndexSelector(protease_index[0] + '-' + protease_index[1])
    protease_selector = ResidueIndexSelector('1-306,317-622')
    #protease_selector = ResidueIndexSelector('1-306,308-613')
    return protease_selector

def add_coordinate_constraint(xtal_pose):
    coord_cst_gen = CoordinateConstraintGenerator()
    coord_cst_gen.set_reference_pose(xtal_pose)
    coord_cst_gen.set_residue_selector(get_protease_selector(pose))
    add_coord_cst = AddConstraints()
    add_coord_cst.add_generator(coord_cst_gen)
    return add_coord_cst

def apply_constraints(args, pose, score_function):
    if args.cst:
        add_fa_constraints_from_cmdline(pose, score_function)
    if args.enzdescst:
        enzdes_cst = AddOrRemoveMatchCsts()
        enzdes_cst.set_cst_action(ADD_NEW)
        enzdes_cst.apply(pose)

def create_residue_selector(site):
    mutation_selector = ResidueIndexSelector(site)
    substrate_selector = NotResidueSelector(get_protease_selector(pose))
    interface_selector = InterGroupInterfaceByVectorSelector()
    interface_selector.group1_selector(mutation_selector)
    interface_selector.group2_selector(NotResidueSelector(mutation_selector))
    movable_selector = OrResidueSelector(interface_selector, substrate_selector)
    repacking_selector = AndResidueSelector(movable_selector, NotResidueSelector(mutation_selector))
    static_selector = NotResidueSelector(movable_selector)

    interface_selector_2 = InterGroupInterfaceByVectorSelector()
    interface_selector_2.group1_selector(movable_selector)
    interface_selector_2.group2_selector(static_selector)
    minimization_selector = OrResidueSelector(movable_selector, interface_selector_2)
    # c_term_selector = ResidueIndexSelector('200-306,516-622')
    # minimization_selector = AndResidueSelector(min_selector_temp, NotResidueSelector(c_term_selector))

    return mutation_selector, repacking_selector, movable_selector, static_selector, minimization_selector
    
def create_fast_relax_mover(score_function, aa, mutation_selector, repacking_selector, movable_selector, \
        static_selector, minimization_selector):
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

    # task factory for the reference pose
    ref_task_factory = TaskFactory()
    ref_task_factory.push_back(IncludeCurrent())
    repack_2 = RestrictToRepackingRLT()
    ref_task_factory.push_back(OperateOnResidueSubset(repack_2, movable_selector))
    prevent_2 = PreventRepackingRLT()
    ref_task_factory.push_back(OperateOnResidueSubset(prevent_2, static_selector))

    # move map for both the mutated pose and the reference pose
    move_map = MoveMap()
    move_map.set_bb(True)
    minimization_vector = minimization_selector.apply(pose)
    move_map.set_chi(minimization_vector)
    move_map.set_jump(1, True) # allow substrate to make rigid body movement

    # fast relax mover for the mutated pose
    fast_relax = FastRelax()
    fast_relax.set_scorefxn(score_function)
    fast_relax.set_task_factory(task_factory)
    fast_relax.set_movemap(move_map)

    # fast relax mover for the reference pose
    ref_fast_relax = FastRelax()
    ref_fast_relax.set_scorefxn(score_function)
    ref_fast_relax.set_task_factory(ref_task_factory)
    ref_fast_relax.set_movemap(move_map)
    return fast_relax, ref_fast_relax

def apply_fast_relax(fast_relax, pose, add_coord_cst, n_decoy):
    for decoy in range(n_decoy):
        pose_copy = Pose(pose) # deep copy pose
        add_coord_cst.apply(pose_copy)
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
    rmsd_metric.set_residue_selector(get_protease_selector(pose))
    # Use the RMSDMetirc to calculate the RMSD of 'pose'
    rmsd = rmsd_metric.calculate(pose)
    return rmsd

def mutagenesis(score_function, site, aa, pose, xtal_pose, add_coord_cst, prefix):
    # create mutation, repacking, static, minimization residues selector
    mutation_selector, repacking_selector, movable_selector, static_selector, \
            minimization_selector = create_residue_selector(site)
    # create relax task factories
    fast_relax, ref_fast_relax = create_fast_relax_mover(score_function, aa, mutation_selector, \
            repacking_selector, movable_selector, static_selector, minimization_selector)
    # if the mutated pdb file exists, load this file
    if isfile(prefix + '_' + str(site) + aa + '.pdb'):
        point_mutated_pose = pose_from_pdb(prefix + '_' + str(site) + aa + '.pdb')
    # if not, apply the FastRelax mover to a copy of the pose
    else:
        point_mutated_pose = apply_fast_relax(fast_relax, pose, add_coord_cst, n_decoy=10)
        calculate_rmsd(point_mutated_pose, xtal_pose)
        point_mutated_pose.dump_pdb(prefix + '_' + str(site) + aa + '.pdb')
    # if the reference pdb file exists, load this file
    if isfile(prefix + '_' + str(site) + aa + '.ref.pdb'):
        ref_pose = pose_from_pdb(prefix + '_' + str(site) + aa + '.ref.pdb')
    # if not, apply the FastRelax mover to a copy of the pose
    else:
        ref_pose = apply_fast_relax(ref_fast_relax, pose, add_coord_cst, n_decoy=10)
        calculate_rmsd(ref_pose, xtal_pose)
        ref_pose.dump_pdb(prefix + '_' + str(site) + aa + '.ref.pdb')
    return point_mutated_pose, ref_pose

def calculate_ddG(score_function, point_mutated_pose, ref_pose, site):
    # calculate delta delta total energy
    dd_total_energy = calculate_energy(score_function, point_mutated_pose) \
            - calculate_energy(score_function, point_mutated_pose, score_type='coordinate_constraint') \
            - calculate_energy(score_function, ref_pose) \
            + calculate_energy(score_function, ref_pose, score_type='coordinate_constraint')
    # calculate delta delta shell energy
    shell_selector = ResidueIndexSelector('168,191,50,189,190,165,167,173,192,142,166,170,44,49,52,54,164,181,40,187,140,144,163,172,25,27,42')
    dd_shell_energy = calculate_energy(score_function, point_mutated_pose, selector=shell_selector) \
            - calculate_energy(score_function, point_mutated_pose, selector=shell_selector, score_type='coordinate_constraint') \
            - calculate_energy(score_function, ref_pose, selector=shell_selector) \
            + calculate_energy(score_function, ref_pose, selector=shell_selector, score_type='coordinate_constraint')
    # calculate delta delta substrate energy
    substrate_selector = NotResidueSelector(get_protease_selector(pose))
    dd_substrate_energy = calculate_energy(score_function, point_mutated_pose, selector=substrate_selector) \
            - calculate_energy(score_function, point_mutated_pose, selector=substrate_selector, score_type='coordinate_constraint') \
            - calculate_energy(score_function, ref_pose, selector=substrate_selector) \
            + calculate_energy(score_function, ref_pose, selector=substrate_selector, score_type='coordinate_constraint')
    # calculate delta delta residue energy
    mutation_selector = ResidueIndexSelector(site)
    dd_residue_energy = calculate_energy(score_function, point_mutated_pose, selector=mutation_selector) \
            - calculate_energy(score_function, point_mutated_pose, selector=mutation_selector, score_type='coordinate_constraint') \
            - calculate_energy(score_function, ref_pose, selector=mutation_selector) \
            + calculate_energy(score_function, ref_pose, selector=mutation_selector, score_type='coordinate_constraint')
    # calculate delta interaction energy
    protease_selector = get_protease_selector(pose)
    substrate_selector = NotResidueSelector(protease_selector)
    delta_interaction_energy = calculate_interaction_energy(score_function, point_mutated_pose, protease_selector, substrate_selector) \
            - calculate_interaction_energy(score_function, ref_pose, protease_selector, substrate_selector)
    # calculate delta constraint energy
    delta_constraint_energy = calculate_energy(score_function, point_mutated_pose, score_type='atom_pair_constraint') + \
            calculate_energy(score_function, point_mutated_pose, score_type='angle_constraint') + \
            calculate_energy(score_function, point_mutated_pose, score_type='dihedral_constraint') - \
            calculate_energy(score_function, ref_pose, score_type='atom_pair_constraint') - \
            calculate_energy(score_function, ref_pose, score_type='angle_constraint') - \
            calculate_energy(score_function, ref_pose, score_type='dihedral_constraint')
    return dd_total_energy, dd_shell_energy, dd_substrate_energy, dd_residue_energy, \
            delta_interaction_energy, delta_constraint_energy

def write_csv(site, aa, dd_total_energy, dd_shell_energy, dd_substrate_energy, \
        dd_residue_energy, delta_interaction_energy, delta_constraint_energy):
    with open(str(site) + '_' + aa + '.dat', 'w') as pf:
        pf.write(str(round(dd_total_energy, 3)) + ',' + str(round(dd_shell_energy, 3)) + ',' + \
                str(round(dd_substrate_energy, 3)) + ',' + str(round(dd_residue_energy, 3)) + ',' + \
                str(round(delta_interaction_energy, 3)) + ',' + str(round(delta_constraint_energy, 3)))

if __name__ == '__main__':
    args = parse_arguments()
    initialize_pyrosetta(args)
    score_function = set_symmetric_score_function()
    pose, xtal_pose = create_pose(args)
    fold_tree = set_fold_tree(args.fold_tree)
    pose.fold_tree(fold_tree)
    xtal_pose.fold_tree(fold_tree)
    apply_symmetry(pose)
    apply_symmetry(xtal_pose)
    add_coord_cst = add_coordinate_constraint(xtal_pose)
    apply_constraints(args, pose, score_function)
    prefix = basename(args.pdb).replace('.pdb', '').replace('.gz', '')
    point_mutated_pose, ref_pose = mutagenesis(score_function, args.site, args.aa, pose, xtal_pose, \
            add_coord_cst, prefix)
    dd_total_energy, dd_shell_energy, dd_substrate_energy, dd_residue_energy, delta_interaction_energy, \
            delta_constraint_energy = calculate_ddG(score_function, point_mutated_pose, ref_pose, args.site)
    write_csv(args.site, args.aa, dd_total_energy, dd_shell_energy, dd_substrate_energy, dd_residue_energy, \
            delta_interaction_energy, delta_constraint_energy)
