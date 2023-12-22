import loos
from loos import pyloos
import argparse
from pathlib import Path
import pickle
import numpy as np
from spyrmsd.rmsd import symmrmsd
from collections import defaultdict
import shlex

def atomic_group_to_adjacency(ag: loos.AtomicGroup):
    n = len(ag)
    id_to_order = {at.id(): i for i, at in enumerate(ag)}
    adjacency = np.zeros((n, n), dtype=int)
    # maybe there's a more efficient way to do this, but this seemed easier
    for i, at in enumerate(ag):
        # assuming that the atom IDs don't have a gap and start at 1
        bond_ids = at.getBonds()
        index_in_grp = [id_to_order[bond_id] for bond_id in bond_ids]
        adjacency[i][index_in_grp] = 1
    return adjacency


def align_by_bonds(A: loos.AtomicGroup, B: loos.AtomicGroup):
    a_id_order = {at.id(): i for i, at in enumerate(A)}
    b_id_order = {at.id(): i for i, at in enumerate(B)}


def get_atomic_nums(model_path:Path, informat='pdb'):
    from openbabel import openbabel
    mol = openbabel.OBMol()
    conv = openbabel.OBConversion()
    conv.SetInFormat(informat)
    conv.ReadFile(mol, str(model_path))
    atomic_nums = []
    for at in openbabel.OBMolAtomIter(mol):
        atomic_nums.append(at.GetAtomicNum())
    return np.array(atomic_nums, dtype=int)


def purify_expt_multiconf(atomic_group, preferred_loc='A', ligand_sel=None):
    purified_group = loos.AtomicGroup()
    if ligand_sel:
        ligand_ag = loos.selectAtoms(atomic_group, ligand_sel)
    else:
        ligand_ag = loos.AtomicGroup()
    for i in atomic_group:
        alt = i.altLoc()
        # residues that have no alternative confs will have altLoc() == ''
        if alt == '':
            purified_group.append(i)
        elif alt == preferred_loc:
            if not ligand_ag.contains(i):
                purified_group.append(i)
    return purified_group



def get_altloc_subgroups(ref_ag):
    ref_sub_groups = defaultdict(loos.AtomicGroup)
    if ref_ag[0].altLoc() == '':
        altloc = 'A'
        for i in ref_ag:
            ref_sub_groups[altloc].append(i)
    else:
        for i in ref_ag:
            ref_sub_groups[i.altLoc()].append(i)

def receptor_p_from_ligand(receptor_dir, ligand_path):
    return receptor_dir.joinpath(*ligand_path.parts[-2:]).with_suffix('.pdb')


parser = argparse.ArgumentParser()
parser.add_argument('ref_struct', type=Path, help='path to referece structure.')
parser.add_argument('sel_string', type=str, help='File or selection string to align and compute RMSD over.')
parser.add_argument('receptor_dir', type=Path, help='Path to dir of structures to compute RMSD using. will be rglobbed for pdbs.')
parser.add_argument('ligand_poses', type=Path, help='Path to ligand pose lists pickle')
parser.add_argument('--ligand-sel', type=str, help='Optionally select ligand for comparison.', default=None)
parser.add_argument('--ligand-mutual-sel', type=str, default='!hydrogen', help='If comparing ligand confs, subset' 
                    'reference ligand and docked ligand using this selection.')
parser.add_argument('--spyrmsd', action=argparse.BooleanOptionalAction, default=True,
                   help='If thrown, compute spyrmsd using software from Meli and Biggin, DOI:10.1186/s13321-020-00455-2')
# argv = shlex.split("--ligand-sel 'resname == \"MR3\"' " \
#       'xtals/1-methylpyrrole/2OU0.pdb ' \
#       'loos-pocket-sel-ca.txt ' \
#       't4l-2/binding/resect-1.0/tica-500-msm-2000-k-75-nframes-20/receptor ' \
#       't4l-2/binding/resect-1.0/tica-500-msm-2000-k-75-nframes-20/extracted_scores/12xsmina/1-methylpyrrole.pickle')

# args = parser.parse_args(argv)
args = parser.parse_args()

ref_struct = loos.createSystem(str(args.ref_struct))
# purify_expt_multiconf(ref_struct)
outdir = args.receptor_dir.parent/f'receptor-rmsds/{args.ligand_poses.stem}'
if not outdir.is_dir():
    outdir.mkdir(parents=True, exist_ok=True)
with args.ligand_poses.open('rb') as f:
    ligand_poseps = pickle.load(f)
psel = Path(args.sel_string)
if psel.is_file():
    sel = psel.read_text()
else:
    sel = args.sel_string

# this should land on a PDB except in the case of incorrect paths.
model = loos.createSystem(str(receptor_p_from_ligand(args.receptor_dir, ligand_poseps[0][0])))
samples_pocket = loos.selectAtoms(model, sel)
ref_pocket = purify_expt_multiconf(loos.selectAtoms(ref_struct, sel))

if args.ligand_sel:
    ligand_rmsds = []
    ligand_full = loos.createSystem(str(ligand_poseps[0][0].with_suffix('.pdb')))
    ligand_subset = loos.selectAtoms(ligand_full, args.ligand_mutual_sel)
    
    ref_ligand_subset = loos.selectAtoms(loos.selectAtoms(ref_struct, args.ligand_sel), args.ligand_mutual_sel)
    
    for at in ligand_subset:
        at.resname(ref_ligand_subset[0].resname())
    ref_lig_altloc_groups = defaultdict(loos.AtomicGroup)
    # deal with altlocs in ref.
    if ref_ligand_subset[0].altLoc() == '':
        altloc = 'A'
        if args.spyrmsd:
            for at in ref_ligand_subset:
                at.altLoc(altloc)
                ref_lig_altloc_groups[altloc].append(at)
        else:  # with regular rmsd calculations, ligand atoms need to be in order.
            reorder_atoms = ligand_subset.atomOrderMapFrom(ref_ligand_subset)
            for i in reorder_atoms:
                at = ref_ligand_subset[i]
                at.altLoc(altloc)
                ref_lig_altloc_groups[altloc].append(at)
    else:
        # We will need to sort twice. 
        # First extricate altlocs into different AGs, 
        # storing these in defaultdict.
        for i in ref_ligand_subset:
            ref_lig_altloc_groups[i.altLoc()].append(i)
        if not args.spyrmsd:  # If using spyrmsd, symmetry ops should resolve divergent atom orders.
        # Now reorder those AGs and overwrite the original, potentially out of order groups.
            for altloc in ref_lig_altloc_groups:
                original = ref_lig_altloc_groups[altloc]
                reorder_atoms = ligand_subset.atomOrderMapFrom(original)
                reordered = loos.AtomicGroup()
                for index in reorder_atoms:
                    reordered.append(original[index])
                ref_lig_altloc_groups[altloc] = reordered
    
           
    
    if args.spyrmsd:
        if not ligand_subset.hasBonds():
            ligand_subset.findBonds(1.85)
        # have to write ligand pdb to file so that obabel can 'convert' it.
        ligand_subset_path = args.ref_struct.parent / ('pose-' + args.ref_struct.name)
        ligand_subset_path.write_text(str(loos.PDB.fromAtomicGroup(ligand_subset)))
        pose_adjacency = atomic_group_to_adjacency(ligand_subset)
        ligpose_anums = get_atomic_nums(ligand_subset_path)
        if not ref_ligand_subset.hasBonds():
            for altloc in ref_lig_altloc_groups:
                ref_lig_altloc_groups[altloc].findBonds(1.85)
        ref_ligand_path = args.ref_struct.parent / ('ligand-' + args.ref_struct.name)
        ref_ligand_path.write_text(str(loos.PDB.fromAtomicGroup(ref_lig_altloc_groups['A'])))
        ref_anums = get_atomic_nums(ref_ligand_path) 
        ref_adjacency = atomic_group_to_adjacency(ref_lig_altloc_groups['A'])
    
    ref_lig_residues = ref_lig_altloc_groups['A'].splitByMolecule()
    if len(ref_lig_residues) != 1:
        print('Error: multiple ligand molecules selected.', args.ref_struct, f"'{args.ligand_sel}'")
        exit(1)

rmsds = []
for state_ix, state_lig_poseps in enumerate(ligand_poseps):
    if not state_lig_poseps[0].is_file():
        print('pose file', state_lig_poseps[0])
        raise FileNotFoundError
    vtraj = pyloos.VirtualTrajectory(
        *(pyloos.Trajectory(str(receptor_p_from_ligand(args.receptor_dir, lp)), model, subset=sel) 
          for lp in state_lig_poseps)
        ) 
    if args.ligand_sel:
        ligand_vtraj = pyloos.VirtualTrajectory(
                *(pyloos.Trajectory(str(lp.with_suffix('.pdb')), ligand_full, subset=args.ligand_mutual_sel) 
                  for lp in state_lig_poseps)
            )   
    for i, _ in enumerate(vtraj):
        xform = loos.XForm(ref_pocket.superposition(samples_pocket))
        ref_struct.applyTransform(xform)
        rmsd = samples_pocket.rmsd(ref_pocket)
        rmsds.append(rmsd)
        # print(state_lig_poseps[i])
        if args.ligand_sel:
            # this will be coordinate updated ligand_full subsetted by args.ligand_mutual_sel.
            ligand_pose = ligand_vtraj[i]
            if args.spyrmsd:
                altloc_rmsds = []
                coords = ligand_pose.getCoords()
                for altloc in ref_lig_altloc_groups:
                    ref_lig_altloc = ref_lig_altloc_groups[altloc]
                    coords_ref = ref_lig_altloc.getCoords()
                    lig_rmsd = symmrmsd(
                        coords_ref,
                        coords,
                        ref_anums,
                        ligpose_anums,
                        ref_adjacency,
                        pose_adjacency
                    )
                    # print('spyrmsd:', lig_rmsd)
                    altloc_rmsds.append(lig_rmsd)
                ligand_rmsds.append(min(altloc_rmsds))
            else:
                lig_rmsd = min(ref_lig_altloc_groups[altloc].rmsd(ligand_pose) 
                            for altloc in ref_lig_altloc_groups)
                ligand_rmsds.append(lig_rmsd)

# ref_pocket = loos.selectAtoms(vtraj._reference, sel)
# rmsds = [loos.selectAtoms(frame, sel).rmsd(ref_pocket) for frame in vtraj]
s = args.ref_struct.stem + '.npy'
np.save(str(outdir/s), rmsds)
if args.ligand_sel:
    if args.spyrmsd:
        outdir = args.receptor_dir.parent/f'noreorder-ligand-spyrmsds/{args.ligand_poses.stem}'
    else:
        outdir = args.receptor_dir.parent/f'ligand-rmsds/{args.ligand_poses.stem}'
    if not outdir.is_dir():
        outdir.mkdir(parents=True, exist_ok=True)
    s = args.ref_struct.stem + '.npy'
    np.save(str(outdir/s), ligand_rmsds)