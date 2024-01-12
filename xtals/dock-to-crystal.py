from jug import TaskGenerator
from vina import Vina
import subprocess as sp

from pathlib import Path


@TaskGenerator
def dock_vina(box_center, box_size, exhaustiveness, receptor_path, ligand_path, output_path):
    v = Vina(sf_name='vina', cpu=exhaustiveness)
    v.set_receptor(str(receptor_path))
    v.set_ligand_from_file(str(ligand_path))
    v.compute_vina_maps(center=box_center, box_size=box_size)
    v.dock(exhaustiveness=exhaustiveness)
    v.write_poses(str(output_path), n_poses=1, overwrite=True)
    return True


# @TaskGenerator
def dock_smina(box_center, box_size, exhaustiveness, receptor_path, ligand_path, output_path):
    return sp.run(['smina', '--receptor', str(receptor_path), '--ligand', str(ligand_path),
                   '--center_x', f'{box_center[0]}',
                   '--center_y', f'{box_center[1]}',
                   '--center_z', f'{box_center[2]}',
                   '--size_x', f'{box_size[0]}',
                   '--size_y', f'{box_size[1]}',
                   '--size_z', f'{box_size[2]}',
                   '--exhaustiveness', f'{exhaustiveness}',
                   '--cpu', '1',
                   '--num_modes', '1',
                   '--out', str(output_path)])


receptor_p = Path('n-butylbenzene/186L.pdbqt')
ligands_dir = Path('../all_ligands')
out_dir_p = Path('186L-docks')
box_center = (0, 0, 0)
box_size = (12, 12, 12)
exhaust = 32
out_dir_p.mkdir(exist_ok=True, parents=True)
outdict = {}
for ligp in ligands_dir.rglob('*.pdbqt'):
    output_pose_p = out_dir_p/ligp.name
    print(ligp)
    dock_smina(
        box_center,
        box_size,
        exhaust,
        receptor_p,
        ligp,
        output_pose_p
    )
