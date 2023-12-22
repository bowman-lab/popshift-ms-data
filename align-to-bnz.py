from pathlib import Path
from pymol import cmd
from openbabel import pybel


ligand_dir = Path('all_ligands')
ligand_ps = list(ligand_dir.glob('*.pdbqt'))

refname = 'isobutylbenzene'
ref = cmd.load(f'xtals/{refname}/184L.pdb')
ref_res_name = 'I4B'
cmd.select('targ', f'resn {ref_res_name}')
cmd.create(refname, 'targ')
for ligp in ligand_ps:
    lign = ligp.stem
    if lign == refname:
        continue
    cmd.load(str(ligp), lign)
    cmd.align(lign, refname)
    out_pdb_p = ligp.with_name(ligp.stem + '-aligned.pdb')
    cmd.save(str(out_pdb_p), lign)
    # the next is needed, LOL.
    lig_pdb = next(pybel.readfile('pdb', str(out_pdb_p)))
    out_pdb_p.with_suffix('.sdf').write_text(lig_pdb.write('sdf'))
