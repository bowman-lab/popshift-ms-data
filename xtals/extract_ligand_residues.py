import loos
from pathlib import Path
import json

with open('roster.json') as f:
    roster = json.load(f)

cut = 4.0
bump = 0.25
bondscut = 1.65
ignores = set(['apo'])
assoc = {}
lign_assoc = {}
errors = False
for ligname, pdbs in roster.items():
    if ligname in ignores:
        print('ignoring', ligname)
        continue
    else:
        ligp = Path(ligname)
        for pdb_p in ligp.iterdir():
            model = loos.selectAtoms(loos.createSystem(str(pdb_p)), '!(hydrogen || resname == "HOH")')
            model.findBonds(bondscut)
            mols = model.splitByMolecule()
            a99 = loos.selectAtoms(mols[0], 'resid == 99')
            contacting = []
            workingcut = cut
            while len(contacting) < 1:
                contacting = [x[0].resname() for x in mols[1:]
                              if x.contactWith(workingcut, a99)]
                workingcut += bump

            if len(contacting) == 1:
                pdb_id = pdb_p.stem
                assoc[pdb_id] = contacting[0]
                lign_assoc[ligname] = contacting[0]
                print(ligname, pdb_id, 'workingcut:', workingcut)
            else:
                errors = True
                print(pdb_p, 'workingcut:', workingcut, contacting)
            break

if not errors:
    print('no errors!')
    with open('ligand-assoc.json', 'w') as f:
        json.dump(assoc, f, indent=4)
    with open('ligand-name-resname.json', 'w') as f:
        json.dump(lign_assoc, f, indent=4)
