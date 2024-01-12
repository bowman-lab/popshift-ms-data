from pathlib import Path
import requests
import json


rcsb_url = 'https://files.rcsb.org/view/'
roster = {}
just_atoms = False
with open('roster.txt') as f:
    for line in f:
        split = line.split()
        roster[split[0]] = split[1:]

with open('roster.json', 'w') as f:
    json.dump(roster, f)

for ligname, pdb_ids in roster.items():
    ligdir = Path(ligname)
    if not ligdir.is_dir():
        ligdir.mkdir()
    for pdb_id in pdb_ids:
        pdb_p = ligdir / f'{pdb_id}.pdb'
        if not pdb_p.exists():
            # Get PDB text from RCSB website
            pdb_req = requests.get(f'{rcsb_url}{pdb_id}.pdb')
            try:
                pdb_req.raise_for_status()
                pdb_text = pdb_req.text
                if just_atoms:
                    with pdb_p.open() as pdb_f:
                        # Extract ATOM and TER lines
                        for line in pdb_text.split('\n'):
                            if line[:4] in ['ATOM', 'TER ']:
                                pdb_f.write(f'{line}\n')
                else:
                    pdb_p.write_text(pdb_text)
            except:
                print(pdb_id, ': request failed.')