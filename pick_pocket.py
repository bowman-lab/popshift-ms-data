import loos
import numpy as np

refstruct = '187L.pdb'
ligand_name = 'PXY'
contact_distances = [4, 5]  # angstroms

m = loos.createSystem(refstruct)
m_l = loos.selectAtoms(m, f'resname != "{ligand_name}"')
l = loos.selectAtoms(m, f'resname == "{ligand_name}"')
for c in contact_distances:
    w = m_l.within(c, l)
    res = w.splitByResidue()
    resids = [i[0].resid() for i in res]
    np.savetxt(f'pocket-resids-{c}.txt', resids, fmt='%i')
