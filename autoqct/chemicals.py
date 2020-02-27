import numpy as np

# Calculate the length of bond between two atoms
def bondLength(atom1, atom2):
    return np.sqrt(np.sum((atom1-atom2)**2))

# Calculate the angle between three atoms
def angles(atom1, atom2, atom3):
    v_1 = atom1 - atom2
    v_2 = atom3 - atom2

    v_1mag = np.sqrt(np.sum(v_1**2))
    v_1norm = v_1 / v_1mag

    v_2mag = np.sqrt(np.sum(v_2**2))
    v_2norm = v_2 / v_2mag

    return np.degrees(np.arccos(np.dot(v_1norm, v_2norm)))


# Given 2 hydrogens and 1 oxygen does this fit the bond lengths and angle criteria?
def isWater(atomH1, atomO, atomH2):
    l1 = bondLength(atomH1, atomO)
    l2 = bondLength(atomH2, atomO)
    theta = angles(atomH1, atomO, atomH2)
    # np.abs(1.043 - l1) <= 0.1 and np.abs(1.043 - l2) <= 0.1 and np.abs(theta-105.5) <= 10:
    return np.abs(1.043 - l1) <= 0.15 \
	   and np.abs(1.043 - l2) <= 0.15 \
	   and np.abs(theta - 105.5) <= 15

# Do atoms fit hydroxide bond distance?
def isHydroxide(atomH, atomO):
    l = bondLength(atomH, atomO)
    return np.abs(0.95-l) < 0.5)

# Do clusters fit the clustering criteria radius to water hydrogen atoms
def isClustered(ionX, waterH, clusterLength):
    for waterPair in waterH:
        if    bondLength(ionX, waterPair[0]) > clusterLength \
	   or bondLength(ionX, waterPair[1]) > clusterLength:
            return False

    return True

def isIon(ion):
    ##### Could Need Work #####
    return ion == 'Cl' or 'F' or 'Br' or 'I'

def isProton(protonHydrogen):
    ##### Could Need Work #####
    return True

