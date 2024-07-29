from expasy import aminoacid_mw, monoisotopic
from re import findall as refindall
from collections import defaultdict


def protein_weight(protein: str) -> float:
    # Sum the molecular weights of amino acids
    mass = sum(aminoacid_mw[aa] for aa in protein)
    # Adjust for peptide bonds
    mass -= (len(protein) - 1) * 18.010565
    return mass

def cal_b_y_ion_mass(peptide):
    """
    Takes the peptide sequence and returns a dataframe of b and y ion masses.
    :param peptide: The peptide sequence
    :return: DataFrame with b and y ion masses
    """
    peptide = peptide.upper()
    lst = []
    # Compute b and y ion masses
    for i in range(1, len(peptide) + 1):
        b_ion = peptide[:i]
        y_ion = peptide[i-1:]

        b_mass = protein_weight(b_ion) + 1.007825 - 18.010565
        y_mass = protein_weight(y_ion) + 19.01839 - 18.010565  # 19.01839 accounts for the proton and H2O

        lst.append(b_mass)
        lst.append(y_mass)

    return lst