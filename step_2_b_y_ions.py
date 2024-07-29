from expasy import aminoacid_mw, monoisotopic
from re import findall as refindall
from collections import defaultdict


def protein_weight(protein: str) -> float:
    # Sum the molecular weights of amino acids
    mass = sum(aminoacid_mw[aa] for aa in protein)
    # Adjust for peptide bonds
    mass -= (len(protein) - 1) * 18.010565
    return mass

