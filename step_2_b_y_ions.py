from expasy import aminoacid, monoisotopic
from re import findall as refindall

def molecular_weight(molecule: str) -> float:
    """
    the function to calculate the molecular weight of the individual amino acid

    :param molecule: string
    :return: mass of the molecule
    """
    return sum(
        monoisotopic[atom] * int(num or '1')
        for atom, num in refindall(r'([A-Z][a-z]*)(\d*)', molecule)
    )
