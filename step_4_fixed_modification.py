from step_1_fasta_cleavage import *
def apply_fixed_modifications(sequence, modifications):
    """
    Apply fixed modifications to the specified amino acids in the sequence.

    :param sequence: The amino acid sequence as a string.
    :param modifications: A dictionary where keys are amino acids and values are mass changes (positive or negative).
    :return: A list of tuples where each tuple contains the amino acid and its total mass (base mass + modification).
    """


    # Initialize list to store modified masses
    modified_sequence = []

    # Iterate through each amino acid in the sequence
    for aa in sequence:
        # Get the original mass of the amino acid
        base_mass = protein_weight(sequence)

        # Apply modification if it exists; otherwise, no change (modification is 0 by default)
        modified_mass = base_mass + modifications.get(aa, 0)

        # Store the amino acid and its modified mass
        modified_sequence.append((aa, modified_mass))

    return modified_sequence

