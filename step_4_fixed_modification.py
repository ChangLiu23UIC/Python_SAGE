from step_1_fasta_cleavage import *
def apply_fixed_modifications(sequence, modifications_dict):
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
        base_mass = protein_weight(sequence)
        modified_mass = base_mass + modifications_dict.get(aa)
        modified_sequence.append((aa, modified_mass))

    return modified_sequence


def apply_variable_modifications(fixed_sequence, variable_modifications):
    """
    Apply variable modifications to specific amino acids and generate all possible combinations.

    :param fixed_sequence:
    :param variable_modifications:
    :return:
    """
    variable_positions = []

    for i, (aa, mass) in enumerate(fixed_sequence):
        if aa in variable_modifications:
            variable_positions.append((i, variable_modifications[aa]))

    if variable_positions:
        all_combinations = []
        modification_options = []
        for i, mods in variable_positions:
            mod_options = [(0, fixed_sequence[i][1])]
            for mod in mods:
                mod_options.append((mod, fixed_sequence[i][1] + mod))
            modification_options.append(mod_options)

        # Generate all combinations of modified and unmodified amino acids
        for combination in itertools.product(*modification_options):
            modified_sequence = fixed_sequence.copy()
            for (i, _), (mod, new_mass) in zip(variable_positions, combination):
                modified_sequence[i] = (modified_sequence[i][0], new_mass)
            all_combinations.append(modified_sequence)

        return all_combinations
    else:
        # If no variable modifications, return the fixed sequence as is
        return [fixed_sequence]
