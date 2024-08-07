from expasy import expasy_rules
import re


def peptide_cleavage(rule:str, sequence:str, mis_cleavage:int = 2, min_len:int = 7, max_len:int = 50):
    """
    Method for peptide cleavage with the expasy rule selected with regular expression. Default min len is 7 and max len
    is 50 with miscleavage of 2.
    :param rule:
    :param sequence:
    :return:
    """
    exp_rule = expasy_rules[rule]
    cut_end_pos_list = [i.start() for i in re.finditer(exp_rule, sequence)] + [(len(sequence) - 1)]
    cut_start_pos_list = [0] + [i + 1 for i in cut_end_pos_list[:-1]]

    # pieces_list = {"start_pos": cut_start_pos_list, "end_pos": cut_end_pos_list}

    result = [
        sequence[cut_start_pos_list[i]:cut_end_pos_list[j]]
        for i in range(len(cut_start_pos_list))
        for j in range(i, min(i + mis_cleavage + 1, len(cut_start_pos_list)))
        if min_len <= len(sequence[cut_start_pos_list[i]:cut_end_pos_list[j]]) <= max_len
    ]

    return result