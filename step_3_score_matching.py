import numpy as np
import heapq


def insert_with_limit(new_element, limit=5):
    """
    Heap tree modification of the elements
    Call this function in a for loop for the tree generation
    :param new_element:
    :param limit:
    :return:
    """
    heap = []
    if len(heap) < limit:
        heapq.heappush(heap, new_element)
    else:
        if new_element > heap[0]:
            heapq.heapreplace(heap, new_element)
    return heap


def calculate_spectral_angle(theo_arr, maml_arr):
    """
    Calculate the angle similarity between two vectors(arraies)
    :param theo_arr:
    :param maml_arr:
    :return:
    """
    dot_product = np.dot(theo_arr, maml_arr)
    norm_theoretical = np.linalg.norm(theo_arr)
    norm_observed = np.linalg.norm(maml_arr)

    if norm_observed == 0 or norm_theoretical == 0:
        return 0.0

    similarity = dot_product / (norm_theoretical * norm_observed)
    return 1.0 - (2.0 * np.arccos(np.clip(similarity, -1.0, 1.0)) / np.pi)

def integrate_peaks(scores, spectral_angles, settings):
    """

    :param scores:
    :param spectral_angles:
    :param settings:
    :return:
    """
    best_rt = np.argmax(scores)
    best_score = scores[best_rt]

    if best_score == 0.0:
        return None

    threshold = best_score * 0.50
    left = np.searchsorted(scores[:best_rt], threshold, side='right')
    right = np.searchsorted(scores[best_rt:], threshold, side='left') + best_rt

    if settings['integration'] == 'sum':
        peak_area = np.sum(scores[left:right + 1])
    elif settings['integration'] == 'apex':
        peak_area = scores[best_rt]

    return {
        'rt': best_rt,
        'score': best_score,
        'peak_area': peak_area,
        'spectral_angle': spectral_angles[best_rt]
    }


if __name__ == '__main__':
    import time
    from step_1_fasta_cleavage import *

    start_time = time.time()


    # Example usage
    rule = 'trypsin'
    filename = 'human_1_rev.fasta'
    uniprot_sequences = run_fasta_to_class_parallel(rule, filename)

    end_time = time.time()


    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.2f} seconds")