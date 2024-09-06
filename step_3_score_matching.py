import numpy as np


def calculate_spectral_angle(theoretical_dist, mzml_dict):
    dot_product = np.dot(theoretical_dist, mzml_dict)
    norm_theoretical = np.linalg.norm(theoretical_dist)
    norm_observed = np.linalg.norm(mzml_dict)

    if norm_observed == 0 or norm_theoretical == 0:
        return 0.0

    similarity = dot_product / (norm_theoretical * norm_observed)
    return 1.0 - (2.0 * np.arccos(np.clip(similarity, -1.0, 1.0)) / np.pi)

def integrate_peaks(scores, spectral_angles, settings):
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