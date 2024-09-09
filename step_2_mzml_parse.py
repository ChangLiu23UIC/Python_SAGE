import math
import itertools
import pandas as pd
from b_y_ion_utility import *
from pyteomics import mzml
import lxml


def read_mzml_file(mzml_file):
    """
    This will read the mzml file and return the MS1 and MS2 spectra.
    :param mzml_file:
    :return: spec_ms1_list (MS1 spectra), spec_ms2_list (MS2 spectra)
    """
    with mzml.MzML(mzml_file) as reader:
        spec_ms1_list = {}
        spec_ms2_list = {}  # Initialize once outside the loop to store all MS2 spectra

        for spectrum in reader:
            scan_id = spectrum.get('id')
            ms_level = spectrum.get('ms level')

            # Accessing m/z and intensity arrays
            mz_array = spectrum['m/z array']
            intensity_array = spectrum['intensity array']

            # Extract scan ID (making sure it's formatted correctly)
            scan_id = scan_id.split(' ')[-1].split("=")[1]

            # Create a tuple of (m/z, intensity)
            spec = list(zip(mz_array, intensity_array))

            if ms_level == 2:
                # Access precursor information from the spectrum
                if 'precursorList' in spectrum:
                    precursor = spectrum['precursorList']['precursor'][0]
                    precursor_mz = precursor['selectedIonList']['selectedIon'][0]['selected ion m/z']
                    precursor_scan_id = precursor['spectrumRef'].split(' ')[-1].split("=")[1]

                    # Add precursor scan id and mz to MS2 spectrum list
                    if precursor_scan_id not in spec_ms2_list:
                        spec_ms2_list[precursor_scan_id] = {}

                    spec_ms2_list[precursor_scan_id][precursor_mz] = spec

            elif ms_level == 1:
                # Store MS1 spectrum
                spec_ms1_list[scan_id] = spec

    return spec_ms1_list, spec_ms2_list


if __name__ == '__main__':
    import time

    start_time = time.time()
    spec1, spec2 = read_mzml_file("OV3-DMSO-n3-F9.mzML")
    end_time = time.time()

    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.2f} seconds")