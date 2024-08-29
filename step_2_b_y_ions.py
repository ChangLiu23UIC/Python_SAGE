import math
import itertools
import pandas as pd
from b_y_ion_utility import *
from pyteomics import mzml
import lxml

def read_mzml_file(mzml_file):
    with mzml.MzML(mzml_file) as reader:
        for spectrum in reader:
            scan_id = spectrum.get('id')
            ms_level = spectrum.get('ms level')

            # Accessing m/z and intensity arrays
            mz_array = spectrum['m/z array']
            intensity_array = spectrum['intensity array']

            # # Print some information about the spectrum
            # print(f"Spectrum ID: {scan_id}")
            # print(f"MS Level: {ms_level}")
            # print(f"Retention Time: {retention_time} minutes")
            # print(f"Number of peaks: {len(mz_array)}")

            for mz, intensity in zip(mz_array, intensity_array):
                print(f"m/z: {mz}, Intensity: {intensity}")
                print(len(mz_array))

if __name__ == '__main__':
    import time

    start_time = time.time()
    read_mzml_file("OV3-DMSO-n3-F9.mzML")
    end_time = time.time()

    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.2f} seconds")