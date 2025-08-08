__version__ = '0.1.0'

from .utils import *
from .parsing import *
from .pipeline import *

"""
SynDesign - End-to-End Design Tool for Saturation Prime Editing Powered by DeepPrime

=========================
Modular pipeline version 

Modules:
- sd_utils:         Common utilities and constants.
- sd_parsing:       COSMIC and ClinVar mutation parsers.
- sd_pegrna:        Input sequence generation for pegRNAs.
- sd_run:           Multiprocessing and model prediction.
- sd_plotting:      Plotting and result formatting.
- sd_main_runner:   Example orchestration script.

Usage:
    from sd_pegrna import make_dp_input_v2
    from sd_run import run_sd_parallel
    from sd_plotting import create_html_output_table

Author: PJM
"""

