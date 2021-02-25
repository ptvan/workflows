import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 3             
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, facecolor='white')