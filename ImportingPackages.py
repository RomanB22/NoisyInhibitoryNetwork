"""File for importing packages and configure the BRIAN2 simulator to use in the HPC cluster

"""
import matplotlib
matplotlib.use('Agg') 
import sys
import brian2 as b2
import numpy as np
from scipy.optimize import root_scalar, fmin, curve_fit, minimize
import time
import os, glob
import matplotlib.pyplot as plt
import scipy.stats
import itertools as it
import re # Module for treating with regular expressions
import seaborn as sns
from scipy import spatial

TigerFish = True
if TigerFish: 	# For Tigerfish use only
	cache_dir = os.path.expanduser('~/.cython/brian-pid-{}'.format(os.getpid()))
	b2.prefs.codegen.runtime.cython.cache_dir = cache_dir
	b2.prefs.codegen.runtime.cython.multiprocess_safe = False
