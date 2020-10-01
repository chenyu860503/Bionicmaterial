# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 21:26:05 2020

@author: USER
"""
import sys
import numpy as np
import glob
import os 
from pathlib import Path
filepath = 'C:/Users/USER/Desktop/BionicMaterial/script/dump/'
dumpfile = 'dump.Si.*'
bonddump = 'bond.dump.*'

atoms_files = sorted(glob.glob(filepath+dumpfile))
atoms_files = print(glob.glob(filepath+dumpfile))