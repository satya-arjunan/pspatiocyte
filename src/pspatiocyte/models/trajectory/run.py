#!/usr/bin/env python3

import subprocess
import csv

subprocess.run(['mpirun', '-np', '8', 'trajectory'], stdout=subprocess.PIPE)
subprocess.run(['python3', '../../scripts/gather_coordinates.py'])
