# batchrun.py --- 
# 
# Filename: batchrun.py
# Description: 
# Author: Subhasis Ray
# Created: Mon May 26 15:52:10 2025 (+0530)
# 

# Code:
"""Simple script to run all python files starting with 'test_' inside
'core' folder"""

import os
import subprocess

ddir = os.path.join(os.path.dirname(__file__), 'core')

if __name__ == '__main__':
    for fname in os.listdir(ddir):
        if fname.startswith('test_') and fname.endswith('.py'):
            fname = os.path.join(ddir, fname)
            print('#' * 70, '\n#\tStarting', fname)
            print('=' * 70)
            subprocess.run(['python', fname], check=True)

# 
# batchrun.py ends here
