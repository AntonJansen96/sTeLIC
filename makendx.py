#!/usr/bin/env python3

from science.utility import createIndexFile
createIndexFile('charmm-gui-closed.pdb', 'chap.ndx', groups=['resid 225 and name CA'])
