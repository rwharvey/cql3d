#!/usr/bin/python

import sys, os

for fn in sys.argv[2:]:
    print fn
    os.popen('dop '+fn+' tmpmpi.f '+sys.argv[1])
    os.popen('mv tmpmpi.f '+fn)
    
