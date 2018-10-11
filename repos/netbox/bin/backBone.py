#!/usr/bin/python

import os
import sys
import env

classpath = env.getClassPath()
cmd = ""
for arg in sys.argv[1:]:
    cmd += arg + " "
os.system ("java -Xmx1192M -classpath %s org.mskcc.netbox.script.BackBone %s" % (classpath, cmd))
