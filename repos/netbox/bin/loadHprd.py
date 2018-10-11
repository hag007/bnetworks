#!/usr/bin/python

import os
import sys
import env

if len(sys.argv) == 2:
	classpath = env.getClassPath()
	os.system ("java -classpath %s org.mskcc.netbox.script.ImportHprd %s" % (classpath, sys.argv[1]))
else:
	print ("command line usage:  %s <interactions_file>" % (sys.argv[0]))

