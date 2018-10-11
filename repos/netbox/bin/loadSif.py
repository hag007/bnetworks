#!/usr/bin/python

import os
import sys
import env

if len(sys.argv) == 3:
	classpath = env.getClassPath()
	os.system ("java -classpath %s org.mskcc.netbox.script.ImportSif %s %s" % (classpath, sys.argv[1], sys.argv[2]))
else:
	print ("command line usage:  %s <interactions_file.sif> <data_source_name>" % (sys.argv[0]))

