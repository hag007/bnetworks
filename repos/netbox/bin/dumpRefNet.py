#!/usr/bin/python

import os
import sys
import env

classpath = env.getClassPath()
os.system ("java -Xmx1192M -classpath %s org.mskcc.netbox.script.DumpReferenceNetwork" % (classpath))
