import glob
import os
import sys

# Gets the Global ClassPath
def getClassPath():
	separator = ":"
	if os.name == "nt":
		separator = ";"
	netBoxHome = os.environ.get("NETBOX_HOME")
	if netBoxHome == None:
		print ("NETBOX_HOME Environment Variable is not set.  Aborting.")
		sys.exit()
	else:
		classpath = ""
		# use os.path.join so that this works across operating systems.
		temp = os.path.join ("%s" % (netBoxHome), "lib", "*.jar")
		jars = glob.glob(temp)
		for jar in jars:
			classpath += "%s%s" % (jar, separator)
		return classpath
