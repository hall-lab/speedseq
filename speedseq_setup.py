#!/usr/bin/env python
"""Installer for speedseq: lightweight, flexible, and open source pipeline that identifies genomic variation 

https://github.com/cc2qe/speedseq

Handles installation of:
 
- Required third party software

Requires: 
Linux packages: cmake gcc-c++ gcc git make python27 python-devel python-yaml ncurses-devel zlib-devel 

Run speedseq_install.py -h for usage.
"""
import argparse
import os
import shutil
import subprocess
import sys
import shlex
import tempfile

class PACMAN(object):
	"""
	A class to deal with finding/using the linux specific package manager
	"""
	def __init__(self, quiet):
		self.isYum     = False
		self.isApt_get = False
		self.quiet     = quiet
	def check_installer(self):
		if not self.quiet:
			print "\nLooking for Linux package installer...\n"
		if (which("yum")):
			if not self.quiet:
				print "\nYum found...\n"
			self.isYum = True
		if (which("apt-get")):
			if not self.quiet:
				print "\nApt-get found...\n"
			self.isApt_get = True
		if (not self.isYum and not self.isApt_get):
			sys.stderr('Linux package installer (yum/apt-get) cannot be found, make sure it is located in your $PATH')
	def install(self):
		if (self.isYum):
			#centos/redhat - yum
			subprocess.call(shlex.split("sudo yum update"))
			packages = ["cmake", "gcc-c++", "gcc", "git", "make", "python27", "python-devel", \
			"python-yaml", "ncurses-devel", "zlib-devel"]
			for i in range(len(packages)):
				subprocess.call(shlex.split("sudo yum install " + packages[i]))
		elif (self.isApt_get):
			#ubuntu - apt-get
			subprocess.call(shlex.split("sudo apt-get update"))
			subprocess.call(shlex.split("sudo apt-get install build-essential cmake gpp gcc git " + \
			"make python2.7 python-dev python-yaml ncurses-dev zlib1g-dev"))

class INSTALLER(object):
	"""
	A class to deal with installing each of the software pieces in speedseq
	"""
	def __init__(self, software, quiet):
		self.name         = software
		self.quiet        = quiet
		self.isInstalled  = False
		self.notInstalled = True
		self.update       = False

	def download(self, method, url):
		self.filename = url.split('/')[-1]
		if (method == "curl"):
			self.dlcmd = "curl -OL " + url 
		elif method == "wget":
			self.dlcmd = "wget " + url 
		elif method == "git":
			self.dlcmd = "git clone --recursive " + url
		if not self.quiet:
			print "\nDownloading " + self.name + "...\n"
		subprocess.call(shlex.split(self.dlcmd))

	def unpack(self, method):
		if (method == "tar"):
			self.unpackcmd = "tar -xvf " + self.filename
			#get the directory things were unpacked into, tar has 2 extentions, use basename
		elif (method == "unzip"):
			self.unpackcmd = "unzip " + self.filename
			#get the directory things were unpacked into, zip has 1 extension, use basename
		if not self.quiet:
			print "\nUnpacking " + self.name + "...\n"
		subprocess.call(shlex.split(self.unpackcmd))

	def install(self, method, installdir):
		#This method is messy and very hardcoded
		self.installdir = installdir
		if (method == "make"):
			self.installcmd = "make -C " + self.installdir
		elif (method == "confmake"):
			self.installcmd = "./configure && sudo make && sudo make install"
		elif (method == "perl"):
			#This only applies to VEP
			self.installcmd = "perl " + self.installdir + "/INSTALL.pl -a ac -s homo_sapiens -y GRCh37"
		elif (method == "python2.7"):
			#This only applies to gemini, uses pip install
			self.installcmd = "python2.7 gemini_install.py /usr/local " + self.installdir + " && gemini update"
		if not self.quiet:
			print "\nInstalling " + self.name + "...\n"
		subprocess.call(self.installcmd, shell=True)

	def cp_bin(self, source, target):
		self.copycmd = "sudo cp -r " + os.getcwd() + "/" + source + " " + target
		if not self.quiet:
			print "\nCopying " + source + " from " + self.name + " to target bin ...\n"
		subprocess.call(self.copycmd, shell=True)

	def check_install(self, exe):
		if which(exe):
			print self.name + " is installed...\n"
			self.notInstalled = False
			self.isInstalled  = True 
		else:
			print self.name + " is not installed...\n"
			self.notInstalled = True
			self.isInstalled  = False

	def get_update(self):
		needInput = True
		while needInput:
			s = raw_input(self.name + " was found.\nIt is recommended to install/update because speedseq may not" + \
				" work with the installed version of " + self.name  + ".\nDo you want to install/update? [y/N]\n")
			s = s.lower()
			if ((s == "n") or (s == "no")):
				print "\nNot installing/updating and could lead to potential problems in speedseq" + self.name + "...\nContinuing anyway...\n"
				self.update = False
				needInput = False
			elif ((s == "y") or (s == "yes")):
				print "\nInstalling/updating " + self.name + "\n"
				self.update = True
				needInput = False
			else:
				print "\nUnrecognized input, please input yes or no [y/N]"
				needInput = True

def make():
	subprocess.call('make -k all', shell=True)

def which(prgm):
	for path in os.environ["PATH"].split(":"):
		if os.path.exists(path + "/" + prgm):
			return path + "/" + prgm
	return None

def main(args):
	# Make temp directory if it doesn't exist
        try:
		os.stat(args.tempdir)
		keep = True
        except:
		args.tempdir = tempfile.mkdtemp(dir="./")
		keep = False
	print "Changing directory to " + args.tempdir + " and installing speedseq...\n"
	if os.path.isdir(args.tempdir):
		prevdir = os.getcwd()
		os.chdir(args.tempdir)
	else:
		raise OSError("cd: " + args.tempdir + ": No such file or directory")

	#Linux install
	packageManager = PACMAN(args.quiet)
	packageManager.check_installer()
	packageManager.install()
	check_dependencies()

	#bwa install
	bwa = INSTALLER("bwa", args.quiet)
	bwa.check_install("bwa")
	if (bwa.isInstalled):
		bwa.get_update()
	if (bwa.notInstalled or bwa.update):
		url = "http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.8.tar.bz2"
		bwa.download("curl", url)
		bwa.unpack("tar")
		bwa.install("make", "bwa-0.7.8")
		bwa.cp_bin("bwa-0.7.8/bwa", args.targetbin)

	#parallel install
	parallel = INSTALLER("parallel", args.quiet)
	parallel.check_install("parallel")
	if (parallel.isInstalled):
		parallel.get_update()
	if (parallel.notInstalled or parallel.update):
		url = "http://ftp.gnu.org/gnu/parallel/parallel-20100424.tar.bz2"
		parallel.download("curl", url)
		parallel.unpack("tar")
		parallel.install("confmake", "parallel-20100424")
		parallel.cp_bin("parallel-20100424/src/parallel", args.targetbin)

	#sambamba install
	sambamba = INSTALLER("sambamba", args.quiet)
	sambamba.check_install("sambamba")
	if (sambamba.isInstalled):
		sambamba.get_update()
	if (sambamba.notInstalled or sambamba.update):
		url = "https://github.com/lomereiter/sambamba/releases/download/v0.4.7/sambamba_v0.4.7_centos5.tar.bz2"
		sambamba.download("curl", url)
		sambamba.unpack("tar")
		sambamba.cp_bin("sambamba_v0.4.7", args.targetbin + "/sambamba")

	#VEP install
	vep = INSTALLER("vep", args.quiet)
	vep.check_install("variant_effect_predictor.pl")
	if (vep.isInstalled):
		vep.get_update()
	if (vep.notInstalled or vep.update):
		url = "https://github.com/Ensembl/ensembl-tools/archive/release/76.zip"
		vep.download("curl", url)
		vep.unpack("unzip")
		vep.install("perl", "ensembl-tools-release-76/scripts/variant_effect_predictor")
		vep.cp_bin("ensembl-tools-release-76/scripts/variant_effect_predictor/variant_effect_predictor.pl", args.targetbin)
		vep.cp_bin("Bio", args.targetbin + "/Bio")
	
	#BEDtools install
	bedtools = INSTALLER("bedtools", args.quiet)
	bedtools.check_install("bedtools")
	if (bedtools.isInstalled):
		bedtools.get_update()
	if (bedtools.notInstalled or bedtools.update):
		url = "https://github.com/arq5x/bedtools2/releases/download/v2.20.1/bedtools-2.20.1.tar.gz"
		bedtools.download("curl", url)
		bedtools.unpack("tar")
		bedtools.install("make", "bedtools2-2.20.1")
		bedtools.cp_bin("bedtools2-2.20.1/bin/*", args.targetbin)

    	# #gemini install
	# gemini = INSTALLER("gemini", args.quiet)
	# gemini.check_install("gemini")
	# if (gemini.isInstalled):
	# 	gemini.get_update()
	# if (gemini.notInstalled or gemini.update):
	# 	url = "https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py"
	# 	gemini.download("wget", url)
	# 	gemini.install("python2.7", "/usr/local/share/gemini")
	# 	gemini.cp_bin("/usr/local/gemini/bin", args.targetbin)
	
	print "Checking installations...\n"
	bwa.check_install("bwa")
	parallel.check_install("parallel")
	sambamba.check_install("sambamba")
	vep.check_install("variant_effect_predictor.pl")
	# gemini.check_install("gemini")

	print "Cleaning up...\n"
	if not keep:
		os.chdir(prevdir)
		shutil.rmtree(args.tempdir)

	print "Compiling embedded software...\n"
	make()

def check_dependencies():
		"""Ensure required tools for installation are present.
		"""
		print "Checking required dependencies..."
		for cmd, url in [("git", "http://git-scm.com/"),
			("wget", "http://www.gnu.org/software/wget/"),
			("curl", "http://curl.haxx.se/"),
			("python2.7", "http://www.python.org/"),
			("make", "https://www.gnu.org/software/make/"),
			("cmake", "http://www.cmake.org/")]:
			try:
				retcode = subprocess.call([cmd, "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			except OSError:
				retcode = 127
			if retcode == 127:
				raise OSError("speedseq requires %s for installation (%s)" % (cmd, url))
			else:
				print " %s found" % cmd

if __name__ == "__main__":
		parser = argparse.ArgumentParser(description="Automated installer for speedseq.")
		parser.add_argument("--tempdir", help="Temp directory to install 3rd party software tools [./temp]",
						nargs='?', type=os.path.abspath, default=None)
		parser.add_argument("--targetbin", help="Directory to install the binaries [./bin]",
						nargs='?', type=os.path.abspath, default="./bin")
		parser.add_argument("--quiet", '-q', help="Determines the verbosity of installation",
						default=False, action='store_true')
		if len(sys.argv) == 0:
			parser.print_help()
		else:
			main(parser.parse_args())
