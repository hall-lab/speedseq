#!/usr/bin/env python
"""Installer for speedseq: lightweight, flexible, and open source pipeline that identifies genomic variation 

https://github.com/cc2qe/speedseq

Handles installation of:
 
- Required third party software

Requires: 
Linux packages: cmake gcc-c++ gcc git make python27 python-devel python-yaml ncurses-devel zlib-devel 
Bioinformatics software: BWA, FREEBAYES, GEMINI (bedtools, samtools, pybedtools), LUMPY, PARALLEL, SAMBAMBA, SAMBLASTER, SNPEFF, VCFLIB

Run speedseq_install.py -h for usage.
"""
import argparse
import os
import shutil
import subprocess
import sys
import shlex

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
		elif (method == "python2.7"):
			#This only applies to gemini, uses pip install
			self.installcmd =  "python2.7 gemini_install.py /usr/local " + self.installdir + " && gemini update"
		else:
			self.installcmd = 'echo "snpEff is java"'
		if not self.quiet:
			print "\nInstalling " + self.name + "...\n"
		subprocess.call(self.installcmd, shell=True)

	def cp_bin(self, source, target):
		if os.path.isfile(source):
			self.copycmd = "sudo cp "  + os.getcwd() + "/" + source + " " + target
		elif os.path.isdir(source):
				self.copycmd = "sudo cp -r "  + os.getcwd() + "/" + source + "/* " + target
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

def which(prgm):
	for path in os.environ["PATH"].split(":"):
		if os.path.exists(path + "/" + prgm):
			return path + "/" + prgm
	return None

def main(args):
	print "Changing directory to " + args.targetdir + " and installing speedseq...\n"
	if os.path.isdir(args.targetdir):
		os.chdir(args.targetdir)
	else:
		raise OSError("cd: " + args.targetdir + ": No such file or directory")
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
		url = "http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.6a.tar.bz2"
		bwa.download("curl", url)
		bwa.unpack("tar")
		bwa.install("make", "bwa-0.7.6a")
		bwa.cp_bin("bwa-0.7.6a/bwa", args.targetbin)
		
	#freebayes install
	freebayes = INSTALLER("freebayes", args.quiet)
	freebayes.check_install("freebayes")
	if (freebayes.isInstalled):
		freebayes.get_update()
	if (freebayes.notInstalled or freebayes.update):
		url="git://github.com/ekg/freebayes"
		freebayes.download("git", url)
		freebayes.install("make", "freebayes")
		freebayes.cp_bin("freebayes/bin", args.targetbin)
    	#gemini install
	gemini = INSTALLER("gemini", args.quiet)
	gemini.check_install("gemini")
	if (gemini.isInstalled):
		gemini.get_update()
	if (gemini.notInstalled or gemini.update):
		url = "https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py"
		gemini.download("wget", url)
		gemini.install("python2.7", "/usr/local/share/gemini")
		gemini.cp_bin("/usr/local/gemini/bin", args.targetbin)
	#lumpy install
	lumpy = INSTALLER("lumpy", args.quiet)
	lumpy.check_install("lumpy")
	if (lumpy.isInstalled):
		lumpy.get_update()
	if (lumpy.notInstalled or lumpy.update):
		url="https://github.com/arq5x/lumpy-sv/archive/0.2.1.tar.gz"
		lumpy.download("curl", url)
		lumpy.unpack("tar") 
		lumpy.install("make", "lumpy-sv-0.2.1")
		lumpy.cp_bin("lumpy-sv-0.2.1/bin", args.targetbin)
		lumpy.cp_bin("lumpy-sv-0.2.1/scripts", args.targetbin)
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
	sambamba = INSTALLER("sambamba_v0.4.6", args.quiet)
	sambamba.check_install("sambamba_v0.4.6")
	if (sambamba.isInstalled):
		sambamba.get_update()
	if (sambamba.notInstalled or sambamba.update):
		url = "https://github.com/lomereiter/sambamba/releases/download/v0.4.6-beta/sambamba_v0.4.6-beta_centos5-x86_64.tar.bz2"
		sambamba.download("curl", url)
		sambamba.unpack("tar")
		sambamba.cp_bin("sambamba_v0.4.6", args.targetbin)
	#samblaster install	
	samblaster = INSTALLER("samblaster", args.quiet)
	samblaster.check_install("samblaster")
	if (samblaster.isInstalled):
		samblaster.get_update()
	if (samblaster.notInstalled or samblaster.update):
		url = "https://github.com/GregoryFaust/samblaster/archive/0.1.14.tar.gz"
		samblaster.download("curl", url)
		samblaster.unpack("tar")
		samblaster.install("make", "samblaster-0.1.14")
		samblaster.cp_bin("samblaster-0.1.14/samblaster", args.targetbin)
	#snpeff install
	snpeff = INSTALLER("snpeff", args.quiet)
	snpeff.check_install("snpeff")
	if (snpeff.isInstalled):
		snpeff.get_update()
	if (snpeff.notInstalled or snpeff.update):
		url = "http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip"
		snpeff.download("wget", url)
		snpeff.unpack("unzip")
		snpeff.install(None, "snpEff")
		snpeff.cp_bin("snpEff/snpEff.config", args.targetbin)
		snpeff.cp_bin("snpEff/snpEff.jar", args.targetbin)
		snpeff.cp_bin("snpEff/scripts", args.targetbin)
	#vcflib install
	vcflib = INSTALLER("vcflib", args.quiet)
	if (vcflib.isInstalled):
		vcflib.get_update()
	if (vcflib.notInstalled or vcflib.update):
		url = "https://github.com/ekg/vcflib"
		vcflib.download("git", url)
		vcflib.install("make", "vcflib")
		vcflib.cp_bin("vcflib/tabixpp/bgzip", args.targetbin)
		vcflib.cp_bin("vcflib/tabixpp/tabix", args.targetbin)
		vcflib.cp_bin("vcflib/bin", args.targetbin)
	
	print "Checking installations...\n"
	bwa.check_install("bwa")
	freebayes.check_install("freebayes")
	gemini.check_install("gemini")
	lumpy.check_install("lumpy")
	parallel.check_install("parallel")
	sambamba.check_install("sambamba_v0.4.6")
	samblaster.check_install("samblaster")
	snpeff.check_install("snpeff")
	vcflib.check_install("bgzip")


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
		parser.add_argument("targetdir", help="Directory to install 3rd party software tools",
						nargs='?', type=os.path.abspath, default=os.getcwd())
		parser.add_argument("targetbin", help="Directory to install the binaries into",
						nargs='?', type=os.path.abspath, default="/usr/local/bin/")
		parser.add_argument("--quiet", '-q', help="Determines the verbosity of installation",
						default=False, action='store_true')
		if len(sys.argv) == 0:
			parser.print_help()
		else:
			main(parser.parse_args())
