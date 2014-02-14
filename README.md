#speedseq         
-------------------------------

**Current version:** 0.0.1a

Created by Colby Chiang, Ryan Layer, Greg G Faust, Michael R Lindberg, Aaron R Quinlan, and Ira M Hall.

Current support for Linux only

##Summary
--------------
The ``speedseq`` suite is a lightweight, flexible, and open source pipeline that identifies
genomic variation (Structural Variants, INDELs, and Single Nucleotide Variants). 
There are two modes of analysis supported: 

**1.) Identification of variants in a single sample.**

**2.) Comparison of two samples, i.e. tumor and matched normal**

##Constitutive Pipeline Tools (Required)
------------------------------------------

1.) [BWA](http://bio-bwa.sourceforge.net/)

2.) [FREEBAYES](https://github.com/ekg/freebayes)

3.) [GEMINI](http://gemini.readthedocs.org)

4.) [LUMPY](https://github.com/arq5x/lumpy-sv)

- gnu scientific library

5.) [PARALLEL](http://www.gnu.org/software/parallel/)

6.) [SAMBAMBA](https://github.com/lomereiter/sambamba)

7.) [SAMBLASTER](https://github.com/GregoryFaust/samblaster)

8.) [SNPEFF](http://snpeff.sourceforge.net/)

9.) [VCFLIB](https://github.com/ekg/vcflib)


##Installation
----------------

There is an automatic (coming soon) and manual installation process for ``speedseq``.

The following are required for both installations:
- cmake
- g++
- gcc
- git
- make
- python27
- python-devel
- python-yaml
- ncurses-devel
- zlib-devel

A Linux package manager can be used to obtain these by with the command:

~~~~~~~~~~~~~~~~~~~
sudo yum update
sudo yum -y install cmake gcc-c++ gcc git make python27 python-devel python-yaml ncurses-devel zlib-devel 
~~~~~~~~~~~~~~~~~~~
 
or 

~~~~~~~~~~~~~~~~~~~
sudo apt-get update
sudo apt-get install build-essential cmake gpp gcc git make python2.7 python-dev python-yaml ncurses-dev zlib1g-dev 
~~~~~~~~~~~~~~~~~~~

###Automatic installation (coming soon)

``speedseq`` can be installed with the following commands:
~~~~~~~~~~~~~~~~~~
	git clone https://github.com/cc2qe/speedseq
	cd speedseq
	sudo python setup.py install
	sudo cp -r bin/* /usr/local/bin/
~~~~~~~~~~~~~~~~~~

###Manual installation

The following instructions for installation assumes that the required tools are not installed.  
It is recommended that the specified versions of each tool is used for this release of ``speedseq``.  
The use of unspecified versions of any pipeline component is not guaranteed to work. 

``speedseq`` can be installed with the following commands: 
~~~~~~~~~~~~~~~~~~
	git clone https://github.com/cc2qe/speedseq
	cd speedseq
	sudo cp -r bin/* /usr/local/bin/
~~~~~~~~~~~~~~~~~~

####Obtain each of the pipeline tools and install:
	
####1.) BWA

``bwa`` can be installed and used by ``speedseq`` with the following commands: 
~~~~~~~~~~~~~~~~~~
	curl -OL http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.6a.tar.bz2
	tar -xvf bwa-0.7.6a.tar.bz2
	cd bwa-0.7.6a
	make
	sudo cp bwa /usr/local/bin
~~~~~~~~~~~~~~~~~~
-
####2.) FREEBAYES

``freebayes`` can be installed and used by ``speedseq`` with the following commands: 
~~~~~~~~~~~~~~~~~~~
	git clone --recursive git://github.com/ekg/freebayes.git
	cd freebayes
	make
	sudo cp -r bin/* /usr/local/bin/
~~~~~~~~~~~~~~~~~~~
-
####3.) GEMINI

``gemini`` can be installed and used by ``speedseq`` with the following commands: 
~~~~~~~~~~~~~~~~~~~
	wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py	
	#or curl -OL https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py > gemini_install.py
	sudo python2.7 gemini_install.py /usr/local /usr/local/share/gemini
	export PATH=$PATH:/usr/local/gemini/bin
	# it would be wise to add the above line to your ``.bashrc`` or ``.bash_profile``
	gemini update
~~~~~~~~~~~~~~~~~~~
-
####4.) LUMPY

``lumpy-sv`` can be installed and used by ``speedseq`` with the following commands:

- gnu scientific library
~~~~~~~~~~~~~~~~~~~
	curl -OL ftp://ftp.gnu.org/gnu/gsl/gsl-1.9.tar.gz
	tar -xvf gsl-1.9.tar.gz
	cd gsl-1.9
	./configure && sudo make && sudo make install
~~~~~~~~~~~~~~~~~~~

Now that software dependencies have been met, install ``lumpy-sv``:
~~~~~~~~~~~~~~~~~~~
	curl -OL https://github.com/arq5x/lumpy-sv/archive/v0.1.5.tar.gz
	tar -xvf v0.1.5.tar.gz
	cd lumpy-sv-0.1.5
	make 
	sudo cp -r bin/* /usr/local/bin/
	sudo cp -r scripts/* /usr/local/bin/
~~~~~~~~~~~~~~~~~~~
-
####5.) PARALLEL

``parallel`` can be installed and used by ``speedseq`` with the following commands:
~~~~~~~~~~~~~~~~~~~
	curl -OL http://ftp.gnu.org/gnu/parallel/parallel-20100424.tar.bz2
	tar -xvf parallel-20100424.tar.bz2
	cd parallel-20100424
	./configure && sudo make && sudo make install
	sudo cp src/parallel /usr/local/bin/
~~~~~~~~~~~~~~~~~~~
-
####6.) SAMBAMBA

``sambamba`` can be installed and used by ``speedseq`` by: 
~~~~~~~~~~~~~~~~~~
	curl -OL https://github.com/lomereiter/sambamba/releases/download/v0.4.4/sambamba_v0.4.4_centos5.tar.bz2
	tar -xvf sambamba_v0.4.4_centos5.tar.bz2 
	sudo cp sambamba_v0.4.4 /usr/local/bin/
~~~~~~~~~~~~~~~~~~
-
####7.) SAMBLASTER

``samblaster`` can be installed and used by ``speedseq`` with the following commands: 
~~~~~~~~~~~~~~~~~~
	git clone git://github.com/GregoryFaust/samblaster.git
	cd samblaster
	make
	sudo cp samblaster /usr/local/bin/
~~~~~~~~~~~~~~~~~~
-
####8.) SNPEFF

``snpeff`` can be installed and used by ``speedseq`` with the following commands: 
~~~~~~~~~~~~~~~~~~
	wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
	unzip snpEff_latest_core.zip
	cd snpEff
	sudo scp *jar /usr/local/bin/
	sudo scp snpEff.config /usr/local/bin
	sudp scp -r scripts/* /usr/local/bin/
~~~~~~~~~~~~~~~~~~
-
####9.) VCFLIB

``vcflib`` can be installed and used by ``speedseq`` with the following commands: 
~~~~~~~~~~~~~~~~~~
	git clone --recursive  https://github.com/ekg/vcflib
	cd vcflib
	make
	sudo cp bin/* /usr/local/bin/
~~~~~~~~~~~~~~~~~~
-
**For alternative installations and release issues for any of the above tools please consult the website/creator.**

##Example Usage
----------------------

Use ``speedseq aln`` to align and dedupe the dataset.
~~~~~~~~~~~~~~~~~~
	speedseq aln -o NA12878 human_g1k_v37.fasta NA12878.fq.gz
~~~~~~~~~~~~~~~~~~

Use ``speedseq var`` to call SNPs and indels on a single sample.
~~~~~~~~~~~~~~~~~~
	speedseq var -o NA12878 -w annotations/ceph18.b37.include.2014-01-15.bed human_g1k_v37.fasta NA12878.bam
~~~~~~~~~~~~~~~~~~

Use ``speedseq lumpy`` to call structural variants.
~~~~~~~~~~~~~~~~~~
	speedseq lumpy -o NA12878 -x annotations/ceph18.b37.exclude.2014-01-15.bed -B NA12878.bam -D NA12878.discordants.bam -S NA12878.splitters.bam
~~~~~~~~~~~~~~~~~~

Use ``speedseq somatic`` to call SNPs and indels on a tumor/normal pair.
~~~~~~~~~~~~~~~~~~
	speedseq somatic -o TCGA-B6-A0I6 -w annotations/ceph18.b37.include.2014-01-15.bed human_g1k_v37.fasta TCGA-B6-A0I6.normal.bam TCGA-B6-A0I6.tumor.bam
~~~~~~~~~~~~~~~~~~




