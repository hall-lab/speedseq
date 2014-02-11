#speedseq         
-------------------------------

**Current version:** 0.0.1a

Created by Colby Chiang, Ryan Layer, Greg G Faust, Michael R Lindberg, Aaron R Quinlan, and Ira M Hall.

Current support is for Linux support only

##Summary
--------------
The ``speedseq`` suite is a lightweight, flexible, and open source pipeline that identifies
genomic variation (Structural Variants, INDELs, and Single Nucleotide Variants). 
There are two modes of analysis supported: 

1.) Identification of variants in a single sample.

2.) Comparison of two samples, i.e. tumor and matched normal

##Constitutive Pipeline Tools (Required)
------------------------------------------

1.) [BEDTOOLS](https://github.com/arq5x/bedtools)

2.) [BWA](https://github.com/lh3/bwa)

3.) [FREEBAYES](https://github.com/ekg/freebayes)

4.) [LUMPY](https://github.com/arq5x/lumpy-sv)

5.) [SAMBAMBA](https://github.com/lomereiter/sambamba)

6.) [TABIX](https://github.com/samtools/tabix)

7.) [SNPEFF](https://github.com/CBMi-BiG/snpEff)

8.) [VCFLIB](https://github.com/ekg/vcflib)

9.) [GEMINI](https://github.com/arq5x/gemini)

##Installation
----------------

``speedq`` can be installed with the following commands: 
~~~~~~~~~~~~~~~~~~
	git clone https://github.com/cc2qe/speedseq
	cd speedseq
	make
	sudo scp -r bin/* /usr/local/bin/
~~~~~~~~~~~~~~~~~~

python setup.py or Makefile coming soon

###Manual installation

The following installation process assumes that none of the above required tools are installed.  
It is recommended that recent versions of tools are used in this release of ``speedseq`` meaning that
the use of any other versions cannot be guaranteed to work. 

**gpp** and **cmake** are required for compilation


####1.) BEDTOOLS
	
``bedtools`` can be installed and used by ``speedseq`` with the following commands: 
~~~~~~~~~~~~~~~~~~
	git clone https://github.com/arq5x/bedtools
	cd bedtools
	make
	sudo scp -r bin/* /usr/local/bin/
~~~~~~~~~~~~~~~~~~
	
For alternative installations, release issues, and unmentioned dependencies, please consult the tool website/creator.

####2.) BWA

``bwa`` can be installed and used by ``speedseq`` with the following commands: 
~~~~~~~~~~~~~~~~~~
	git clone https://github.com/lh3/bwa
	cd bwa
	sudo scp bwa /usr/local/bin
~~~~~~~~~~~~~~~~~~
	
For alternative installations, release issues, and unmentioned dependencies, please consult the tool website/creator.

####3.) FREEBAYES

``freebayes`` can be installed and used by ``speedseq`` with the following commands: 
~~~~~~~~~~~~~~~~~~~
	git clone https://github.com/ekg/freebayes
	cd freebayes
	make
	sudo scp -r bin/* /usr/local/bin/
~~~~~~~~~~~~~~~~~~~
	
For alternative installations, release issues, and unmentioned dependencies, please consult the tool website/creator.

####4.) LUMPY

``lumpy-sv`` can be installed and used by ``speedseq`` with the following commands:
~~~~~~~~~~~~~~~~~~~
	git clone https://github.com/arq5x/lumpy-sv
	cd lumpy-sv
	make 
	sudo scp -r bin/* /usr/local/bin/
	sudo scp -r scripts/* /usr/local/bin/
~~~~~~~~~~~~~~~~~~~
	
For alternative installations, release issues, and unmentioned dependencies, please consult the tool website/creator.


####5.) SAMBAMBA

``sambamba`` can be installed and used by ``speedseq`` by: 

* Go to https://github.com/lomereiter/sambamba/releases

* Download sambamba_v0.4.4_centos5.tar.bz2
~~~~~~~~~~~~~~~~~~
	tar -xvf sambamba_v0.4.4_centos5.tar.bz2
	sudo scp sambamba_v0.4.4 /usr/local/bin/
~~~~~~~~~~~~~~~~~~

For alternative installations, release issues, and unmentioned dependencies, please consult the tool website/creator.

####6.) TABIX

``tabix`` can be installed and used by ``speedseq`` with the following commands: 
~~~~~~~~~~~~~~~~~~
	git clone --recursive  git://github.com/samtools/tabix
	cd tabix
	make
	sudo scp tabix /usr/local/bin/
	sudo scp bgzip /usr/local/bin/
~~~~~~~~~~~~~~~~~~
	
For alternative installations, release issues, and unmentioned dependencies, please consult the tool website/creator.


####7.) SNPEFF

``snpeff`` can be installed and used by ``speedseq`` with the following commands: 
~~~~~~~~~~~~~~~~~~
	wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip/download
	unzip snpEff_latest_core.zip
	cd snpEff
	sudo scp *jar /usr/local/bin/
	sudo scp snpEff.config /usr/local/bin
	sudp scp -r scripts/* /usr/local/bin/
~~~~~~~~~~~~~~~~~~

For alternative installations, release issues, and unmentioned dependencies, please consult the tool website/creator.


####8.) VCFLIB

``vcflib`` can be installed and used by ``speedseq`` with the following commands: 
~~~~~~~~~~~~~~~~~~
	git clone --recursive  https://github.com/ekg/vcflib
	cd vcflib
	make
	sudo scp
~~~~~~~~~~~~~~~~~~

For alternative installations, release issues, and unmentioned dependencies, please consult the tool website/creator.


####9.) GEMINI

``gemini`` can be installed and used by ``speedseq`` with the following commands: 
~~~~~~~~~~~~~~~~~~
	sudo yum -y install python27 git gcc gcc-c++ zlib-devel
	git clone https://github.com/arq5x/gemini
	cd gemini
	sudo python setup.py install
	sudo python gemini/install-data.py /usr/local/share/
~~~~~~~~~~~~~~~~~~

For alternative installations, release issues, and unmentioned dependencies, please consult the tool website/creator.


##Example Use
----------------------

Use ``speedseq aln`` to align and dedupe the dataset.
~~~~~~~~~~~~~~~~~~
	speedseq aln -o NA12878 -t 24 -T temp_NA12878 human_g1k_v37.fasta NA12878.fq.gz
~~~~~~~~~~~~~~~~~~


Use ``speedseq somatic`` to call SNPs and Indels.
~~~~~~~~~~~~~~~~~~
	speedseq somatic -o NA12878.som -w annotations/ceph18.b37.include.2014-01-15.bed -t 24 human_g1k_v37.fasta NA12878.bam &> NA12878.som.log

	tabix -p vcf NA12878.som.annot.vcf.gz
	
	# somatic variant calling with min_alt_fraction 0.01 and snpEff annotation.
	speedseq somatic -A -F 0.01 -o NA12878.som.F01 -w annotations/ceph18.b37.include.2014-01-15.bed -t 24 /shared/genomes/b37/full/human_g1k_v37.fasta NA12878.bam &> NA12878.som.F01.log
	
	# make a vcf with min qual 1
	zcat NA12878.som.F01.vcf.gz | awk '$1~"^#" || $6>=1' | bgzip -c > NA12878.som.F01.Q1.vcf.gz
~~~~~~~~~~~~~~~~~~

Use ``speedseq lumpy`` to call structural variants.
~~~~~~~~~~~~~~~~~~
	speedseq lumpy -o NA12878.lumpy -x annotations/ceph18.b37.exclude.2014-01-15.bed -T temp_NA12878 -B NA12878.bam -D NA12878.discordants.bam -S NA12878.splitters.bam &> NA12878.lumpy.log
~~~~~~~~~~~~~~~~~~








