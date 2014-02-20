#speedseq         

**Current version:** 0.0.1a

Created by Colby Chiang, Ryan Layer, Greg G Faust, Michael R Lindberg, Aaron R Quinlan, and Ira M Hall.

Current support for Linux only

##Table of Contents

1. [Summary](#summary)
2. [Constitutive Pipeline Tools](#constitutive-pipeline-tools-required)
3. [Installation](#installation)
  * [Automatic Installation](#automatic-installation-coming-soon)
  * [Manual Installation](#manual-installation)
4. [Usage](#usage)
5. [Annotations](#annotations)
6. [Example Workflows](#example-workflows)


##Summary

The `speedseq` suite is a lightweight, flexible, and open source pipeline that identifies
genomic variation (single nucleotide variants (SNVs), indels, and structural variants (SVs)).

##Constitutive Pipeline Tools (Required)

1. [BWA](http://bio-bwa.sourceforge.net/)
2. [FREEBAYES](https://github.com/ekg/freebayes)
3. [GEMINI](http://gemini.readthedocs.org)
4. [LUMPY](https://github.com/arq5x/lumpy-sv)
  * gnu scientific library
5. [PARALLEL](http://www.gnu.org/software/parallel/)
6. [SAMBAMBA](https://github.com/lomereiter/sambamba)
7. [SAMBLASTER](https://github.com/GregoryFaust/samblaster)
8. [SNPEFF](http://snpeff.sourceforge.net/)
9. [VCFLIB](https://github.com/ekg/vcflib)

##Installation

There is an automatic (coming soon) and manual installation process for `speedseq`.

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

`speedseq` can be installed with the following commands:
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

`speedseq` can be installed with the following commands: 
~~~~~~~~~~~~~~~~~~
	git clone https://github.com/cc2qe/speedseq
	cd speedseq
	sudo cp -r bin/* /usr/local/bin/
~~~~~~~~~~~~~~~~~~

####Obtain each of the pipeline tools and install:
	
####1) BWA

`bwa` can be installed and used by ``speedseq`` with the following commands: 
~~~~~~~~~~~~~~~~~~
	curl -OL http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.6a.tar.bz2
	tar -xvf bwa-0.7.6a.tar.bz2
	cd bwa-0.7.6a
	make
	sudo cp bwa /usr/local/bin
~~~~~~~~~~~~~~~~~~
-
####2) FREEBAYES

``freebayes`` can be installed and used by ``speedseq`` with the following commands: 
~~~~~~~~~~~~~~~~~~~
	git clone --recursive git://github.com/ekg/freebayes.git
	cd freebayes
	make
	sudo cp -r bin/* /usr/local/bin/
~~~~~~~~~~~~~~~~~~~
-
####3) GEMINI

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
####4) LUMPY

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
####5) PARALLEL

``parallel`` can be installed and used by ``speedseq`` with the following commands:
~~~~~~~~~~~~~~~~~~~
	curl -OL http://ftp.gnu.org/gnu/parallel/parallel-20100424.tar.bz2
	tar -xvf parallel-20100424.tar.bz2
	cd parallel-20100424
	./configure && sudo make && sudo make install
	sudo cp src/parallel /usr/local/bin/
~~~~~~~~~~~~~~~~~~~
-
####6) SAMBAMBA

``sambamba`` can be installed and used by ``speedseq`` by: 
~~~~~~~~~~~~~~~~~~
	curl -OL https://github.com/lomereiter/sambamba/releases/download/v0.4.4/sambamba_v0.4.4_centos5.tar.bz2
	tar -xvf sambamba_v0.4.4_centos5.tar.bz2 
	sudo cp sambamba_v0.4.4 /usr/local/bin/
~~~~~~~~~~~~~~~~~~
-
####7) SAMBLASTER

``samblaster`` can be installed and used by ``speedseq`` with the following commands: 
~~~~~~~~~~~~~~~~~~
	git clone git://github.com/GregoryFaust/samblaster.git
	cd samblaster
	make
	sudo cp samblaster /usr/local/bin/
~~~~~~~~~~~~~~~~~~
-
####8) SNPEFF

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
####9) VCFLIB

``vcflib`` can be installed and used by ``speedseq`` with the following commands: 
~~~~~~~~~~~~~~~~~~
	git clone --recursive  https://github.com/ekg/vcflib
	cd vcflib
	make
	sudo cp bin/* /usr/local/bin/
~~~~~~~~~~~~~~~~~~
-
**For alternative installations and release issues for any of the above tools please consult the website/creator.**

##Usage

`speedseq` is a modular pipeline with four components: [`aln`](#speedseq-aln), [`var`](#speedseq-var), [`somatic`](#speedseq-somatic), and [`lumpy`](#speedseq-lumpy).

-
###speedseq aln

`speedseq aln` takes paired-end fastq sequences as input, and produces a duplicate-marked, sorted, indexed BAM file that can be processed with other `speedseq` modules. Currently, `speedseq aln` does not support single-end reads.

Internally, `speedseq aln` runs the following steps to produce [three output BAM files](#output):

1. Alignment with [BWA-MEM](http://bio-bwa.sourceforge.net/)
2. Duplicate marking with [samblaster](https://github.com/GregoryFaust/samblaster)
3. Discordant-read and split-read extraction with [samblaster](https://github.com/GregoryFaust/samblaster)
4. Position sorting with [sambamba](https://github.com/lomereiter/sambamba)
5. BAM indexing with [sambamba](https://github.com/lomereiter/sambamba)

~~~~~~~~~~~~~~~~~~
usage:   speedseq aln [options] <reference.fa> <in1.fq> [in2.fq]
~~~~~~~~~~~~~~~~~~

#####Positional arguments

~~~~~~~~~~~~~~~~~~
reference.fa	genome reference fasta file (indexed with bwa) (required)
in1.fq          paired-end fastq file. if -p flag is used then expected to be
                  an interleaved paired-end fastq file, and in2.fq may be omitted.
                  (may be gzipped) (required)
in2.fq	        paired-end fastq file. (may be gzipped) (required)
~~~~~~~~~~~~~~~~~~

#####Alignment options

These options determine the behavior of BWA-MEM
~~~~~~~~~~~~~~~~~~
-o STR		output prefix that will be  [default: in1.fq]
-R              read group header line such as "@RG\tID:libraryname\tSM:samplename" (required)
-p       	first fastq file consists of interleaved paired-end sequences
-t INT          number of threads to use [default: 1]
-T DIR          temp directory [default: ./temp]
~~~~~~~~~~~~~~~~~~

#####Samblaster options

These options determine the behavior of `samblaster`
~~~~~~~~~~~~~~~~~~
-i              include duplicates in splitters and discordants
                  (default: exclude duplicates)
-c INT		maximum number of split alignments for a read to be
                  included in splitter file [default: 2]
-m INT		minimum non-overlapping base pairs between two alignments
                  for a read to be included in splitter file [default: 20]
~~~~~~~~~~~~~~~~~~

#####Global options

~~~~~~~~~~~~~~~~~~
-K FILE         path to speedseq.config file (default: same directory as speedseq)
-h              show help message
~~~~~~~~~~~~~~~~~~

####Output

`speedseq aln` produces three sorted, indexed BAM files (plus their corresponding .bai index files):

* `outprefix.bam`
  * The full, duplicate-marked, sorted BAM file for the library. This file may serve as input for `speedseq var`, `speedseq somatic`, and `speedseq lumpy`.
* `outprefix.splitters.bam`
  * This BAM file contains split reads called by the BWA-MEM alignment of the library. This file excludes duplicate reads by default, but they will be included if the `-i` flag is specified as a `speedseq aln` command line parameter.
* `outprefix.discordants.bam`
  * This BAM file contains discordant read-pairs called by the BWA-MEM alignment of the library. They may be discordant by strand orientation, intrachromosomal distance, or interchromosomal mapping. This file excludes duplicate reads by default, but they will be included if the `-i` flag is specified as a `speedseq aln` command line parameter.

-
###speedseq var

`speedseq var` runs [freebayes](https://github.com/ekg/freebayes) one or more BAM files.

~~~~~~~~~~~~~~~~~~
usage:   speedseq var [options] <reference.fa> <input1.bam> [input2.bam [...]]
~~~~~~~~~~~~~~~~~~

#####Options

~~~~~~~~~~~~~~~~~~
-o STR          output prefix [default: input1.bam]
-w FILE         BED file of windowed genomic intervals. For human genomes,
                  we recommend using the annotations/ceph18.b37.include.2014-01-15.bed
                  BED file to parallelize the variant calling. This BED file excludes
                  regions of the genome where the coverage in the CEPH1463 pedigree
                  was greater than twice the mode coverage plus 5 standard deviations.
                  We believe these extremely high depth regions are areas of misassembly
                  in the GRCh37 human reference genome in which variant calling is
                  time-consuming and error-prone.
-t INT          number of threads to use [default: 1]
-T DIR          temp directory [default: ./temp]
-A BOOL         annotate the vcf with snpEff (true or false) [default: true]
-K FILE         path to speedseq.config file [default: same directory as speedseq]
-h              show help message
~~~~~~~~~~~~~~~~~~

####Output

`speedseq var` produces a single bgzipped VCF file that is indexed with `tabix` and optionally annotated with [SnpEff](http://snpeff.sourceforge.net/):

* `outprefix.vcf.gz`

-
###speedseq somatic

`speedseq somatic` runs [freebayes](https://github.com/ekg/freebayes) on a tumor/normal pair of BAM files

~~~~~~~~~~~~~~~~~~
usage:   speedseq somatic [options] <reference.fa> <normal.bam> <tumor.bam>
~~~~~~~~~~~~~~~~~~

#####Options

~~~~~~~~~~~~~~~~~~
-o STR           output prefix [default: tumor.bam]
-w FILE          BED file of windowed genomic intervals. For human genomes,
                   we recommend using the annotations/ceph18.b37.include.2014-01-15.bed
                   BED file to parallelize the variant calling. This BED file excludes
                   regions of the genome where the coverage in the CEPH1463 pedigree
                   was greater than twice the mode coverage plus 5 standard deviations.
                   We believe these extremely high depth regions are areas of misassembly
                   in the GRCh37 human reference genome in which variant calling is
                   time-consuming and error-prone.
-t INT           number of threads to use [default: 1]
-F FLOAT         require at least this fraction of observations supporting
                   an alternate allele within a single individual in order
                   to evaluate the position [0.05]
-C INT           require at least this count of observations supporting
                   an alternate allele within a single individual in order
                   to evaluate the position [2]
-q FLOAT         minimum variant quality to output [default: 1]
-T DIR           temp directory [./temp]
-A BOOL          annotate the vcf with snpEff (true or false) (default: true)
-K FILE          path to speedseq.config file (default: same directory as speedseq)
-h               show help message
~~~~~~~~~~~~~~~~~~

####Output

`speedseq somatic` produces a single bgzipped VCF file that is indexed with `tabix` and optionally annotated with [SnpEff](http://snpeff.sourceforge.net/):

* `outprefix.vcf.gz`

-
###speedseq lumpy

`speedseq lumpy` runs [lumpy-sv](https://github.com/arq5x/lumpy-sv) on one or more BAM files

#####LUMPY options
~~~~~~~~~~~~~~~~~~
-B FILE          full BAM file(s) (comma separated) (required)
                   example: -B in1.bam,in2.bam,in3.bam
-S FILE          split reads BAM file(s) (comma separated, order same as in -B) (required)
                   example: -B in1.splitters.bam,in2.splitters.bam,in3.splitters.bam
-D FILE          discordant reads BAM files(s) (comma separated, order same as in -B) (required)
-o STR           output prefix [fullBam.bam]
-x FILE          BED file to exclude
-m INT           minimum weight for a call [default: 4]
-r FLOAT         trim threshold [default: 1e-10]
-L INT           read length [default: auto-calculate]
-T DIR           temp directory [default: ./temp]
~~~~~~~~~~~~~~~~~~

The flags `-s` and `-p` are automatically generated using the defaults below, but may be overridden by the user by explicitly defining them using the following format.

~~~~~~~~~~~~~~~~~~
-s STR           lumpy split read parameters [auto]
                    bam_file:<splitreads.bam>,
                    back_distance:<20>,
                    min_mapping_threshold:<20>,
                    weight:<1>,
                    id:<11>,
                    min_clip:<20>

-p STR           lumpy discordant read parameters [auto]
                    bam_file:<discreads.bam>,
                    histo_file:<auto>,
                    mean:<auto>,
                    stdev:<auto>,
                    read_length:<auto>,
                    min_non_overlap:<read_length>,
                    discordant_z:<5>,
                    back_distance:<20>,
                    min_mapping_threshold:<20>,
                    weight:<1>,
                    id:<10>
~~~~~~~~~~~~~~~~~~

#####Global options

~~~~~~~~~~~~~~~~~~
-K FILE          path to speedseq.config file (default: same directory as speedseq)
-h               show help message
~~~~~~~~~~~~~~~~~~

####Output

`speedseq lumpy` produces a BEDPE file.

* `outprefix.bedpe`

The tab-delimited BEDPE file has the following structure:
~~~~~~~~~~~~~~~~~~
1. chromosome 1
2. interval 1 start
3. interval 1 end
4. chromosome 2
5. interval 2 start
6. interval 2 end
7. id
8. evidence set score
9. strand 1
10. strand 2
11. type
12. id of samples containing evidence for this breakpoint
13. strand configurations observed in the evidence set
~~~~~~~~~~~~~~~~~~

Example:
~~~~~~~~~~~~~~~~~~
chr1    34971904    34971945    chr1    34976002    34976043    0x7f9eb0917210  0.0110386   +   -   TYPE:DELETION   IDS:11,1    STRANDS:+-,1
~~~~~~~~~~~~~~~~~~

##Annotations

For human genome alignment using the GRCh37 build, we recommend using the `annotations/ceph18.b37.include.2014-01-15.bed` BED file to parallelize the variant calling (`speedseq var` and `speedseq somatic`). This BED file excludes regions of the genome where the coverage in the CEPH1463 pedigree was greater than twice the mode coverage plus 5 standard deviations. We believe these extremely high depth regions are areas of misassembly in the GRCh37 human reference genome in which variant calling is time-consuming and error-prone.

Additionally, the regions in `annotations/ceph18.b37.include.2014-01-15.bed` are variable-width windows which each contain approximately the same coverage depth in the CEPH1463 pedigree. This ensures that the parallelization of Freebayes uses approximately the same amount of time per region.

In the `speedseq lumpy` module, we recommend excluding the genomic regions in the `annotations/ceph18.b37.exclude.2014-01-15.bed` BED file. These regions represent the complement of those in `annotations/ceph18.b37.include.2014-01-15.bed`.

##Example Workflows

###Call variants on a single sample

1. Use `speedseq aln` to produce a sorted, duplicate-marked, BAM alignment from paired-end fastq data.

  ~~~~~~~~~~~~~~~~~~
  speedseq aln -o NA12878 -R "@RG\tID:NA12878.S1\tSM:NA12878" \
      human_g1k_v37.fasta NA12878.1.fq.gz NA12878.2.fq.gz
  ~~~~~~~~~~~~~~~~~~

  Note: if using an interleaved paired-end fastq file, use the `-p` flag

  ~~~~~~~~~~~~~~~~~~
  speedseq aln -p -o NA12878 -R "@RG\tID:NA12878.S1\tSM:NA12878" \
      human_g1k_v37.fasta NA12878.interleaved.fq.gz
  ~~~~~~~~~~~~~~~~~~

2. Use `speedseq var` to call SNVs and indels on a single sample.

  ~~~~~~~~~~~~~~~~~~
  speedseq var -o NA12878 \
      -w annotations/ceph18.b37.include.2014-01-15.bed \
      human_g1k_v37.fasta NA12878.bam
  ~~~~~~~~~~~~~~~~~~

3. Use `speedseq lumpy` to call structural variants.

  ~~~~~~~~~~~~~~~~~~
  speedseq lumpy -o NA12878 \
      -x annotations/ceph18.b37.exclude.2014-01-15.bed \
      -B NA12878.bam \
      -D NA12878.discordants.bam \
      -S NA12878.splitters.bam
  ~~~~~~~~~~~~~~~~~~

###Call variants on a single sample sequenced with multiple libraries

1. Use `speedseq aln` to produce a sorted, duplicate-marked, BAM alignment of each library.

  ~~~~~~~~~~~~~~~~~~
  speedseq aln -o NA12878_S1 -R "@RG\tID:NA12878.S1\tSM:NA12878" \
      human_g1k_v37.fasta NA12878.S1.1.fq.gz NA12878.S1.2.fq.gz

  speedseq aln -o NA12878_S2 -R "@RG\tID:NA12878.S2\tSM:NA12878" \
      human_g1k_v37.fasta NA12878.S2.1.fq.gz NA12878.S2.2.fq.gz

  speedseq aln -o NA12878_S3 -R "@RG\tID:NA12878.S3\tSM:NA12878" \
      human_g1k_v37.fasta NA12878.S3.1.fq.gz NA12878.S3.2.fq.gz
  ~~~~~~~~~~~~~~~~~~

2. Use `speedseq var` to call SNVs and indels. `speedseq var` will automatically recognize libraries with the same SM readgroup tag to be the same sample.

  ~~~~~~~~~~~~~~~~~~
  speedseq var -o NA12878 \
      -w annotations/ceph18.b37.include.2014-01-15.bed \
      human_g1k_v37.fasta NA12878_S1.bam NA12878_S2.bam NA12878_S3.bam
  ~~~~~~~~~~~~~~~~~~

3. `speedseq lumpy` does not currently handle BAM files made from multiple libraries but we plan to add this functionality in the future.


###Call variants on multiple samples

1. Use `speedseq aln` to produce sorted, duplicate-marked, BAM alignments for each sample.

  ~~~~~~~~~~~~~~~~~~
  speedseq aln -o NA12877 -R "@RG\tID:NA12877.S1\tSM:NA12877" \
      human_g1k_v37.fasta NA12877.1.fq.gz NA12877.2.fq.gz

  speedseq aln -o NA12878 -R "@RG\tID:NA12878.S1\tSM:NA12878" \
      human_g1k_v37.fasta NA12878.1.fq.gz NA12878.2.fq.gz

  speedseq aln -o NA12879 -R "@RG\tID:NA12879.S1\tSM:NA12879" \
      human_g1k_v37.fasta NA12879.1.fq.gz NA12879.2.fq.gz
  ~~~~~~~~~~~~~~~~~~

2. Use `speedseq var` to call SNVs and indels on multiple samples.

  ~~~~~~~~~~~~~~~~~~
  speedseq var -o cephtrio \
      -w annotations/ceph18.b37.include.2014-01-15.bed \
      human_g1k_v37.fasta NA12877.bam NA12878.bam NA12879.bam
  ~~~~~~~~~~~~~~~~~~

3. Use `speedseq lumpy` to call structural variants on multiple samples.

  ~~~~~~~~~~~~~~~~~~
  speedseq lumpy -o cephtrio \
      -x annotations/ceph18.b37.exclude.2014-01-15.bed \
      -B NA12877.bam,NA12878.bam,NA12879.bam \
      -D NA12877.discordants.bam,NA12878.discordants.bam,NA12879.discordants.bam \
      -S NA12877.splitters.bam,NA12878.splitters.bam,NA12879.splitters.bam
  ~~~~~~~~~~~~~~~~~~

###Call variants on a tumor/normal pair

1. Use `speedseq aln` to produce sorted, duplicate-marked, BAM alignments for the tumor/normal pair

  ~~~~~~~~~~~~~~~~~~
  speedseq aln -p -o TCGA-B6-A0I6.normal \
      -R "@RG\tID:TCGA-B6-A0I6-10A-01D-A128-09\tSM:TCGA-B6-A0I6-10A-01D-A128-09" \
      human_g1k_v37.fasta TCGA-B6-A0I6-10A-01D-A128-09.interleaved.fq.gz

  speedseq aln -p -o TCGA-B6-A0I6.tumor \
      -R "@RG\tID:TCGA-B6-A0I6-01A-11D-A128-09\tSM:TCGA-B6-A0I6-01A-11D-A128-09" \
      human_g1k_v37.fasta TCGA-B6-A0I6-01A-11D-A128-09.interleaved.fq.gz
  ~~~~~~~~~~~~~~~~~~

2. Use `speedseq somatic` to call SNVs and indels on the tumor/normal pair.
  ~~~~~~~~~~~~~~~~~~
  speedseq somatic -o TCGA-B6-A0I6 \
      -w annotations/ceph18.b37.include.2014-01-15.bed \
      -F 0.05 \
      -q 1 \
      human_g1k_v37.fasta TCGA-B6-A0I6.normal.bam TCGA-B6-A0I6.tumor.bam
  ~~~~~~~~~~~~~~~~~~

3. Use `speedseq lumpy` to call structural variants on the tumor/normal pair.
  ~~~~~~~~~~~~~~~~~~
  speedseq lumpy -o TCGA-B6-A0I6 \
      -x annotations/ceph18.b37.exclude.2014-01-15.bed \
      -B TCGA-B6-A0I6.normal.bam,TCGA-B6-A0I6.tumor.bam \
      -D TCGA-B6-A0I6.normal.discordants.bam,TCGA-B6-A0I6.tumor.discordants.bam \
      -S TCGA-B6-A0I6.normal.splitters.bam,TCGA-B6-A0I6.tumor.splitters.bam
  ~~~~~~~~~~~~~~~~~~


