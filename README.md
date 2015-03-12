# SpeedSeq         

SpeedSeq is a flexible and open-source framework to rapidly identify genomic variation.

Chiang, C, R M Layer, G G Faust, M R Lindberg, D B Rose, E P Garrison, G T Marth, A R Quinlan, and I M Hall. 2014. SpeedSeq: Ultra-Fast Personal Genome Analysis and Interpretation. bioRxiv. [doi:10.1101/012179](http://dx.doi.org/10.1101/012179).

##Table of Contents

1. [Quick start](#quick-start)
2. [Installation](#installation)

4. [Usage](#usage)
5. [Annotations](#annotations)
6. [Example Workflows](#example-workflows)

## Quick start

1. Install
	```
	git clone --recursive https://github.com/cc2qe/speedseq
	cd speedseq
	make
	```

2. Align from FASTQ
	```
	bin/speedseq align \
		-o sample \
		-R "@RG\tID:id\tSM:samplename\tLB:lib" \
		human_g1k_v37.fasta \
		in1.fq.gz \
		in2.fq.gz
	```

3. Detect SNVs and indels
	```
	bin/speedseq var \
		-o sample \
		human_g1k_v37.fasta \
		sample.bam
	```

4. Detect structural variants
	```
	bin/speedseq sv \
		-o sample \
		-B sample.bam \
		-S sample.splitters.bam \
		-D sample.discordants.bam \
		-R human_g1k_v37.fasta \
	```

## Installation

#### Prerequisites
* Python 2.7
	* numpy
	* pysam
	* scipy
* ROOT (required if running CNVnator)
* Variant Effect Predictor (required if annotating VCF files)

```
git clone --recursive https://github.com/cc2qe/speedseq
cd speedseq
make
```


## Components

* [BWA](http://bio-bwa.sourceforge.net/)
* [SAMBLASTER](https://github.com/GregoryFaust/samblaster)
* [Sambamba](https://github.com/lomereiter/sambamba)
* [FreeBayes](https://github.com/ekg/freebayes)
* [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html)
* [LUMPY](https://github.com/arq5x/lumpy-sv)
* [SVtyper](https://github.com/cc2qe/svtyper)
* [CNVnator](http://sv.gersteinlab.org/cnvnator/)
* [GEMINI](https://github.com/arq5x/gemini)
* [GNU Parallel](http://www.gnu.org/software/parallel/)

##Installation

Preparing genome files

bwa mem index

Current support for Linux only

There is an automatic and manual installation process for SpeedSeq.

The following are required for both installations:

* cmake
* g++ 4.3 or later
* gcc
* git
* make
* python2.7
* python-devel
* python-yaml
* ncurses-devel
* zlib-devel
* numpy
* pip

A Linux package manager can be used to obtain these by with the command:

```
# general purpose tools
sudo yum update
sudo yum -y install cmake gcc-c++ gcc git make python27 python-devel python-yaml ncurses-devel zlib-devel numpy python-pip

# PERL modules for VEP
sudo yum -y install "perl(Archive::Extract)" "perl(CGI)" "perl(DBI)" "perl(Time::HiRes)" "perl(Archive::Tar)" "perl(Archive::Zip)"
```
 
or 

```
sudo apt-get update
sudo apt-get -y install build-essential cmake gpp gcc git make python2.7 python-dev python-yaml ncurses-dev zlib1g-dev python-numpy python-pip
```

###Automatic installation

SpeedSeq can be installed with the following commands:
```
git clone --recursive https://github.com/cc2qe/speedseq
cd speedseq
python speedseq_setup.py
sudo cp -r bin/* /usr/local/bin/
```

Configure the paths to the SpeedSeq dependencies by modifying the [speedseq.config](bin/speedseq.config) file. The [speedseq.config](bin/speedseq.config) file should reside in the same directory as the SpeedSeq executable. By default, SpeedSeq attempts to source the dependencies from the $PATH.

Note that the optional component CNVnator cannot be automatically installed (see the [CNVnator installation instructions](#cnvnator)). SpeedSeq uses a multi-threaded implementation of CNVnator v0.3, which can be for in [src/cnvnator](src/cnvnator)

####Genome
We recommend using the GRCh37 human genome for SpeedSeq, available here:

* ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
* ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai

###Manual installation

The following instructions for installation assumes that the required tools are not installed. If you already have the SpeedSeq components on your system, then you can simply link them to the proper paths in the [speedseq.config](bin/speedseq.config] file.
It is recommended that the specified versions of each tool are used for this release of SpeedSeq.  
The use of unspecified versions of any pipeline component is not guaranteed to work. 

SpeedSeq can be installed with the following commands: 
```
git clone --recursive https://github.com/cc2qe/speedseq
cd speedseq
make
sudo cp -r bin/* /usr/local/bin/
```

#####Configuration

Configure the paths to the SpeedSeq dependencies by modifying the [speedseq.config](bin/speedseq.config) file. The [speedseq.config](bin/speedseq.config) file should reside in the same directory as the SpeedSeq executable. By default, SpeedSeq attempts to source the dependencies from the $PATH.

#####Genome

We recommend using the GRCh37 human genome for SpeedSeq, available here:

* ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
* ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai

####Obtain each of the external pipeline tools and install:
	
#### Sambamba
https://github.com/lomereiter/sambamba

Sambamba can be installed and used by SpeedSeq by: 
```
curl -OL https://github.com/lomereiter/sambamba/releases/download/v0.4.6-beta/sambamba_v0.4.6-beta_centos5-x86_64.tar.bz2
tar -xvf sambamba_v0.4.6-beta_centos5-x86_64.tar.bz2
sudo cp sambamba_v0.4.6 /usr/local/bin/
```

#### CNVnator
http://sv.gersteinlab.org/cnvnator/

CNVnator can be installed and used by SpeedSeq with the following commands:

1. Install the ROOT package (http://root.cern.ch/drupal/)
   ```
   # install dependencies
   sudo yum -y install libX11-devel libXpm-devel libXft-devel libXext-devel

   # get ROOT
   curl -OL ftp://root.cern.ch/root/root_v5.34.20.source.tar.gz
   tar -xvf root_v5.34.20.source.tar.gz
   cd root
   ./configure
   make
   cd ..
   sudo mv root /usr/local
   ```

2. Add the following line to `~/.bashrc`
   ```
   source /usr/local/root/bin/thisroot.sh
   ```

   Then run the following for reload `.bashrc`
   ```
   . ~/.bashrc
   ```

3. Navigate to SpeedSeq git directory and compile the multi-threaded implementation of CNVnator v0.3
   ```
   cd speedseq
   make cnvnator-multi
   sudo cp bin/* /usr/local/bin
   ```

4. Create chromosomes directory from genome fasta file

   CNVnator requires a directory containing each chromosome as a separate fasta file. Run the following commands to generate this directory.
   ```
   mkdir chroms
   cd chroms
   cat human_g1k_v37.fasta | awk 'BEGIN { CHROM="" } { if ($1~"^>") CHROM=substr($1,2); print $0 > CHROM".fa" }'
   ```

   Then set the path to the `chroms` directory as `CNVNATOR_CHROMS_DIR` in [speedseq.config](bin/speedseq.config)

#### VEP
http://www.ensembl.org/info/docs/tools/vep/index.html

VEP can be installed and used by SpeedSeq with the following commands:

1. Install required PERL modules
   ```
   sudo yum -y install "perl(Archive::Extract)" "perl(CGI)" "perl(DBI)" "perl(Time::HiRes)" "perl(Archive::Tar)" "perl(Archive::Zip)"
   ```

2. Download the software and install
   ```
   curl -OL https://github.com/Ensembl/ensembl-tools/archive/release/76.zip
   unzip 76.zip
   perl ensembl-tools-release-76/scripts/variant_effect_predictor/INSTALL.pl -a ac -s homo_sapiens -y GRCh37 --CACHEDIR ~/.vep
   ```

3. Copy files to a directory in $PATH
   ```
   sudo cp ensembl-tools-release-76/scripts/variant_effect_predictor/variant_effect_predictor.pl /usr/local/bin
   sudo cp -r ensembl-tools-release-76/scripts/variant_effect_predictor/Bio /usr/local/bin
   ```

4. Edit [speedseq.config](bin/speedseq.config) to set correct paths for the `VEP` and `VEP_CACHE_DIR` variables.

#### GEMINI
https://github.com/arq5x/gemini

GEMINI can be installed and used by SpeedSeq with the following commands: 
```
wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py	
# or curl -OL https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py > gemini_install.py
sudo python2.7 gemini_install.py /usr/local /usr/local/share/gemini
export PATH=$PATH:/usr/local/gemini/bin
# it would be wise to add the above line to your ``.bashrc`` or ``.bash_profile``
gemini update
```

#### BEDTools
https://github.com/arq5x/bedtools2

BEDTools can be installed and used by SpeedSeq with the following commands:
```
curl -OL https://github.com/arq5x/bedtools2/releases/download/v2.20.1/bedtools-2.20.1.tar.gz
tar -xvf bedtools-2.20.1.tar.gz
make -C bedtools2-2.20.1
sudo cp bedtools2-2.20.1/bin/* /usr/local/bin
```

#### GNU Parallel
http://www.gnu.org/software/parallel/

Parallel can be installed and used by SpeedSeq with the following commands:
```
curl -OL http://ftp.gnu.org/gnu/parallel/parallel-20100424.tar.bz2
tar -xvf parallel-20100424.tar.bz2
cd parallel-20100424
./configure && sudo make && sudo make install
sudo cp src/parallel /usr/local/bin/
```

**For alternative installations and release issues for any of the above tools please consult the website/creator.**

##Usage

SpeedSeq is a modular pipeline with four components: [`align`](#speedseq-align), [`var`](#speedseq-var), [`somatic`](#speedseq-somatic), and [`sv`](#speedseq-sv).

* [`speedseq align`](#speedseq-align)
  * Take paired-end fastq sequences as input, and produce a duplicate-marked, sorted, indexed BAM file that can be processed with other speedseq modules. Currently, speedseq align does not support single-end reads.
* [`speedseq var`](#speedseq-var)
  * Run FreeBayes one or more BAM files
* [`speedseq somatic`](#speedseq-somatic)
  * Run FreeBayes on a tumor/normal pair of BAM files
* [`speedseq sv`](#speedseq-sv)
  * Run LUMPY on one or more BAM files, with breakend genotyping and read-depth calculation.

###speedseq align

`speedseq align` takes paired-end fastq sequences as input, and produces a duplicate-marked, sorted, indexed BAM file that can be processed with other SpeedSeq modules. Currently, `speedseq align` does not support single-end reads.

Internally, `speedseq align` runs the following steps to produce [three output BAM files](#output):

1. Alignment with BWA-MEM
   * Prior to alignment, run `bwa index genome.fasta` on the reference genome
2. Duplicate marking with SAMBLASTER
3. Discordant-read and split-read extraction with SAMBLASTER
4. Position sorting with Sambamba
5. BAM indexing with Sambamba

```
usage:   speedseq align [options] <reference.fa> <in1.fq> [in2.fq]
```

#####Positional arguments

```
reference.fa	genome reference fasta file (indexed with bwa) (required)
in1.fq          paired-end fastq file. if -p flag is used then expected to be
                  an interleaved paired-end fastq file, and in2.fq may be omitted.
                  (may be gzipped) (required)
in2.fq	        paired-end fastq file. (may be gzipped) (required)
```

#####Alignment options

These options determine the behavior of BWA-MEM
```
-o STR          output prefix [default: in1.fq]
-R              read group header line such as "@RG\tID:libraryname\tSM:samplename" (required)
-p              first fastq file consists of interleaved paired-end sequences
-t INT          number of threads to use [default: 1]
-T DIR          temp directory [default: ./temp]
```

#####Samblaster options

These options determine the behavior of SAMBLASTER
```
-i              include duplicates in splitters and discordants
                  (default: exclude duplicates)
-c INT          maximum number of split alignments for a read to be
                  included in splitter file [default: 2]
-m INT          minimum non-overlapping base pairs between two alignments
                for a read to be included in splitter file [default: 20]
```

#####Sambamba options
```
-M              amount of memory in GB to be used for sorting [default: 20]
```

#####Global options

```
-K FILE         path to speedseq.config file (default: same directory as speedseq)
-v              verbose
-h              show help message
```

####Output

`speedseq align` produces three sorted, indexed BAM files (plus their corresponding .bai index files):

* `outprefix.bam`
  * The full, duplicate-marked, sorted BAM file for the library. This file may serve as input for [`speedseq var`](#speedseq-var), [`speedseq somatic`](#speedseq-somatic), and [`speedseq sv`](#speedseq-sv).
* `outprefix.splitters.bam`
  * This BAM file contains split reads called by the BWA-MEM alignment of the library. It may be used as the `-S` flag input to [`speedseq sv`](#speedseq-sv). This file excludes duplicate reads by default, but they will be included if the `-i` flag is specified as a [`speedseq align`](#speedseq-align) command line parameter.
* `outprefix.discordants.bam`
  * This BAM file contains discordant read-pairs called by the BWA-MEM alignment of the library. These reads may be discordant by strand orientation, intrachromosomal distance, or interchromosomal mapping. This BAM file may be used as the `-D` flag input to [`speedseq sv`](#speedseq-sv). This file excludes duplicate reads by default, but they will be included if the `-i` flag is specified as a [`speedseq align`](#speedseq-align) command line parameter.

###speedseq var

`speedseq var` runs FreeBayes one or more BAM files.

```
usage:   speedseq var [options] <reference.fa> <input1.bam> [input2.bam [...]]
```

#####Positional arguments

```
reference.fa    genome reference fasta file
input.bam       BAM file(s) to call variants on. Must have readgroup information,
                  and the SM readgroup tags will be the VCF column header
```

#####Options

```
-o STR          output prefix [default: input1.bam]
-w FILE         BED file of windowed genomic intervals. For human genomes,
                  we recommend using the annotations/ceph18.b37.include.2014-01-15.bed
                  (see Annotations)
-q FLOAT        minimum variant QUAL score to output [1]
-t INT          number of threads to use [default: 1]
-T DIR          temp directory [default: ./temp]
-A              annotate the vcf with VEP
-K FILE         path to speedseq.config file [default: same directory as speedseq]
-v              verbose
-h              show help message
```

####Output

`speedseq var` produces a single indexed VCF file that is optionally annotated with VEP.

* `outprefix.vcf.gz`

###speedseq somatic

`speedseq somatic` runs FreeBayes on a tumor/normal pair of BAM files

```
usage:   speedseq somatic [options] <reference.fa> <normal.bam> <tumor.bam>
```

#####Positional arguments

```
reference.fa      genome reference fasta file
normal.bam        germline BAM file(s) (comma separated BAMs from multiple libraries).
                    Must have readgroup information, and the SM readgroup tag will
                    be the VCF column header
tumor.bam         tumor BAM file(s) (comma separated BAMs for multiple libraries).
                    Must have readgroup information, and the SM readgroup tag will
                    be the VCF column header
```

#####Options

```
-o STR           output prefix [default: tumor.bam]
-w FILE          BED file of windowed genomic intervals. For human genomes,
                   we recommend using the annotations/ceph18.b37.include.2014-01-15.bed
                   (see Annotations)
-t INT           number of threads to use [default: 1]
-F FLOAT         require at least this fraction of observations supporting
                   an alternate allele within a single individual in order
                   to evaluate the position [0.05]
-C INT           require at least this count of observations supporting
                   an alternate allele within a single individual in order
                   to evaluate the position [2]
-n FLOAT         minimum normal log odds ratio for PASS [7]
-u FLOAT         minimum tumor log odds ratio for PASS [7]
-q FLOAT         minimum QUAL score to output non-passing somatic variants [1e-5]
-T DIR           temp directory [./temp]
-A               annotate the vcf with VEP
-K FILE          path to speedseq.config file (default: same directory as speedseq)
-v               verbose
-h               show help message
```

####Output

`speedseq somatic` produces a single indexed VCF file that is optionally annotated with VEP.

* `outprefix.vcf.gz`

###speedseq sv

`speedseq sv` runs LUMPY on one or more BAM files, with optional breakend genotyping by SVtyper, and optional read-depth analysis by CNVnator.

#####Options
```
-B FILE          full BAM file(s) (comma separated) (required)
                   example: -B in1.bam,in2.bam,in3.bam
-S FILE          split reads BAM file(s) (comma separated, order same as in -B) (required)
                   example: -S in1.splitters.bam,in2.splitters.bam,in3.splitters.bam
-D FILE          discordant reads BAM file(s) (comma separated, order same as in -B) (required)
                   example: -D in1.discordants.bam,in2.discordants.bam,in3.discordants.bam
-R FILE          indexed reference genome fasta file (required)
-o STR           output prefix [in1.bam]
-t INT           threads [1]
-x FILE          BED file to exclude
-g               genotype SV breakends with svtyper
-d               calculate read-depth with CNVnator
-A               annotate the vcf with VEP
-m INT           minimum weight for a call [default: 4]
-r FLOAT         trim threshold [0]
-T DIR           temp directory [default: ./temp]
-k               keep temporary files
```

The flags `-s` and `-p` are automatically generated using the defaults below, but may be overridden by the user by explicitly defining them using the following format.

```
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
```

#####Global options

```
-K FILE          path to speedseq.config file (default: same directory as speedseq)
-v               verbose
-h               show help message
```

####Output

`speedseq sv` produces a VCF file and a BEDPE file containing the same information.

* `outprefix.sv.vcf.gz`
* `outprefix.sv.bedpe`

The tab-delimited BEDPE file has the following structure:
```
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
```

Example:
```
chr1    34971904    34971945    chr1    34976002    34976043    0x7f9eb0917210  0.0110386   +   -   TYPE:DELETION   IDS:11,1    STRANDS:+-,1
```

##Annotations

For human genome alignment using the GRCh37 build, we recommend using the [annotations/ceph18.b37.include.2014-01-15.bed](annotations/ceph18.b37.include.2014-01-15.bed) BED file to parallelize the variant calling ([`speedseq var`](#speedseq-var) and [`speedseq somatic`](#speedseq-somatic)). This BED file excludes regions of the genome where the coverage in the CEPH1463 pedigree was greater than twice the mode coverage plus 3 standard deviations. We believe these extremely high depth regions that we excluded are areas of misassembly in the GRCh37 human reference genome in which variant calling is time-consuming and error-prone.

Additionally, the regions in [annotations/ceph18.b37.include.2014-01-15.bed](annotations/ceph18.b37.include.2014-01-15.bed) are variable-width windows which each contain approximately the same coverage depth in the CEPH1463 pedigree, and sorted from highest to lowest depth. This ensures that the parallelization of Freebayes uses approximately the same amount of time per region.

The regions in [annotations/ceph18.b37.exclude.2014-01-15.bed](annotations/ceph18.b37.exclude.2014-01-15.bed) represent the complement of the regions in [annotations/ceph18.b37.include.2014-01-15.bed](annotations/ceph18.b37.include.2014-01-15.bed).

In the [`speedseq sv`](#speedseq-sv) module, we recommend excluding the genomic regions in the [annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed](annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed) BED file. These regions represent the complement of those in [annotations/ceph18.b37.include.2014-01-15.bed](annotations/ceph18.b37.include.2014-01-15.bed) as well as the mitochondrial chromosome.

##Example Workflows

###Call variants on a single sample

1. Use `speedseq align` to produce a sorted, duplicate-marked, BAM alignment from paired-end fastq data.

  ```
  speedseq align -o NA12878 -R "@RG\tID:NA12878.S1\tSM:NA12878" \
      human_g1k_v37.fasta NA12878.1.fq.gz NA12878.2.fq.gz
  ```

  Note: if using an interleaved paired-end fastq file, use the `-p` flag

  ```
  speedseq align -p -o NA12878 -R "@RG\tID:NA12878.S1\tSM:NA12878" \
      human_g1k_v37.fasta NA12878.interleaved.fq.gz
  ```

2. Use `speedseq var` to call SNVs and indels on a single sample.

  ```
  speedseq var -o NA12878 \
      -w annotations/ceph18.b37.include.2014-01-15.bed \
      human_g1k_v37.fasta NA12878.bam
  ```

3. Use `speedseq sv` to call structural variants. The optional `-g` and `-d` flags perform breakend genotyping and read-depth calculation respectively

  ```
  speedseq sv -o NA12878 \
      -x annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed \
      -g \
      -d \
      -B NA12878.bam \
      -D NA12878.discordants.bam \
      -S NA12878.splitters.bam
  ```

###Call variants on a single sample sequenced with multiple libraries

1. Use `speedseq align` to produce a sorted, duplicate-marked, BAM alignment of each library.

  ```
  speedseq align -o NA12878_S1 -R "@RG\tID:NA12878.S1\tSM:NA12878" \
      human_g1k_v37.fasta NA12878.S1.1.fq.gz NA12878.S1.2.fq.gz

  speedseq align -o NA12878_S2 -R "@RG\tID:NA12878.S2\tSM:NA12878" \
      human_g1k_v37.fasta NA12878.S2.1.fq.gz NA12878.S2.2.fq.gz

  speedseq align -o NA12878_S3 -R "@RG\tID:NA12878.S3\tSM:NA12878" \
      human_g1k_v37.fasta NA12878.S3.1.fq.gz NA12878.S3.2.fq.gz
  ```

2. Use `speedseq var` to call SNVs and indels. `speedseq var` will automatically recognize libraries with the same SM readgroup tag to be the same sample.

  ```
  speedseq var -o NA12878 \
      -w annotations/ceph18.b37.include.2014-01-15.bed \
      human_g1k_v37.fasta NA12878_S1.bam NA12878_S2.bam NA12878_S3.bam
  ```

3. Use `speedseq sv` to call structural variants. (For proper library merging, each library from the same sample must have the same SM readgroup tag.)
  ```
  speedseq -sv -o NA12878 \
      -x annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed \
      -B NA12878_S1.bam,NA12878_S2.bam,NA12878_S3.bam
      -S NA12878_S1.splitters.bam,NA12878_S2.splitters.bam,NA12878_S3.splitters.bam
      -D NA12878_S1.discordants.bam,NA12878_S2.discordants.bam,NA12878_S3.discordants.bam
  ```

###Call variants on multiple samples

1. Use `speedseq align` to produce sorted, duplicate-marked, BAM alignments for each sample.

  ```
  speedseq align -o NA12877 -R "@RG\tID:NA12877.S1\tSM:NA12877" \
      human_g1k_v37.fasta NA12877.1.fq.gz NA12877.2.fq.gz

  speedseq align -o NA12878 -R "@RG\tID:NA12878.S1\tSM:NA12878" \
      human_g1k_v37.fasta NA12878.1.fq.gz NA12878.2.fq.gz

  speedseq align -o NA12879 -R "@RG\tID:NA12879.S1\tSM:NA12879" \
      human_g1k_v37.fasta NA12879.1.fq.gz NA12879.2.fq.gz
  ```

2. Use `speedseq var` to call SNVs and indels on multiple samples.

  ```
  speedseq var -o cephtrio \
      -w annotations/ceph18.b37.include.2014-01-15.bed \
      human_g1k_v37.fasta NA12877.bam NA12878.bam NA12879.bam
  ```

3. Use `speedseq sv` to call structural variants on multiple samples.

  ```
  speedseq sv -o cephtrio \
      -x annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed \
      -B NA12877.bam,NA12878.bam,NA12879.bam \
      -D NA12877.discordants.bam,NA12878.discordants.bam,NA12879.discordants.bam \
      -S NA12877.splitters.bam,NA12878.splitters.bam,NA12879.splitters.bam
  ```

###Call variants on a tumor/normal pair

1. Use `speedseq align` to produce sorted, duplicate-marked, BAM alignments for the tumor/normal pair

  ```
  speedseq align -p -o TCGA-B6-A0I6.normal \
      -R "@RG\tID:TCGA-B6-A0I6-10A-01D-A128-09\tSM:TCGA-B6-A0I6-10A-01D-A128-09" \
      human_g1k_v37.fasta TCGA-B6-A0I6-10A-01D-A128-09.interleaved.fq.gz

  speedseq align -p -o TCGA-B6-A0I6.tumor \
      -R "@RG\tID:TCGA-B6-A0I6-01A-11D-A128-09\tSM:TCGA-B6-A0I6-01A-11D-A128-09" \
      human_g1k_v37.fasta TCGA-B6-A0I6-01A-11D-A128-09.interleaved.fq.gz
  ```

2. Use `speedseq somatic` to call SNVs and indels on the tumor/normal pair.
  ```
  speedseq somatic -o TCGA-B6-A0I6 \
      -w annotations/ceph18.b37.include.2014-01-15.bed \
      -F 0.05 \
      -q 1 \
      human_g1k_v37.fasta TCGA-B6-A0I6.normal.bam TCGA-B6-A0I6.tumor.bam
  ```

3. Use `speedseq sv` to call structural variants on the tumor/normal pair.
  ```
  speedseq sv -o TCGA-B6-A0I6 \
      -x annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed \
      -B TCGA-B6-A0I6.normal.bam,TCGA-B6-A0I6.tumor.bam \
      -D TCGA-B6-A0I6.normal.discordants.bam,TCGA-B6-A0I6.tumor.discordants.bam \
      -S TCGA-B6-A0I6.normal.splitters.bam,TCGA-B6-A0I6.tumor.splitters.bam
  ```


