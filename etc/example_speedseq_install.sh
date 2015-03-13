# Install SpeedSeq on blank Amazon Linux machine

# ---------------------------------
# 1. Install prerequisites
# ---------------------------------

# SpeedSeq prerequisites
sudo yum -y update
sudo yum -y install cmake gcc-c++ gcc git make python27 python-devel python-yaml ncurses-devel zlib-devel numpy python-pip
sudo yum -y install "perl(Archive::Extract)" "perl(CGI)" "perl(DBI)" "perl(Time::HiRes)" "perl(Archive::Tar)" "perl(Archive::Zip)"

# CNVnator prerequisites
sudo yum -y install libX11-devel libXpm-devel libXft-devel libXext-devel
curl -OL ftp://root.cern.ch/root/root_v5.34.20.source.tar.gz
tar -xvf root_v5.34.20.source.tar.gz
cd root
./configure
make
cd ..
sudo mv root /usr/local

# Source the file and add ROOT source to ~/.bashrc
source /usr/local/root/bin/thisroot.sh
echo "source /usr/local/root/bin/thisroot.sh" >> ~/.bashrc

# Remove the tarball
rm root_v5.34.20.source.tar.gz

# ---------------------------------
# 2. Install SpeedSeq core components
# ---------------------------------

# Get SpeedSeq (commit 5950b8bb2b5170aaa65837517f6cad20ee574720)
mkdir code
cd code
git clone --recursive https://github.com/cc2qe/speedseq.git
cd speedseq
make

# Add speedseq to $PATH
echo -e "export PATH=\$HOME/code/speedseq/bin:\$PATH" >> ~/.bashrc
. ~/.bashrc

# ---------------------------------
# 3. Install VEP
# ---------------------------------
cd ~
curl -OL https://github.com/Ensembl/ensembl-tools/archive/release/76.zip
unzip 76.zip
perl ensembl-tools-release-76/scripts/variant_effect_predictor/INSTALL.pl \
    -c ~/code/speedseq/annotations/vep_cache \
    -a ac -s homo_sapiens -y GRCh37

cp ensembl-tools-release-76/scripts/variant_effect_predictor/variant_effect_predictor.pl ~/code/speedseq/bin/.
mv Bio ~/code/speedseq/bin/.

# clean up
rm -r 76.zip ensembl-tools-release-76

# ---------------------------------
# 4. Prepare reference genome files
# ---------------------------------

# Get human reference GRCh37
mkdir -p ~/genomes
cd ~/genomes
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai
gunzip human_g1k_v37.fasta.gz
bwa index ~/genomes/human_g1k_v37.fasta

# make the CNVnator chroms
mkdir -p ~/code/speedseq/annotations/cnvnator_chroms
cd ~/code/speedseq/annotations/cnvnator_chroms
cat ~/genomes/human_g1k_v37.fasta | awk 'BEGIN { CHROM="" } { if ($1~"^>") CHROM=substr($1,2); print $0 > CHROM".fa" }'



