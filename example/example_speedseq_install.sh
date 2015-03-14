# Install SpeedSeq on blank Amazon Linux machine
# running on Amazon Linux AMI 2014.09.2 (HVM) - ami-146e2a7c
# m3.large, 2 vCPUs, 7.5 GB memory

# ---------------------------------
# 1. Install prerequisites
# ---------------------------------

# SpeedSeq prerequisites
sudo yum -y update
sudo yum -y install make automake cmake gcc gcc-c++ git ncurses-devel zlib-devel

# Get Python 2.7 and modules
sudo yum install -y python27 python27-devel python27-pip lapack lapack-devel blas blas-devel
sudo pip-2.7 install pysam numpy scipy

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
echo "source /usr/local/root/bin/thisroot.sh" >> ~/.bashrc
source ~/.bashrc

# clean up
rm root_v5.34.20.source.tar.gz

# ---------------------------------
# 2. Install SpeedSeq core components
# ---------------------------------

# Get SpeedSeq (tested on commit cb097f33f37710c6be1ffe9256ee39983a776769)
mkdir code
cd code
git clone --recursive https://github.com/cc2qe/speedseq.git
cd speedseq
make

# Add speedseq to $PATH
echo -e "export PATH=\$HOME/code/speedseq/bin:\$PATH" >> ~/.bashrc
source ~/.bashrc

# ---------------------------------
# 3. Install VEP
# ---------------------------------

# Get required perl modules
sudo yum -y install "perl(Archive::Extract)" "perl(CGI)" "perl(DBI)" "perl(Time::HiRes)" "perl(Archive::Tar)" "perl(Archive::Zip)"

# Download and install VEP
cd ~
curl -OL https://github.com/Ensembl/ensembl-tools/archive/release/76.zip
unzip 76.zip
perl ensembl-tools-release-76/scripts/variant_effect_predictor/INSTALL.pl \
    -c ~/code/speedseq/annotations/vep_cache \
    -a ac -s homo_sapiens -y GRCh37

# Copy executables to SpeedSeq directory
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

# ---------------------------------
# 5. Install GEMINI
# ---------------------------------
wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py
python gemini_install.py /usr/local /usr/local/share/gemini
echo -e "export PATH=\$PATH:/usr/local/gemini/bin" >> ~/.bashrc
source ~/.bashrc






