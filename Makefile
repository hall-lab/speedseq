export MKFILE_DIR = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
TARGET_BIN=bin
ANNOTATIONS_DIR=annotations
SRC=src

SPEEDSEQ_DIR=$(shell pwd)
BWA_DIR=$(SRC)/bwa
SAMBLASTER_DIR=$(SRC)/samblaster
FREEBAYES_DIR=$(SRC)/freebayes
LUMPY_DIR=$(SRC)/lumpy-sv
SVTYPER_DIR=$(SRC)/svtyper
CNVNATOR_DIR=$(SRC)/cnvnator
TABIX_DIR=$(SRC)/tabix
VAWK_DIR=$(SRC)/vawk
MBUFFER_DIR=$(SRC)/mbuffer
PARALLEL_DIR=$(SRC)/parallel
BAMKIT_DIR=$(SRC)/bamkit

# all
all:
	@echo "" > $(MKFILE_DIR)/install.log
	@echo "Installing align module..." >> $(MKFILE_DIR)/install.log
	$(MAKE) align
	@echo "Done." >> $(MKFILE_DIR)/install.log

	@echo "Installing var and somatic modules..." >> $(MKFILE_DIR)/install.log
	$(MAKE) var
	@echo "Done." >> $(MKFILE_DIR)/install.log

	@echo "Installing sv module..." >> $(MKFILE_DIR)/install.log
	$(MAKE) sv
	@echo "Done." >> $(MKFILE_DIR)/install.log

	@echo "Installing realign module..." >> $(MKFILE_DIR)/install.log
	$(MAKE) realign
	@echo "Done." >> $(MKFILE_DIR)/install.log

	@echo "Installation successful" >> $(MKFILE_DIR)/install.log

# modules
align: bwa sambamba samblaster parallel config

var: freebayes tabix vawk parallel config

somatic: var

sv: lumpy sambamba samblaster vawk bamkit tabix svtyper cnvnator-multi config

realign: bwa sambamba samblaster parallel mbuffer bamkit config

# autogenerate speedseq.config
config:
	@echo "" > $(TARGET_BIN)/speedseq.config
	@echo "SPEEDSEQ_HOME=$(MKFILE_DIR)" >> $(TARGET_BIN)/speedseq.config
	@echo "" >> $(TARGET_BIN)/speedseq.config
	@echo "# general" >> $(TARGET_BIN)/speedseq.config
	@echo "SAMBAMBA=$(MKFILE_DIR)/$(TARGET_BIN)/sambamba" >> $(TARGET_BIN)/speedseq.config
	@echo "BGZIP=$(MKFILE_DIR)/$(TARGET_BIN)/bgzip" >> $(TARGET_BIN)/speedseq.config
	@echo "TABIX=$(MKFILE_DIR)/$(TARGET_BIN)/tabix" >> $(TARGET_BIN)/speedseq.config
	@echo "VAWK=$(MKFILE_DIR)/$(TARGET_BIN)/vawk" >> $(TARGET_BIN)/speedseq.config
	@echo "PARALLEL=$(MKFILE_DIR)/$(TARGET_BIN)/parallel" >> $(TARGET_BIN)/speedseq.config
	@echo "PYTHON=`which python2.7`" >> $(TARGET_BIN)/speedseq.config

	@echo "" >> $(TARGET_BIN)/speedseq.config
	@echo "# align" >> $(TARGET_BIN)/speedseq.config
	@echo "BWA=$(MKFILE_DIR)/$(TARGET_BIN)/bwa" >> $(TARGET_BIN)/speedseq.config
	@echo "SAMBLASTER=$(MKFILE_DIR)/$(TARGET_BIN)/samblaster" >> $(TARGET_BIN)/speedseq.config

	@echo "" >> $(TARGET_BIN)/speedseq.config
	@echo "# var/somatic" >> $(TARGET_BIN)/speedseq.config
	@echo "FREEBAYES=$(MKFILE_DIR)/$(TARGET_BIN)/freebayes" >> $(TARGET_BIN)/speedseq.config
	@echo "VEP=$(MKFILE_DIR)/$(TARGET_BIN)/variant_effect_predictor.pl" >> $(TARGET_BIN)/speedseq.config
	@echo "VEP_CACHE_DIR=$(MKFILE_DIR)/$(ANNOTATIONS_DIR)/vep_cache" >> $(TARGET_BIN)/speedseq.config

	@echo "" >> $(TARGET_BIN)/speedseq.config
	@echo "# sv" >> $(TARGET_BIN)/speedseq.config
	@echo "LUMPY=$(MKFILE_DIR)/$(TARGET_BIN)/lumpy" >> $(TARGET_BIN)/speedseq.config
	@echo "LUMPYEXPRESS=$(MKFILE_DIR)/$(TARGET_BIN)/lumpyexpress" >> $(TARGET_BIN)/speedseq.config
	@echo "PAIREND_DISTRO=$(MKFILE_DIR)/$(TARGET_BIN)/pairend_distro.py" >> $(TARGET_BIN)/speedseq.config
	@echo "SVTYPER=$(MKFILE_DIR)/$(TARGET_BIN)/svtyper" >> $(TARGET_BIN)/speedseq.config
	@echo "BAMGROUPREADS=$(MKFILE_DIR)/$(TARGET_BIN)/bamgroupreads.py" >> $(TARGET_BIN)/speedseq.config
	@echo "BAMFILTERRG=$(MKFILE_DIR)/$(TARGET_BIN)/bamfilterrg.py" >> $(TARGET_BIN)/speedseq.config
	@echo "BAMLIBS=$(MKFILE_DIR)/$(TARGET_BIN)/bamlibs.py" >> $(TARGET_BIN)/speedseq.config

	@echo "" >> $(TARGET_BIN)/speedseq.config
	@echo "# CNVnator" >> $(TARGET_BIN)/speedseq.config
	@echo "CNVNATOR_WRAPPER=$(MKFILE_DIR)/$(TARGET_BIN)/cnvnator_wrapper.py" >> $(TARGET_BIN)/speedseq.config
	@echo "CNVNATOR_MULTI=$(MKFILE_DIR)/$(TARGET_BIN)/cnvnator-multi" >> $(TARGET_BIN)/speedseq.config
	@echo "ANNOTATE_RD=$(MKFILE_DIR)/$(TARGET_BIN)/annotate_rd.py" >> $(TARGET_BIN)/speedseq.config
	@echo "CNVNATOR_CHROMS_DIR=$(MKFILE_DIR)/$(ANNOTATIONS_DIR)/cnvnator_chroms" >> $(TARGET_BIN)/speedseq.config

	@echo "" >> $(TARGET_BIN)/speedseq.config
	@echo "# realign" >> $(TARGET_BIN)/speedseq.config
	@echo "BAMTOFASTQ=$(MKFILE_DIR)/$(TARGET_BIN)/bamtofastq.py" >> $(TARGET_BIN)/speedseq.config
	@echo "MBUFFER=$(MKFILE_DIR)/$(TARGET_BIN)/mbuffer" >> $(TARGET_BIN)/speedseq.config
	@echo "BAMHEADRG=$(MKFILE_DIR)/$(TARGET_BIN)/bamheadrg.py" >> $(TARGET_BIN)/speedseq.config
	@echo "BAMCLEANHEADER=$(MKFILE_DIR)/$(TARGET_BIN)/bamcleanheader.py" >> $(TARGET_BIN)/speedseq.config

# applications
bwa:
	$(MAKE) -C $(BWA_DIR)
	cp $(BWA_DIR)/bwa $(TARGET_BIN)

sambamba:
	cp $(SRC)/sambamba $(TARGET_BIN)

samblaster:
	$(MAKE) -C $(SAMBLASTER_DIR)
	cp $(SAMBLASTER_DIR)/samblaster $(TARGET_BIN)

freebayes:
	$(MAKE) -C $(FREEBAYES_DIR)
	cp $(FREEBAYES_DIR)/bin/freebayes $(TARGET_BIN)

lumpy:
	$(MAKE) -C $(LUMPY_DIR)
	cp $(LUMPY_DIR)/scripts/pairend_distro.py $(TARGET_BIN)
	cp $(LUMPY_DIR)/bin/lumpy $(TARGET_BIN)
	cp $(LUMPY_DIR)/bin/lumpyexpress $(TARGET_BIN)
	cp $(LUMPY_DIR)/scripts/vcfToBedpe $(TARGET_BIN)

svtyper:
	cp $(SVTYPER_DIR)/svtyper $(TARGET_BIN)

cnvnator-multi:
ifeq ($(ROOTSYS),)
	@echo -e  "\nWARNING: CNVnator not compiled because the ROOT package is not installed."
	@echo "Please see the README for instructions on manually installing ROOT."
else
	$(MAKE) -C $(CNVNATOR_DIR)
	cp $(CNVNATOR_DIR)/bin/cnvnator-multi $(TARGET_BIN)
	cp $(CNVNATOR_DIR)/bin/cnvnator_wrapper.py $(TARGET_BIN)
	cp $(CNVNATOR_DIR)/bin/cnvnator2VCF.pl $(TARGET_BIN)
	cp $(CNVNATOR_DIR)/bin/annotate_rd.py $(TARGET_BIN)
endif

tabix:
	$(MAKE) -C $(TABIX_DIR)
	cp $(TABIX_DIR)/tabix $(TARGET_BIN)
	cp $(TABIX_DIR)/bgzip $(TARGET_BIN)

vawk:
	cp $(VAWK_DIR)/vawk $(TARGET_BIN)

mbuffer:
	cd $(MBUFFER_DIR); ./configure --prefix=$(shell pwd)
	$(MAKE) -C $(MBUFFER_DIR)
	cp $(MBUFFER_DIR)/mbuffer $(TARGET_BIN)

parallel:
	cd $(PARALLEL_DIR); ./configure --prefix=$(shell pwd)
	$(MAKE) -C $(PARALLEL_DIR)
	cp $(PARALLEL_DIR)/src/parallel $(TARGET_BIN)

bamkit:
	cp $(BAMKIT_DIR)/bamtofastq.py $(TARGET_BIN)
	cp $(BAMKIT_DIR)/bamheadrg.py $(TARGET_BIN)
	cp $(BAMKIT_DIR)/bamgroupreads.py $(TARGET_BIN)
	cp $(BAMKIT_DIR)/bamfilterrg.py $(TARGET_BIN)
	cp $(BAMKIT_DIR)/bamcleanheader.py $(TARGET_BIN)
	cp $(BAMKIT_DIR)/bamlibs.py $(TARGET_BIN)

clean:
	rm -f \
		bin/bgzip \
		bin/sambamba \
		bin/cnvnator \
		bin/cnvnator2VCF.pl \
		bin/cnvnator_wrapper.py \
		bin/freebayes \
		bin/lumpy \
		bin/lumpyexpress \
		bin/pairend_distro.py \
		bin/samblaster \
		bin/svtyper \
		bin/tabix \
		bin/vawk \
		bin/vcfToBedpe \
		bin/bwa \
		bin/mbuffer \
		bin/parallel \
		bin/bamtofastq.py \
		bin/bamheadrg.py \
		bin/bamgroupreads.py \
		bin/bamfilterrg.py \
		bin/bamcleanheader.py \
		bin/bamlibs.py \
		bin/cnvnator-multi \
		bin/annotate_rd.py
	$(MAKE) -C $(BWA_DIR) clean
	$(MAKE) -C $(SAMBLASTER_DIR) clean
	$(MAKE) -C $(FREEBAYES_DIR) clean
	$(MAKE) -C $(LUMPY_DIR) clean
	$(MAKE) -C $(CNVNATOR_DIR) clean
	$(MAKE) -C $(TABIX_DIR) clean
	$(MAKE) -C $(MBUFFER_DIR) clean
	$(MAKE) -C $(PARALLEL_DIR) clean
