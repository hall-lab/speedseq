VERSION=0.0.1

TARGET_BIN=bin
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
SVTOOLS_DIR=$(SRC)/svtools
MBUFFER_DIR=$(SRC)/mbuffer
SCRIPTS_DIR=$(SRC)/scripts

all:	bwa sambamba samblaster freebayes lumpy svtyper tabix vawk svtools mbuffer scripts cnvnator-multi

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

svtools:
	cp $(SVTOOLS_DIR)/bedpeToBed12 $(TARGET_BIN)
	cp $(SVTOOLS_DIR)/bedpeToVcf $(TARGET_BIN)
	cp $(SVTOOLS_DIR)/splitReadSamToBedpe $(TARGET_BIN)
	cp $(SVTOOLS_DIR)/splitterToBreakpoint $(TARGET_BIN)
	cp $(SVTOOLS_DIR)/vcfToBedpe $(TARGET_BIN)
	cp $(SVTOOLS_DIR)/lumpyToBedpe $(TARGET_BIN)

mbuffer:
	cd $(MBUFFER_DIR); ./configure --prefix=$(shell pwd)
	$(MAKE) -C $(MBUFFER_DIR)
	cp $(MBUFFER_DIR)/mbuffer $(TARGET_BIN)

scripts:
	cp $(SCRIPTS_DIR)/bamtofastq.py $(TARGET_BIN)
	cp $(SCRIPTS_DIR)/bamheadrg.py $(TARGET_BIN)
	cp $(SCRIPTS_DIR)/bamgroupreads.py $(TARGET_BIN)
	cp $(SCRIPTS_DIR)/bamfilterrg.py $(TARGET_BIN)
	cp $(SCRIPTS_DIR)/bamcheck.py $(TARGET_BIN)

clean:
	rm -f \
		bin/bedpeToBed12 \
		bin/bedpeToVcf \
		bin/bgzip \
		bin/cnvnator \
		bin/cnvnator2VCF.pl \
		bin/cnvnator_wrapper.py \
		bin/freebayes \
		bin/lumpy \
		bin/pairend_distro.py \
		bin/samblaster \
		bin/splitReadSamToBedpe \
		bin/splitterToBreakpoint \
		bin/svtyper \
		bin/tabix \
		bin/vawk \
		bin/vcfToBedpe \
		bin/bwa \
		bin/lumpyToBedpe \
		bin/mbuffer \
		bin/bamtofastq.py \
		bin/bamheadrg.py \
		bin/bamgroupreads.py \
		bin/bamfilterrg.py \
		bin/bamcheck.py \
		bin/cnvnator-multi \
		bin/annotate_rd.py
	$(MAKE) -C $(BWA_DIR) clean
	$(MAKE) -C $(SAMBLASTER_DIR) clean
	$(MAKE) -C $(FREEBAYES_DIR) clean
	$(MAKE) -C $(LUMPY_DIR) clean
	$(MAKE) -C $(CNVNATOR_DIR) clean
	$(MAKE) -C $(TABIX_DIR) clean
	$(MAKE) -C $(MBUFFER_DIR) clean
