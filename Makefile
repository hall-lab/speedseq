VERSION=0.0.1

TARGET_BIN=bin
SRC=src

SAMBLASTER_DIR=$(SRC)/samblaster
FREEBAYES_DIR=$(SRC)/freebayes
LUMPY_DIR=$(SRC)/lumpy-sv
SVTYPER_DIR=$(SRC)/svtyper
CNVNATOR_DIR=$(SRC)/cnvnator
TABIX_DIR=$(SRC)/tabix
VAWK_DIR=$(SRC)/vawk
SVTOOLS_DIR=$(SRC)/svtools

all:	samblaster freebayes lumpy svtyper tabix vawk svtools cnvnator

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

cnvnator:
	$(MAKE) -C $(CNVNATOR_DIR)
	cp $(CNVNATOR_DIR)/bin/cnvnator $(TARGET_BIN)
	cp $(CNVNATOR_DIR)/bin/cnvnator_wrapper.py $(TARGET_BIN)
	cp $(CNVNATOR_DIR)/bin/cnvnator2VCF.pl $(TARGET_BIN)

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

clean:
	rm -f bin/bedpeToBed12 bin/bedpeToVcf bin/bgzip bin/cnvnator bin/cnvnator2VCF.pl bin/cnvnator_wrapper.py bin/freebayes bin/lumpy bin/pairend_distro.py bin/samblaster bin/splitReadSamToBedpe bin/splitterToBreakpoint bin/svtyper bin/tabix bin/vawk bin/vcfToBedpe
	$(MAKE) -C $(SAMBLASTER_DIR) clean
	$(MAKE) -C $(FREEBAYES_DIR) clean
	$(MAKE) -C $(LUMPY_DIR) clean
	$(MAKE) -C $(CNVNATOR_DIR) clean
	$(MAKE) -C $(LUMPY_DIR) clean
	$(MAKE) -C $(TABIX_DIR) clean

