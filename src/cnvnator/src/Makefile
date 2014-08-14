VERSION	  = v0.3
#ROOTFLAGS = -pthread -m64
ROOTFLAGS = -m64 -O3
LIBS      = -lz
ROOTLIBS  = -L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d \
		-lGpad -lTree -lRint -lMatrix -lPhysics \
		-lMathCore -lThread -lGui

CXX    = g++ $(ROOTFLAGS) -DCNVNATOR_VERSION=\"$(VERSION)\"
SAMDIR = samtools
INC    = -I$(ROOTSYS)/include -I$(SAMDIR)
SAMLIB = $(SAMDIR)/libbam.a

OBJDIR = obj
OBJS   = $(OBJDIR)/cnvnator.o  \
	 $(OBJDIR)/HisMaker.o  \
	 $(OBJDIR)/AliParser.o \
	 $(OBJDIR)/Genotyper.o \
	 $(OBJDIR)/Interval.o  \
	 $(OBJDIR)/Genome.o \
	 $(OBJDIR)/TreeBuilder.o

DISTRIBUTION = $(PWD)/CNVnator_$(VERSION).zip
TMPDIR	     =  /tmp
CNVDIR	     = CNVnator_$(VERSION)
MAINDIR	     = $(TMPDIR)/$(CNVDIR)
SRCDIR	     = $(MAINDIR)/src

all: cnvnator

cnvnator: $(OBJS)
	$(CXX) -o $@ $(OBJS) $(SAMLIB) $(LIBS) $(ROOTLIBS)

$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(OBJDIR)
	$(CXX) $(INC) -c $< -o $@

clean:
	rm -f $(OBJS)

distribution: clean all
	@echo Creating directory ...
	@rm -rf $(MAINDIR)
	@rm -f  $(DISTRIBUTION)
	@mkdir  $(MAINDIR)
	@mkdir  $(SRCDIR)
	@echo Copying files ...
	@cp *.hh *.cpp  $(SRCDIR)
	@cp Makefile    $(SRCDIR)
	@cp -r samtools $(SRCDIR)
	@rm -f $(SRCDIR)/samtools/samtools
	@rm -f $(SRCDIR)/samtools/*.o
	@rm -f $(SRCDIR)/samtools/*/*.o
	@rm -f $(SRCDIR)/samtools/*/*/*.o
	@rm -f $(SRCDIR)/samtools/*.a
	@rm -f $(SRCDIR)/samtools/*/*.a
	@rm -f $(SRCDIR)/samtools/*/*/*.a
	@rm -fr $(SRCDIR)/samtools/.svn
	@rm -fr $(SRCDIR)/samtools/*/.svn
	@rm -fr $(SRCDIR)/samtools/*/*/.svn
	@cp README          $(MAINDIR)
	@cp CITATION        $(MAINDIR)
	@cp license.rtf     $(MAINDIR)
	@cp cnvnator2VCF.pl $(MAINDIR)
	@echo Zipping ...
	@ln -s $(MAINDIR)
	@zip -qr $(DISTRIBUTION) $(CNVDIR)
	@rm $(CNVDIR)
