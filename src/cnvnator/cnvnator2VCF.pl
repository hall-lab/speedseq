#!/usr/bin/perl

use warnings;
use strict;

my $usage = "\tcnvnator2VCF.pl [-prefix prefix] file.calls [genome_dir]\n";

my ($file,$dir,$prefix) = ("","","");
my $n_arg = scalar(@ARGV);
foreach (my $i = 0;$i < $n_arg;$i++) {
    my $arg = $ARGV[$i];
    if ($arg eq "-prefix") {
	if ($i + 1 < $n_arg) { $prefix = $ARGV[++$i]; }
	else                 { die "Not enough arguments.\n"; }
    }
    if    (length($file) <= 0) { $file = $arg; }
    elsif (length($dir)  <= 0) { $dir  = $arg; }
    else {
	print STDERR "Too many arguments.\n";
	last;
    }
}

if (length($file) <= 0) {
    print STDERR $usage,"\n";
    exit;
}

open(FILE,$file) or die "Can't open file ",$file,".\n";
print STDERR "Reading calls ...\n";
my ($pop_id) = split(/\./,$file);
print '##fileformat=VCFv4.0',"\n";
print '##fileDate='.`date '+%Y%m%d'`;
print '##reference=1000GenomesPhase3_decoy-GRCh37',"\n";
print '##source=CNVnator',"\n";
print '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',"\n";
print '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',"\n";
print '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',"\n";
print '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',"\n";
print '##INFO=<ID=natorRD,Number=1,Type=Float,Description="Normalized RD">',"\n";
print '##INFO=<ID=natorP1,Number=1,Type=Float,Description="e-val by t-test">',"\n";
print '##INFO=<ID=natorP2,Number=1,Type=Float,Description="e-val by Gaussian tail">',"\n";
print '##INFO=<ID=natorP3,Number=1,Type=Float,Description="e-val by t-test (middle)">',"\n";
print '##INFO=<ID=natorP4,Number=1,Type=Float,Description="e-val by Gaussian tail (middle)">',"\n";
print '##INFO=<ID=natorQ0,Number=1,Type=Float,Description="Fraction of reads with 0 mapping quality">',"\n";
print '##INFO=<ID=natorPE,Number=1,Type=Integer,Description="Number of paired-ends support the event">',"\n";
print '##INFO=<ID=SAMPLES,Number=.,Type=String,Description="Sample genotyped to have the variant">',"\n";
print '##ALT=<ID=DEL,Description="Deletion">',"\n";
print '##ALT=<ID=DUP,Description="Duplication">',"\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
my ($prev_chrom,$chrom_seq,$count) = ("","",0);
while (my $line = <FILE>) {
    my ($type,$coor,$len,$rd,$p1,$p2,$p3,$p4,$q0,$pe) = split(/\s+/,$line);
    my ($chrom,$start,$end) = split(/[\:\-]/,$coor);
    my $isDel = ($type eq "deletion");
    my $isDup = ($type eq "duplication");
    if ($isDup) {
    } elsif ($isDel) {
    } else {
	print STDERR "Skipping unrecognized event type '",$type,"'.\n";
	next;
    }
    if ($chrom ne $prev_chrom) {
	$chrom_seq  = parseChrom($chrom,$dir);
	$prev_chrom = $chrom;
    }
    $count++;
    my $id = "";
    if (length($prefix) > 0) { $id = $prefix."_"; }
    $id .= "CNVnator_";
    if    ($isDel) { $id .= "del_"; }
    elsif ($isDup) { $id .= "dup_"; }
    $id .= $count;
    print $chrom,"\t",$start,"\t",$id,"\t";
    if ($start < length($chrom_seq)) {
	print substr($chrom_seq,$start - 1,1),"\t";
    } else {
	print "N\t";
    }
    if    ($isDel) { print "<DEL>"; }
    elsif ($isDup) { print "<DUP>"; }
    print "\t.\tPASS\t";
    my $INFO = "END=".$end;
    if    ($isDel) {
	$INFO .= ";SVTYPE=DEL";
	$INFO .= ";SVLEN=-".int($len);
    } elsif ($isDup) {
	$INFO .= ";SVTYPE=DUP";
	$INFO .= ";SVLEN=".int($len);
    }
    $INFO   .= ";IMPRECISE";
    if (defined($rd)) { $INFO .= ";natorRD=".$rd; }
    if (defined($p1)) { $INFO .= ";natorP1=".$p1; }
    if (defined($p2)) { $INFO .= ";natorP2=".$p2; }
    if (defined($p3)) { $INFO .= ";natorP3=".$p3; }
    if (defined($p4)) { $INFO .= ";natorP4=".$p4; }
    if (defined($q0)) { $INFO .= ";natorQ0=".$q0; }
    if (defined($pe)) { $INFO .= ";natorPE=".$pe; }
    print $INFO."\n";
}
close(FILE);

exit;

sub parseChrom
{
    my ($chrom,$dir) = @_;
    my $file = $dir."/".$chrom.".fa";
    if (!open(SEQFILE,"$file")) {
        print STDERR "Can't parse sequence for chromosome ",$chrom,".\n";
        return "";
    }
    my ($header,$seq) = ("","");
    while (my $line = <SEQFILE>) {
        chomp($line);
        if (length($line) <= 0) { next; }
        if (substr($line,0,1) eq ">") {
            $header = substr($line,1);
            $seq = "";
        } else { $seq .= $line; }
    }
    close(SEQFILE);
    return $seq;
}
