#!/usr/bin/env python

import tempfile
import argparse, sys, os
import math, time, re
from argparse import RawTextHelpFormatter
from subprocess import Popen, PIPE, STDOUT

__author__ = "Colby Chiang (cc2qe@virginia.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2014-04-28 14:31 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
svgt\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: Compute genotype of structural variants based on breakpoint depth")
    parser.add_argument('-v', '--input_vcf', type=argparse.FileType('r'), default=None, help='VCF input')
    parser.add_argument('-r', '--root', type=str, required=True, help='CNVnator .root histogram file (required)')
    parser.add_argument('-w', '--window', type=int, required=True, help='CNVnator window size (required)')
    parser.add_argument('-s', '--sample', type=str, required=True, help='sample to annotate')
    parser.add_argument('--cnvnator', type=str, default='cnvnator-multi', help='path to cnvnator-multi binary')
    parser.add_argument('-o', '--output_vcf', type=argparse.FileType('w'), default=sys.stdout, help='output VCF to write (default: stdout)')
    parser.add_argument('--debug', action='store_true', help='debugging verbosity')

    # parse the arguments
    args = parser.parse_args()

    # if no input, check if part of pipe and if so, read stdin.
    if args.input_vcf == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        else:
            args.input_vcf = sys.stdin
    # send back the user input
    return args

class Vcf(object):
    def __init__(self):
        self.file_format = 'VCFv4.2'
        #self.fasta = fasta
        # self.reference = ''
        self.sample_list = []
        self.info_list = []
        self.format_list = []
        self.alt_list = []
        self.add_format('GT', 1, 'String', 'Genotype')

    def add_header(self, header):
        for line in header:
            if line.split('=')[0] == '##fileformat':
                self.file_format = line.rstrip().split('=')[1]
            elif line.split('=')[0] == '##reference':
                self.reference = line.rstrip().split('=')[1]
            elif line.split('=')[0] == '##INFO':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_info(*[b.split('=')[1] for b in r.findall(a)])
            elif line.split('=')[0] == '##ALT':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_alt(*[b.split('=')[1] for b in r.findall(a)])
            elif line.split('=')[0] == '##FORMAT':
                a = line[line.find('<')+1:line.find('>')]
                r = re.compile(r'(?:[^,\"]|\"[^\"]*\")+')
                self.add_format(*[b.split('=')[1] for b in r.findall(a)])
            elif line[0] == '#' and line[1] != '#':
                self.sample_list = line.rstrip().split('\t')[9:]

    # return the VCF header
    def get_header(self):
        header = '\n'.join(['##fileformat=' + self.file_format,
                            '##fileDate=' + time.strftime('%Y%m%d')] + \
                               [i.hstring for i in self.info_list] + \
                               [a.hstring for a in self.alt_list] + \
                               [f.hstring for f in self.format_list] + \
                               ['\t'.join([
                                            '#CHROM',
                                            'POS',
                                            'ID',
                                            'REF',
                                            'ALT',
                                            'QUAL',
                                            'FILTER',
                                            'INFO',
                                            'FORMAT'] + \
                                              self.sample_list
                                          )])
        return header

    def add_info(self, id, number, type, desc):
        if id not in [i.id for i in self.info_list]:
            inf = self.Info(id, number, type, desc)
            self.info_list.append(inf)

    def add_alt(self, id, desc):
        if id not in [a.id for a in self.alt_list]:
            alt = self.Alt(id, desc)
            self.alt_list.append(alt)

    def add_format(self, id, number, type, desc):
        if id not in [f.id for f in self.format_list]:
            fmt = self.Format(id, number, type, desc)
            self.format_list.append(fmt)

    def add_sample(self, name):
        self.sample_list.append(name)

    class Info(object):
        def __init__(self, id, number, type, desc):
            self.id = str(id)
            self.number = str(number)
            self.type = str(type)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##INFO=<ID=' + self.id + ',Number=' + self.number + ',Type=' + self.type + ',Description=\"' + self.desc + '\">'

    class Alt(object):
        def __init__(self, id, desc):
            self.id = str(id)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##ALT=<ID=' + self.id + ',Description=\"' + self.desc + '\">'

    class Format(object):
        def __init__(self, id, number, type, desc):
            self.id = str(id)
            self.number = str(number)
            self.type = str(type)
            self.desc = str(desc)
            # strip the double quotes around the string if present
            if self.desc.startswith('"') and self.desc.endswith('"'):
                self.desc = self.desc[1:-1]
            self.hstring = '##FORMAT=<ID=' + self.id + ',Number=' + self.number + ',Type=' + self.type + ',Description=\"' + self.desc + '\">'

class Variant(object):
    def __init__(self, var_list, vcf):
        self.chrom = var_list[0]
        self.pos = int(var_list[1])
        self.var_id = var_list[2]
        self.ref = var_list[3]
        self.alt = var_list[4]
        self.qual = float(var_list[5])
        self.filter = var_list[6]
        self.sample_list = vcf.sample_list
        self.info_list = vcf.info_list
        self.info = dict()
        self.format_list = vcf.format_list
        self.active_formats = list()
        self.gts = dict()
        # make a genotype for each sample at variant
        for i in xrange(len(self.sample_list)):
            s_gt = var_list[9+i].split(':')[0]
            s = self.sample_list[i]
            self.gts[s] = Genotype(self, s, s_gt)
        # import the existing fmt fields
        for i in xrange(len(self.sample_list)):
            s = self.sample_list[i]
            for j in zip(var_list[8].split(':'), var_list[9+i].split(':')):
                self.gts[s].set_format(j[0], j[1])

        self.info = dict()
        i_split = [a.split('=') for a in var_list[7].split(';')] # temp list of split info column
        for i in i_split:
            if len(i) == 1:
                i.append(True)
            self.info[i[0]] = i[1]

    def set_info(self, field, value):
        if field in [i.id for i in self.info_list]:
            self.info[field] = value
        else:
            sys.stderr.write('\nError: invalid INFO field, \"' + field + '\"\n')
            exit(1)

    def get_info(self, field):
        return self.info[field]

    def get_info_string(self):
        i_list = list()
        for info_field in self.info_list:
            if info_field.id in self.info.keys():
                if info_field.type == 'Flag':
                    i_list.append(info_field.id)
                else:
                    i_list.append('%s=%s' % (info_field.id, self.info[info_field.id]))
        return ';'.join(i_list)

    def get_format_string(self):
        f_list = list()
        for f in self.format_list:
            if f.id in self.active_formats:
                f_list.append(f.id)
        return ':'.join(f_list)

    def genotype(self, sample_name):
        if sample_name in self.sample_list:
            return self.gts[sample_name]
        else:
            sys.stderr.write('\nError: invalid sample name, \"' + sample_name + '\"\n')

    def get_var_string(self):
        s = '\t'.join(map(str,[
            self.chrom,
            self.pos,
            self.var_id,
            self.ref,
            self.alt,
            '%0.2f' % self.qual,
            self.filter,
            self.get_info_string(),
            self.get_format_string(),
            '\t'.join(self.genotype(s).get_gt_string() for s in self.sample_list)
        ]))
        return s

class Genotype(object):
    def __init__(self, variant, sample_name, gt):
        self.format = dict()
        self.variant = variant
        self.set_format('GT', gt)

    def set_format(self, field, value):
        if field in [i.id for i in self.variant.format_list]:
            self.format[field] = value
            if field not in self.variant.active_formats:
                self.variant.active_formats.append(field)
                # sort it to be in the same order as the format_list in header
                self.variant.active_formats.sort(key=lambda x: [f.id for f in self.variant.format_list].index(x))
        else:
            sys.stderr.write('\nError: invalid FORMAT field, \"' + field + '\"\n')
            exit(1)

    def get_format(self, field):
        return self.format[field]

    def get_gt_string(self):
        g_list = list()
        for f in self.variant.active_formats:
            if f in self.format:
                if type(self.format[f]) == float:
                    g_list.append('%0.2f' % self.format[f])
                else:
                    g_list.append(self.format[f])
            else:
                g_list.append('.')
        return ':'.join(map(str,g_list))


# primary function
def sv_readdepth(vcf_file, sample, root, window, vcf_out, debug, cnvnator_path):
    in_header = True
    header = []
    vcf = Vcf()

    # make 2 passes through the vcf file, because it's much faster to give a single
    # command to CNVnator.
    coord_list = tempfile.NamedTemporaryFile(dir='./', delete=False)
    for line in vcf_file:
        if in_header:
            if line[0] == '#':
                header.append(line) 
                if line[1] != '#':
                    vcf_samples = line.rstrip().split('\t')[9:]
                continue
            else:
                in_header = False
                vcf.add_header(header)

        v = line.rstrip().split('\t')
        var = Variant(v, vcf)

        if var.get_info('SVTYPE') != 'BND':
            end = var.get_info('END')
            # if the variant is negative size, swap it.
            if int(end) < int(var.pos):
                coord_list.write('%s:%s-%s\n' % (var.chrom, end, var.pos))
            else:
                coord_list.write('%s:%s-%s\n' % (var.chrom, var.pos, end))

    coord_list.write('exit\n')
    coord_list.close()

    p1 = Popen(['cat', coord_list.name], stdout=PIPE)
    cmd = map(str, [cnvnator_path, '-root', root, '-genotype', window])
    p2 = Popen(cmd, stdin=p1.stdout, stdout=PIPE)
    p3 = Popen(['awk', '{ if($1!="Assuming"){print $4} }'], stdin=p2.stdout, stdout=PIPE)
    cn_list = map(float, p3.communicate()[0].split('\n')[:-1])

    # go through the VCF a second time and add the read depth annotations
    vcf_file.seek(0)
    in_header = True
    header = []
    vcf = Vcf()
    i = 0
    for line in vcf_file:
        if in_header:
            if line[0] == '#':
                header.append(line) 
                if line[1] != '#':
                    vcf_samples = line.rstrip().split('\t')[9:]
                continue
            else:
                in_header = False
                vcf.add_header(header)
                vcf.add_format('CN', 1, 'Float', 'Copy number of structural variant segment.')
                vcf_out.write(vcf.get_header() + '\n')

        v = line.rstrip().split('\t')
        var = Variant(v, vcf)
        if var.get_info('SVTYPE') != 'BND':
            var.genotype(sample).set_format('CN', cn_list[i])
            i += 1

        # write the VCF
        vcf_out.write(var.get_var_string() + '\n')

    vcf_out.close()

    # remove temp file
    os.unlink(coord_list.name)
    
    return

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # call primary function
    sv_readdepth(args.input_vcf, args.sample, args.root, args.window, args.output_vcf, args.debug, args.cnvnator)

    # close the files
    args.input_vcf.close()

# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
