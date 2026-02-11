import os
import sys
from contextlib import ExitStack
from pathlib import Path

root_dir = Path(__file__).resolve().parent.parent
sys.path.append(str(root_dir))

from helpers.readers import BEDReader
from helpers.utils import setupVCFReader, getFileName
from helpers.constants import *

# WARNING: LINES SAVED INTO MEMORY WHEN USING FULL ORDERING
# USE CAUTION WHEN RUNNING FULL ORDERING ON LARGE FILES(i.e. setting chrom_only = False)


def sortVCF(vcf_rdr, order, stk, chrom_only = True):
    return_loc = vcf_rdr.cur_loc

    # vcf reader will already error out if the file is not valid, so no need for extensive checks
    if vcf_rdr.path.endswith(".gz"):
        sof_name = os.path.join(DATA_DIR, vcf_rdr.path[:-7] + ".sorted.vcf")
    else:
        sof_name = os.path.join(DATA_DIR, vcf_rdr.path[:-4] + ".sorted.vcf")

    sof = stk.enter_context(open(sof_name, "w")) 

    # copy header to new file
    vcf_rdr.read() 
    while vcf_rdr.raw_line.startswith("#") and not vcf_rdr.end_state:
        sof.write(vcf_rdr.raw_line)
        vcf_rdr.read() 
    # save the file position of the end of the header/metadata
    vcf_rdr.header_end = vcf_rdr.file_obj.tell() 
    

    print(f"Sorting {getFileName(vcf_rdr.path)}")

    for chrom in order:
        chrom_data_lines = []
        while not vcf_rdr.end_state:
            if vcf_rdr.chrom == chrom:
                if chrom_only:
                    sof.write(vcf_rdr.raw_line)
                else:
                    if chrom_data_lines:
                        i = 0
                        while vcf_rdr.pos > chrom_data_lines[i][0]:
                            i += 1
                            if i == len(chrom_data_lines):
                                break
                
                        chrom_data_lines.insert(i, [vcf_rdr.pos, vcf_rdr.raw_line])
                        
                    else:
                        chrom_data_lines.append([vcf_rdr.pos, vcf_rdr.raw_line])
         
            vcf_rdr.read()
        if not chrom_only:
            for lines in chrom_data_lines:
                sof.write(lines[1])
        vcf_rdr._setFilePosition(vcf_rdr.header_end)
        vcf_rdr.end_state = False
        vcf_rdr.read()

    vcf_rdr._setFilePosition(return_loc) # return vcf rdr to original location
    print(f"File sorted.")


def getBEDOrder(bed):
    chrom_order = []
    while not bed.end_state:
        if bed.chrom not in chrom_order:
            chrom_order.append(bed.chrom)
        bed.read()

    return chrom_order


def grabFromDir(file_type: str, dir_name: str):
    path_list = []
    with os.scandir(dir_name) as dir_entries:
        for entry in dir_entries:
            if entry.is_file():
                if entry.path.endswith('.' + file_type) or entry.path.endswith("." + file_type + '.gz'):
                    vcf_list.append(entry.path)

    return path_list


# loads all vcf file paths from a given directory into a list
vcf_list = grabFromDir("vcf", os.path.join(DATA_DIR, "HG001.30x"))

bed_file = "BED_files\\benchmark-catalog-v2.vamos.bed"

# vcf_list = [ # the file name, and the offset amount (eg. [-1, 0] for 1 based inclusive), and bool for whether is is position only
#     #os.path.join(DATA_DIR, "HG001.PAW79146.haplotagged.URfix.atarva.vcf"),
#     os.path.join(DATA_DIR,"HG001.PAW79146.haplotagged.strdust.vcf"),
#     #os.path.join(DATA_DIR,"HG001.PAW79146.haplotagged.URfix.strkit.vcf"), 
#     os.path.join(DATA_DIR,"HG001.PAW79146.haplotagged.URfix.longTR.vcf.gz"),
#     #os.path.join(DATA_DIR,"HG001.PAW79146.haplotagged.URfix.straglr.vcf"),
#     #os.path.join(DATA_DIR,"HG001.PAW79146.haplotagged.URfix.vamos.vcf"),
#     os.path.join(DATA_DIR,"medaka_to_ref.TR.vcf"),
#     #"test_cases.vcf"
#     ]

vcf_rdrs = []
with ExitStack() as stack: 
    # create bed reader and enter the file in stack
    bed = stack.enter_context(BEDReader(os.path.join(DATA_DIR, bed_file)))
    bed.read()

    order = getBEDOrder(bed)
    print(order)

    # create list of vcf reader objects and their new sorted, output file
    for i, vcf in enumerate(vcf_list):
        vcf_rdrs.append(setupVCFReader(vcf=vcf, 
                                        skip_head=False,
                                        stk=stack))
        
        sortVCF(vcf_rdrs[i], order, stack)