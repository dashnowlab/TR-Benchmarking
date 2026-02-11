import os
import itertools 
import Levenshtein as lv
from pathlib import Path

root_dir = Path(__file__).resolve().parent.parent
sys.path.append(str(root_dir))

from helpers.readers import *

# set directory variables for file i/o
PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA_DIR = os.path.join(PROJ_ROOT, 'HG001.PAW79146')

VCF = os.path.join(DATA_DIR, 'HG001.PAW79146.haplotagged.URfix.straglr.vcf')
vr = VCFReader(VCF)

with vr as vr:
    vr.skipMetaData()
    diffs = 0

    while not vr.end_state:
        vr.read()
        good = False
        len_ind = 0
        for i, strs in enumerate(vr.info):
            if "SVLEN" in strs:
                good = True
                len_ind = i

        if good:
            svl = int(vr.info[len_ind].strip("SVLEN=").split(',')[0])
            calc_len = vr.pos_info["end"] - vr.pos_info["start"]
            if svl != calc_len:
                print(f"{vr.pos_info} calc len: {calc_len}, svlen: {svl}")
                diffs += 1

    print(diffs)



