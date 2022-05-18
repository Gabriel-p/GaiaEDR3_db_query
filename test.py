

import os
from astropy.io import ascii


# txt_f = '/home/gabriel/Descargas/all_files.txt'
# txt_f = 'Gaia_EDR3_files_FULL.txt'
# fdata = ascii.read(txt_f, format="no_header", delimiter=' ')
# all_names = []
# for name in fdata['col1']:
#     if not name.endswith('csv.gz'):
#         print(name)
#     if name in all_names:
#         print(name)
#     all_names.append(name)

txt_f = 'frame_ranges.txt'
fdata = ascii.read(txt_f)
all_names = []
for name in fdata['filename']:
    if not name.endswith('csv.gz'):
        print(name)
    if name in all_names:
        print(name)
    all_names.append(name)
breakpoint()

txt_f = 'Gaia_EDR3_files_FULL.txt'
fdata = ascii.read(txt_f, format="no_header")

downl_files = []
for filename in os.listdir("datafiles"):
    fname = filename.replace('_p', '')
    downl_files.append(fname)

# all_files = []
for file in fdata['col1']:
    name = file.split('/')[-1]
    # all_files.append(name)
    if name not in downl_files:
        print(name)