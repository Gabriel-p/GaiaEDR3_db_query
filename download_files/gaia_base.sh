#!/bin/bash

source ~/.bashrc

cd GaiaEDR3/

conda activate py3

python Gaia_EDR3_download.py > out_gaia.txt
