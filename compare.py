

from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt


"""
 -- output format : csv
 SELECT "I/350/gaiaedr3".RA_ICRS,  "I/350/gaiaedr3".DE_ICRS,  "I/350/gaiaedr3".Source,  "I/350/gaiaedr3".Plx,  "I/350/gaiaedr3".pmRA,  "I/350/gaiaedr3".pmDE,  "I/350/gaiaedr3".Gmag, 
 "I/350/gaiaedr3"."BP-RP"
 FROM "I/350/gaiaedr3"
 WHERE 1=CONTAINS(POINT('ICRS',"I/350/gaiaedr3".RA_ICRS,"I/350/gaiaedr3".DE_ICRS), BOX('ICRS', 048.2525, +17.5252, 1., 1.)) 
AND "I/350/gaiaedr3".Gmag<20
"""

data = ascii.read('cluster.csv')
vizier = ascii.read('/home/gabriel/Descargas/result.csv')

plt.scatter(data['RA_ICRS'], data['DE_ICRS'], marker='x', alpha=.3)
plt.scatter(vizier['RA_ICRS'], vizier['DE_ICRS'], alpha=.3)
plt.show()
