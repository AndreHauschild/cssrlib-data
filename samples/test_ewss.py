"""
Emergency Warning Satellite Service (EWSS) sample
"""
from binascii import unhexlify
import numpy as np
import bitstruct as bs
from cssrlib.gnss import time2doy, epoch2time, gpst2time
from cssrlib.ewss import jmaDec, camfDec
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
import matplotlib.pyplot as plt

# QZSS DCR MT43 (212bits)
# QZSS DCX MT44 L-Alert, J-Alert, International (CAMF) (212bits)

# Select test case
#
dataset = 1

if dataset == 1:

    ep = [2025, 8, 21, 7, 0, 0]
    file_sbas = '../data/doy2025-233/233h_sbas.txt'
    prn_ref = 199
    sbas_type = 1  # L1: 0, L5: 1
    nf = 2

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))

nep = 300

cs = jmaDec()
cs.monlevel = 0
cs.time = time

csx = camfDec()
csx.monlevel = 2
csx.time = time

flg_dcr = True  # True if DCR is enabled
flg_dcx = True  # True if DCX is enabled
flg_plt = True

pos = None

params = []

if True:
    if 'sbas' in file_sbas:  # SIS
        dtype = [('wn', 'int'), ('tow', 'float'), ('prn', 'int'),
                 ('type', 'int'), ('marker', 'S2'), ('nav', 'S124')]
        v = np.genfromtxt(file_sbas, dtype=dtype)

    tow = np.unique(v['tow'])
    # Loop over number of epoch from file start
    #
    for ne in range(nep):
        vi = v[(v['tow'] == tow[ne]) & (v['prn'] == prn_ref)]

        vi_dcr = vi[vi['type'] == 43]  # DCR

        if flg_dcr and len(vi_dcr) > 0:
            for vi_ in vi_dcr:
                cs.time = gpst2time(vi_['wn'], vi_['tow'])
                buff = unhexlify(vi_['nav'])
                cs.decode(buff, 14)
                print(cs.gen_msg(cs.dc))

            if cs.dc == 10:  # Typhoon
                pos = cs.pos if pos is None else np.vstack([pos, cs.pos])
                params.append(cs.params)

        vi_dcx = vi[vi['type'] == 44]  # DCX

        if flg_dcx and len(vi_dcx) > 0:
            for vi_ in vi_dcx:
                csx.time = gpst2time(vi_['wn'], vi_['tow'])
                buff = unhexlify(vi_['nav'])

                sdmt, sdm = bs.unpack_from('u1u9', buff, 14)
                # print(f"sdmt={sdmt} sdm={sdm:0b}")
                csx.decode(buff, 24)  # decode CAMF message
                csx.decode_ext(buff, 24+122)  # decode extended message


if flg_plt:
    extent = [128, 134, 30, 35]  # kyusyu, Japan

    fig = plt.figure()

    stamen_terrain = cimgt.GoogleTiles()
    ax = plt.axes(projection=stamen_terrain.crs)
    ax.set_extent(extent)
    ax.add_image(stamen_terrain, 8)

    if pos is not None:
        sm = plt.plot(pos[:, 1], pos[:, 0], 'r.', transform=ccrs.PlateCarree())
    plt.show()
