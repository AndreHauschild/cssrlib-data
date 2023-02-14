"""
 static test for PPP (IGS)
"""
import matplotlib.pyplot as plt
import numpy as np
from os.path import expanduser
import cssrlib.gnss as gn
import sys
from cssrlib.cssrlib import cssr
from cssrlib.gnss import ecef2pos, Nav, time2gpst, timediff
from cssrlib.pppigs import rtkinit, pppigspos
from cssrlib.rinex import rnxdec

l6file = '../data/2021078M.l6'
griddef = '../data/clas_grid.def'
"""
navfile = '../data/SEPT078M.21P'
obsfile = '../data/SEPT078M.21O'

# based on GSI F5 solution
xyz_ref = [-3962108.673,   3381309.574,   3668678.638]
pos_ref = ecef2pos(xyz_ref)
"""

navfile = '~/GNSS_NAV/IGS/2021/BRDC00IGS_R_20210780000_01D_MN.rnx'
#obsfile = '~/GNSS_OBS/VGS/HIGHRATE/2021/078/CHOF00JPN_R_20210781200_15M_01S_MO.rnx'
obsfile = '~/GNSS_OBS/IGS/HIGHRATE/2021/078/CHOF00JPN_S_20210781200_15M_01S_MO.rnx'
xyz_ref = [-3946217.2224, 3366689.3786, 3698971.7536]
pos_ref = ecef2pos(xyz_ref)

cs = cssr()
cs.monlevel = 1
cs.week = 2149
cs.read_griddef(griddef)

dec = rnxdec()
nav = Nav()
nav = dec.decode_nav(expanduser(navfile), nav)
nep = 720

t = np.zeros(nep)
tc = np.zeros(nep)
enu = np.ones((nep, 3))*np.nan
sol = np.zeros((nep, 4))
dop = np.zeros((nep, 4))
ztd = np.zeros((nep, 1))
smode = np.zeros(nep, dtype=int)

if dec.decode_obsh(expanduser(obsfile)) >= 0:

    rr = dec.pos
    rtkinit(nav, dec.pos)
    pos = ecef2pos(rr)
    inet = cs.find_grid_index(pos)

    fc = open(l6file, 'rb')
    if not fc:
        print("L6 messsage file cannot open.")
        sys.exit(-1)

    for ne in range(nep):

        obs = dec.decode_obs()
        week, tow = time2gpst(obs.t)

        cs.decode_l6msg(fc.read(250), 0)
        if cs.fcnt == 5:  # end of sub-frame
            cs.week = week
            cs.decode_cssr(cs.buff, 0)

        if ne == 0:
            t0 = nav.t = obs.t
            t0.time = t0.time//30*30
            cs.time = obs.t
            nav.time_p = t0

        cstat = cs.chk_stat()
        if cstat:
            pppigspos(nav, obs, cs)
            print(ne)

        t[ne] = timediff(nav.t, t0)
        tc[ne] = timediff(cs.time, t0)

        sol = nav.xa[0:3] if nav.smode == 4 else nav.x[0:3]
        enu[ne, :] = gn.ecef2enu(pos_ref, sol-xyz_ref)
        idx = nav.idx_ztd
        ztd[ne] = nav.xa[idx] if nav.smode == 4 else nav.x[idx]
        smode[ne] = nav.smode

    fc.close()
    dec.fobs.close()

fig_type = 1
ylim = 0.4

idx4 = np.where(smode == 4)[0]
idx5 = np.where(smode == 5)[0]
idx0 = np.where(smode == 0)[0]

fig = plt.figure(figsize=[7, 9])

if fig_type == 1:

    lbl_t = ['east[m]', 'north[m]', 'up[m]']
    for k in range(3):
        plt.subplot(4, 1, k+1)
        plt.plot(t[idx0], enu[idx0, k], 'r.')
        plt.plot(t[idx5], enu[idx5, k], 'y.')
        plt.plot(t[idx4], enu[idx4, k], 'g.')

        plt.xticks(np.arange(0, nep+1, step=30))
        plt.ylabel(lbl_t[k])
        plt.grid()
        #plt.axis([0, ne, -ylim, ylim])

    plt.subplot(4, 1, 4)
    plt.plot(t[idx0], ztd[idx0]*1e2, 'r.', markersize=8, label='none')
    plt.plot(t[idx5], ztd[idx5]*1e2, 'y.', markersize=8, label='float')
    plt.plot(t[idx4], ztd[idx4]*1e2, 'g.', markersize=8, label='fix')
    plt.xticks(np.arange(0, nep+1, step=30))
    plt.ylabel('ztd [cm]')
    plt.grid()
    plt.xlabel('time[s]')
    plt.legend()

elif fig_type == 2:

    ax = fig.add_subplot(111)

    #plt.plot(enu[idx0, 0], enu[idx0, 1], 'r.', label='stdpos')
    plt.plot(enu[idx5, 0], enu[idx5, 1], 'y.', label='float')
    plt.plot(enu[idx4, 0], enu[idx4, 1], 'g.', label='fix')

    plt.xlabel('easting [m]')
    plt.ylabel('northing [m]')
    plt.grid()
    plt.axis('equal')
    plt.legend()
    #ax.set(xlim=(-ylim, ylim), ylim=(-ylim, ylim))

plt.show()

if nav.fout is not None:
    nav.fout.close()
