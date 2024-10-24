"""
 static test for PPP (IGS)
"""
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib.dates as md
import numpy as np
from os.path import exists
from sys import stdout, exit

import cssrlib.gnss as gn
from cssrlib.gnss import ecef2pos, Nav
from cssrlib.gnss import time2doy, time2str, timediff, epoch2time
from cssrlib.gnss import rSigRnx
from cssrlib.gnss import sys2str
from cssrlib.peph import atxdec, searchpcv
from cssrlib.peph import peph, biasdec
from cssrlib.pppssr import pppos
from cssrlib.rinex import rnxdec

# Start epoch and number of epochs
#
dataset = 2  # 0: SEPT078M.21O, 1: SEPT1890.23O, 2: SEPT223Y.23O

if dataset == 0:  # SETP078M.21O
    ep = [2021, 3, 19, 12, 0, 0]
    let = 'M'
    xyz_ref = [-3962108.6617, 3381309.5232, 3668678.6410]
elif dataset == 1:  # SETP1890.23O
    ep = [2023, 7, 8, 4, 0, 0]
    let = '0'
    xyz_ref = [-3962108.7063, 3381309.5703, 3668678.6690]
elif dataset == 2:  # SETP223Z.23O
    ep = [2023, 8, 11, 21, 0, 0]
    let = 'Y'
    xyz_ref = [-3962108.7063, 3381309.5703, 3668678.6690]
else:
    print("ERROR: no RINEX data set selected!")
    exit(1)

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))

nep = 900*2

pos_ref = ecef2pos(xyz_ref)

if dataset == 2:
    navfile = '../data/doy{:03d}/BRD400DLR_S_{:04d}{:03d}0000_01D_MN.rnx'\
        .format(doy, year, doy)
    obsfile = '../data/doy{:03d}/SEPT{:03d}{}.{:02d}O'\
        .format(doy, doy, let, year % 2000)
else:
    navfile = '../data/SEPT{:03d}{}.{:02d}P'.format(doy, let, year % 2000)
    obsfile = '../data/SEPT{:03d}{}.{:02d}O'.format(doy, let, year % 2000)

# ac = 'COD0OPSFIN'
ac = 'COD0MGXFIN'

orbfile = '../data/{}_{:4d}{:03d}0000_01D_05M_ORB.SP3'\
    .format(ac, year, doy)

clkfile = '../data/{}_{:4d}{:03d}0000_01D_30S_CLK.CLK'\
    .format(ac, year, doy)

bsxfile = '../data/{}_{:4d}{:03d}0000_01D_01D_OSB.BIA'\
    .format(ac, year, doy)

if not exists(clkfile):
    orbfile = orbfile.replace('COD0OPSFIN', 'COD0OPSRAP')
    clkfile = clkfile.replace('COD0OPSFIN', 'COD0OPSRAP')
    bsxfile = bsxfile.replace('COD0OPSFIN', 'COD0OPSRAP')
if not exists(orbfile):
    orbfile = orbfile.replace('_05M_', '_15M_')

# Define signals to be processed
#
gnss = "GE"
sigs = []
if 'G' in gnss:
    sigs.extend([rSigRnx("GC1C"), rSigRnx("GC2W"),
                 rSigRnx("GL1C"), rSigRnx("GL2W"),
                 rSigRnx("GS1C"), rSigRnx("GS2W")])
if 'E' in gnss:
    sigs.extend([rSigRnx("EC1C"), rSigRnx("EC5Q"),
                 rSigRnx("EL1C"), rSigRnx("EL5Q"),
                 rSigRnx("ES1C"), rSigRnx("ES5Q")])
if 'C' in gnss:
    sigs.extend([rSigRnx("CC2I"), rSigRnx("CC6I"),
                 rSigRnx("CL2I"), rSigRnx("CL6I"),
                 rSigRnx("CS2I"), rSigRnx("CS6I")])

rnx = rnxdec()
rnx.setSignals(sigs)

nav = Nav()
orb = peph()

# Positioning mode
# 0:static, 1:kinematic
#
nav.pmode = 0

# Decode RINEX NAV data
#
nav = rnx.decode_nav(navfile, nav)

# Load precise orbits and clock offsets
#
nav = orb.parse_sp3(orbfile, nav)
nav = rnx.decode_clk(clkfile, nav)

# Load code and phase biases from Bias-SINEX
#
bsx = biasdec()
bsx.parse(bsxfile)

# Load ANTEX data for satellites and stations
#
if time > epoch2time([2022, 11, 27, 0, 0, 0]):
    atxfile = '../data/I20.ATX' if 'COD0MGXFIN' in ac else '../data/igs20.atx'
elif time > epoch2time([2021, 5, 2, 0, 0, 0]):
    atxfile = '../data/M20.ATX' if 'COD0MGXFIN' in ac else '../data/igs14.atx'
else:
    atxfile = '../data/M14.ATX' if 'COD0MGXFIN' in ac else '../data/igs14.atx'

atx = atxdec()
atx.readpcv(atxfile)

# Intialize data structures for results
#
t = np.zeros(nep)
enu = np.ones((nep, 3))*np.nan
sol = np.zeros((nep, 4))
ztd = np.zeros((nep, 1))
smode = np.zeros(nep, dtype=int)

# Logging level
#
nav.monlevel = 1  # TODO: enabled for testing!

# Load RINEX OBS file header
#
if rnx.decode_obsh(obsfile) >= 0:

    # Auto-substitute signals
    #
    rnx.autoSubstituteSignals()

    # Initialize position
    #
    ppp = pppos(nav, rnx.pos, 'test_pppigs.log')
    nav.ephopt = 4  # IGS
    nav.armode = 3

    nav.elmin = np.deg2rad(5.0)
    nav.thresar = 2.0

    # Get equipment information
    #
    nav.fout.write("FileName: {}\n".format(obsfile))
    nav.fout.write("Start   : {}\n".format(time2str(rnx.ts)))
    if rnx.te is not None:
        nav.fout.write("End     : {}\n".format(time2str(rnx.te)))
    nav.fout.write("Receiver: {}\n".format(rnx.rcv))
    nav.fout.write("Antenna : {}\n".format(rnx.ant))
    nav.fout.write("\n")

    if 'UNKNOWN' in rnx.ant or rnx.ant.strip() == "":
        nav.fout.write("ERROR: missing antenna type in RINEX OBS header!\n")

    # Set PCO/PCV information
    #
    nav.sat_ant = atx.pcvs
    nav.rcv_ant = searchpcv(atx.pcvr, rnx.ant,  rnx.ts)
    if nav.rcv_ant is None:
        nav.fout.write("ERROR: missing antenna type <{}> in ANTEX file!\n"
                       .format(rnx.ant))

    # Print available signals
    #
    nav.fout.write("Available signals\n")
    for sys, sigs in rnx.sig_map.items():
        txt = "{:7s} {}\n".format(sys2str(sys),
                                  ' '.join([sig.str() for sig in sigs.values()]))
        nav.fout.write(txt)
    nav.fout.write("\n")

    nav.fout.write("Selected signals\n")
    for sys, tmp in rnx.sig_tab.items():
        txt = "{:7s} ".format(sys2str(sys))
        for _, sigs in tmp.items():
            txt += "{} ".format(' '.join([sig.str() for sig in sigs]))
        nav.fout.write(txt+"\n")
    nav.fout.write("\n")

    # Skip epochs until start time
    #
    obs = rnx.decode_obs()
    while time > obs.t and obs.t.time != 0:
        obs = rnx.decode_obs()

    # Loop over number of epoch from file start
    #
    for ne in range(nep):

        # Set initial epoch
        #
        if ne == 0:
            nav.t = deepcopy(obs.t)
            t0 = deepcopy(obs.t)

        # Call PPP module with IGS products
        #
        ppp.process(obs, orb=orb, bsx=bsx)

        # Save output
        #
        t[ne] = timediff(nav.t, t0)/86400.0

        sol = nav.xa[0:3] if nav.smode == 4 else nav.x[0:3]
        enu[ne, :] = gn.ecef2enu(pos_ref, sol-xyz_ref)

        ztd[ne] = nav.xa[ppp.IT(nav.na)] \
            if nav.smode == 4 else nav.x[ppp.IT(nav.na)]
        smode[ne] = nav.smode

        nav.fout.write("{} {:14.4f} {:14.4f} {:14.4f} "
                       "ENU {:7.3f} {:7.3f} {:7.3f}, 2D {:6.3f}, mode {:1d}\n"
                       .format(time2str(obs.t),
                               sol[0], sol[1], sol[2],
                               enu[ne, 0], enu[ne, 1], enu[ne, 2],
                               np.sqrt(enu[ne, 0]**2+enu[ne, 1]**2),
                               smode[ne]))

        # Log to standard output
        #
        stdout.write('\r {} ENU {:7.3f} {:7.3f} {:7.3f}, 2D {:6.3f}, mode {:1d}'
                     .format(time2str(obs.t),
                             enu[ne, 0], enu[ne, 1], enu[ne, 2],
                             np.sqrt(enu[ne, 0]**2+enu[ne, 1]**2),
                             smode[ne]))

        # Get new epoch, exit after last epoch
        #
        obs = rnx.decode_obs()
        if obs.t.time == 0:
            break

    # Send line-break to stdout
    #
    stdout.write('\n')

    # Close RINEX observation file
    #
    rnx.fobs.close()

    # Close output file
    #
    if nav.fout is not None:
        nav.fout.close()

fig_type = 1
ylim = 1.0

idx4 = np.where(smode == 4)[0]
idx5 = np.where(smode == 5)[0]
idx0 = np.where(smode == 0)[0]

fig = plt.figure(figsize=[7, 9])
fig.set_rasterized(True)

fmt = '%H:%M'

if fig_type == 1:

    lbl_t = ['East [m]', 'North [m]', 'Up [m]']

    for k in range(3):
        plt.subplot(4, 1, k+1)
        plt.plot_date(t[idx0], enu[idx0, k], 'r.')
        plt.plot_date(t[idx5], enu[idx5, k], 'y.')
        plt.plot_date(t[idx4], enu[idx4, k], 'g.')

        plt.ylabel(lbl_t[k])
        plt.grid()
        plt.ylim([-ylim, ylim])
        plt.gca().xaxis.set_major_formatter(md.DateFormatter(fmt))

    plt.subplot(4, 1, 4)
    plt.plot_date(t[idx0], ztd[idx0]*1e2, 'r.', markersize=8, label='none')
    plt.plot_date(t[idx5], ztd[idx5]*1e2, 'y.', markersize=8, label='float')
    plt.plot_date(t[idx4], ztd[idx4]*1e2, 'g.', markersize=8, label='fix')
    plt.ylabel('ZTD [cm]')
    plt.grid()
    plt.gca().xaxis.set_major_formatter(md.DateFormatter(fmt))

    plt.xlabel('Time [HH:MM]')
    plt.legend()

elif fig_type == 2:

    ax = fig.add_subplot(111)

    # plt.plot(enu[idx0, 0], enu[idx0, 1], 'r.', label='stdpos')
    plt.plot(enu[idx5, 0], enu[idx5, 1], 'y.', label='float')
    plt.plot(enu[idx4, 0], enu[idx4, 1], 'g.', label='fix')

    plt.xlabel('Easting [m]')
    plt.ylabel('Northing [m]')
    plt.grid()
    plt.axis('equal')
    plt.legend()
    # ax.set(xlim=(-ylim, ylim), ylim=(-ylim, ylim))

plotFileFormat = 'eps'
plotFileName = '.'.join(('test_pppigs', plotFileFormat))

plt.savefig(plotFileName, format=plotFileFormat, bbox_inches='tight', dpi=300)
# plt.show()
