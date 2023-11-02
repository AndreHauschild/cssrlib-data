#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 static test for PPP (Kepler)
"""
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib.dates as md
import numpy as np
import os
from os.path import basename, expanduser, splitext
from sys import stdout

import cssrlib.gnss as gn
from cssrlib.gnss import ecef2pos, enu2xyz, Nav
from cssrlib.gnss import time2doy, time2str, timediff, epoch2time
from cssrlib.gnss import rSigRnx
from cssrlib.gnss import sys2str
from cssrlib.gnss import uTropoModel
from cssrlib.peph import atxdec, searchpcv
from cssrlib.peph import peph, biasdec
from cssrlib.pppssr import pppos
from cssrlib.rinex import rnxdec


# Start epoch and number of epochs
#
ep = [2010, 3, 21, 0, 0, 0]
xyz_ref = [4186548.9344, 835107.0779, 4723754.0722]
ecc_ref = [0.0, 0.0, 0.0460]

pos_ref = ecef2pos(xyz_ref)
xyz_ref = xyz_ref + enu2xyz(pos_ref)@ecc_ref

pos_ref = ecef2pos(xyz_ref)

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))

nep = int(4*3600/30)

#tropo = 'Niell'
tropo = 'Hopfi'

bdir = expanduser('~/Projects/PPP_Kepler/SIM_OBER')

obsfile = bdir+'/_obs/EPOS_OBSDATA_{}_OBER_100321.rnx.rnxman'.format(tropo)
logfile = splitext(basename(obsfile))[0]+'.log'
pltfile = logfile.replace('.log', '.eps')

# orbit and biases
#
orbfile = bdir+'/_sp3/a0_sim_MEO_OBER_Hopfi.sp3.sp3man'
bsxfile = bdir+'/_bia/HARDWARE_DELAYS_{}_OBER_100321.bsx'.format(tropo)

# Define signals to be processed
#
gnss = "E"
sigs = []
if 'E' in gnss:
    sigs.extend([rSigRnx("EC1C"), rSigRnx("EC5Q"),
                 rSigRnx("EL1C"), rSigRnx("EL5Q"),
                 rSigRnx("ES1C"), rSigRnx("ES5Q")])

rnx = rnxdec()
rnx.setSignals(sigs)

nav = Nav()
orb = peph()

# Positioning mode
# 0:static, 1:kinematic
#
nav.pmode = 0

# Load precise orbits and clock offsets
#
nav = orb.parse_sp3(orbfile, nav)

# Load code and phase biases from Bias-SINEX
#
bsx = biasdec()
bsx.parse(bsxfile)

# Load ANTEX data for satellites and stations
#
atxfile = bdir+'/_atx/PHASCENT1_igs08.atx_LEO_GAL_6COSMIC_sim'

atx = atxdec()
atx.readpcv(atxfile)

# Initialize data structures for results
#
t = np.zeros(nep)
enu = np.ones((nep, 3))*np.nan
sol = np.zeros((nep, 4))
ztd = np.zeros((nep, 1))
smode = np.zeros(nep, dtype=int)

resc = np.ones((nep, gn.uGNSS.MAXSAT, nav.nf))*np.nan
resp = np.ones((nep, gn.uGNSS.MAXSAT, nav.nf))*np.nan

# Logging level
#
nav.monlevel = 0  # TODO: enabled for testing!

# Load RINEX OBS file header
#
if rnx.decode_obsh(obsfile) >= 0:

    # Auto-substitute signals
    #
    rnx.autoSubstituteSignals()

    # Initialize position
    #
    ppp = pppos(nav, rnx.pos, logfile)
    nav.ephopt = 4  # SP3

    nav.armode = 0
    nav.thresar = 2.0

    nav.elmin = np.deg2rad(7.5)
    nav.cnr_min = 0

    # Match process noise of 20cm/2mm for EPOS simulations
    #
    nav.err = [0, 0.000, 0.002]       # [m] sigma

    # Zero process noise for fixed-position model
    #
    #nav.q[0:3] = 0.0

    # Use Hopfield tropo model
    #
    nav.trpModel = uTropoModel.HOPF

    # Do not use receiver antenna corrections
    #
    nav.useRxPco = False

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

        resc[ne, :, :] = nav.resc
        resp[ne, :, :] = nav.resp

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


# Solution index
#
idx4 = np.where(smode == 4)[0]
idx5 = np.where(smode == 5)[0]
idxA = np.where(smode != 0)[0]
idx0 = np.where(smode == 0)[0]

# Residual statistics
#
print()
print("Residual statistics")
print()
for f in range(nav.nf):
    print("Code  {:1d} {:5.2f} +/- {:5.2f} cm"
          .format(f,
                  np.nanmean(resc[idxA, :, f])*1e2,
                  np.nanstd(resc[idxA, :, f])*1e2))
print()
for f in range(nav.nf):
    print("Phase {:1d} {:5.2f} +/- {:5.2f} cm"
          .format(f,
                  np.nanmean(resp[idxA, :, f])*1e2,
                  np.nanstd(resp[idxA, :, f])*1e2))
print()

# Plot results
#
ylim = 1.0

fig = plt.figure(figsize=[7, 9])
fig.set_rasterized(True)

fmt = '%H:%M'

lbl_t = ['East [m]', 'North [m]', 'Up [m]']

for k in range(3):
    plt.subplot(6, 1, k+1)
    plt.plot_date(t[idx5], enu[idx5, k], 'y.')
    plt.plot_date(t[idx4], enu[idx4, k], 'g.')

    plt.ylabel(lbl_t[k])
    plt.grid()
    plt.ylim([-ylim, ylim])
    plt.gca().xaxis.set_major_formatter(md.DateFormatter(fmt))

plt.subplot(6, 1, 4)
plt.plot_date(t[idx5], ztd[idx5]*1e2, 'y.', markersize=8, label='float')
plt.plot_date(t[idx4], ztd[idx4]*1e2, 'g.', markersize=8, label='fix')
plt.ylabel('ZTD [cm]')
plt.grid()
plt.gca().xaxis.set_major_formatter(md.DateFormatter(fmt))
plt.legend()

plt.subplot(6, 1, 5)
for f in range(nav.nf):
    plt.plot_date(t[idx5], resc[idx5, :, f]*1e2, 'y.', markersize=8)
    plt.plot_date(t[idx4], resc[idx4, :, f]*1e2, 'g.', markersize=8)
plt.ylabel('Code [cm]')
plt.grid()
plt.gca().xaxis.set_major_formatter(md.DateFormatter(fmt))

plt.subplot(6, 1, 6)
for f in range(nav.nf):
    plt.plot_date(t[idx5], resp[idx5, :, f]*1e2, 'y.', markersize=8)
    plt.plot_date(t[idx4], resp[idx4, :, f]*1e2, 'g.', markersize=8)
plt.ylabel('Phase [cm]')
plt.grid()
plt.gca().xaxis.set_major_formatter(md.DateFormatter(fmt))

plt.xlabel('Time [HH:MM]')

plotFileFormat = splitext(pltfile)[1][1:]
plt.savefig(pltfile, format=plotFileFormat, bbox_inches='tight', dpi=300)

# Call the plotting script
#
os.system("{}/plot_pppkepler.py {}".format(os.path.dirname(__file__), logfile))
