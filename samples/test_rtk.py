"""
 static test for RTK
"""

import matplotlib.pyplot as plt
import numpy as np
from sys import exit as sys_exit
from sys import stdout

from cssrlib.rinex import rnxdec, sync_obs
from cssrlib.gnss import rSigRnx, time2str, timediff, Nav, ecef2pos, ecef2enu
from cssrlib.peph import atxdec, searchpcv
from cssrlib.rtk import rtkpos
from cssrlib.plot import plot_enu

ngsantfile = '../data/GSI_PCV.TXT'

nav = Nav()

dataset = 3

if dataset == 1:
    bdir = '../data/doy2021-078/'
    navfile = bdir+'SEPT078M.21P'
    obsfile = bdir+'SEPT078M.21O'
    basefile = bdir+'3034078M.21O'
    xyz_ref = [-3962108.673, 3381309.574, 3668678.638]
    nav.rb = [-3959400.631, 3385704.533, 3667523.111]  # GSI 3034 fujisawa
    atxfile = '../data/antex/igs14.atx'
elif dataset == 2:
    bdir = '../data/doy2023-238/'
    navfile = bdir+'SEPT238A.23P'
    obsfile = bdir+'SEPT238A.23O'
    basefile = bdir+'3034238A.23O'
    xyz_ref = [-3962108.7007, 3381309.5532, 3668678.6648]
    nav.rb = [-3959400.6443, 3385704.4948, 3667523.1275]  # GSI 3034 fujisawa
    atxfile = '../data/antex/igs20.atx'

elif dataset == 3:
    bdir = '../data/doy2025-233/'
    navfile = bdir+'233h_rnx.nav'
    obsfile = bdir+'233h_rnx.obs'
    basefile = bdir+'3034233h.25o'
    xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]
    nav.rb = [-3959400.6242, 3385704.4927, 3667523.1257]  # GSI 3034 fujisawa
    atxfile = '../data/antex/igs20.atx'

pos_ref = ecef2pos(xyz_ref)

# Define signals to be processed
#
gnss = "GEJ"  # "GEJ"
sigs = []
sigsb = []
if 'G' in gnss:
    sigs.extend([rSigRnx("GC1C"), rSigRnx("GC2W"),
                 rSigRnx("GL1C"), rSigRnx("GL2W"),
                 rSigRnx("GS1C"), rSigRnx("GS2W")])

if 'E' in gnss:
    sigs.extend([rSigRnx("EC1C"), rSigRnx("EC5Q"),
                 rSigRnx("EL1C"), rSigRnx("EL5Q"),
                 rSigRnx("ES1C"), rSigRnx("ES5Q")])

if 'J' in gnss:
    sigs.extend([rSigRnx("JC1C"), rSigRnx("JC2L"),
                 rSigRnx("JL1C"), rSigRnx("JL2L"),
                 rSigRnx("JS1C"), rSigRnx("JS2L")])

# rover
dec = rnxdec()
dec.setSignals(sigs)

dec.decode_nav(navfile, nav)

# base
decb = rnxdec()
decb.setSignals(sigs)

decb.decode_obsh(basefile)
dec.decode_obsh(obsfile)

decb.autoSubstituteSignals()
dec.autoSubstituteSignals()

nep = 300

t = np.zeros(nep)
enu = np.zeros((nep, 3))
smode = np.zeros(nep, dtype=int)

rtk = rtkpos(nav, dec.pos, 'test_rtk.log')
rr = dec.pos

# Load ANTEX data for satellites and stations
#
atx = atxdec()
atx.readpcv(atxfile)
atx.readngspcv(ngsantfile)

# Set PCO/PCV information
#
nav.rcv_ant = searchpcv(atx.pcvr, dec.ant,  dec.ts)
if nav.rcv_ant is None:
    print("ERROR: missing antenna type <{}> in ANTEX file!".format(dec.ant))
    sys_exit(-1)
nav.rcv_ant_b = searchpcv(atx.pcvr, decb.ant,  dec.ts)
if nav.rcv_ant_b is None:
    print("ERROR: missing antenna type <{}> in ANTEX file!".format(decb.ant))
    sys_exit(-1)

# nav.excl_sat = [20]
# nav.cnr_min_gpy = 20

# Get equipment information
#
print("Rover:")
print("  Receiver:", dec.rcv)
print("  Antenna :", dec.ant)
print()
print("Base:")
print("  Receiver:", decb.rcv)
print("  Antenna :", decb.ant)
print()

for ne in range(nep):
    obs, obsb = sync_obs(dec, decb)
    if ne == 0:
        t0 = nav.t = obs.t

    rtk.process(obs, obsb=obsb)
    t[ne] = timediff(nav.t, t0)/86400.0
    sol = nav.xa[0:3] if nav.smode == 4 else nav.x[0:3]
    enu[ne, :] = ecef2enu(pos_ref, sol-xyz_ref)
    smode[ne] = nav.smode

    nav.fout.write("{} {:14.4f} {:14.4f} {:14.4f} "
                   "ENU {:7.4f} {:7.4f} {:7.4f}, 2D {:6.4f}, mode {:1d}\n"
                   .format(time2str(obs.t),
                           sol[0], sol[1], sol[2],
                           enu[ne, 0], enu[ne, 1], enu[ne, 2],
                           np.sqrt(enu[ne, 0]**2+enu[ne, 1]**2),
                           smode[ne]))

    # Log to standard output
    #
    stdout.write('\r {} ENU {:7.4f} {:7.4f} {:7.4f}, 2D {:6.4f}, mode {:1d}'
                 .format(time2str(obs.t),
                         enu[ne, 0], enu[ne, 1], enu[ne, 2],
                         np.sqrt(enu[ne, 0]**2+enu[ne, 1]**2),
                         smode[ne]))

# Send line-break to stdout
#
stdout.write('\n')

dec.fobs.close()
decb.fobs.close()

fig_type = 1

if fig_type == 1:
    plot_enu(t, enu, smode)
elif fig_type == 2:
    plot_enu(t, enu, smode, figtype=fig_type)

plotFileFormat = 'eps'
plotFileName = '.'.join(('test_rtk', plotFileFormat))

plt.savefig(plotFileName, format=plotFileFormat, bbox_inches='tight', dpi=300)
# plt.show()
