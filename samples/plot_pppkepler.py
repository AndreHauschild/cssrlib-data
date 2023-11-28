#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

  -------------
  | Changelog |
  -------------

    2023/08/24  AHA  Initial version

"""

import argparse
from datetime import datetime
from glob import glob
import matplotlib.pyplot as plt
import matplotlib.dates as md
import matplotlib.ticker as tck
from os.path import expanduser
import numpy as np


def toEpoch(field):

    try:
        return datetime.strptime(field, '%Y-%m-%d %H:%M:%S')
    except:
        return None


# Format time axis
#
def timeAxis(ax):

    # ax.xaxis.set_major_locator(md.HourLocator(interval=3))
    # ax.xaxis.set_minor_locator(md.HourLocator())
    ax.xaxis.set_major_formatter(md.DateFormatter('%H:%M'))

    return ax

# Define RMS errors
#


def rms(v):
    return np.sqrt(np.mean(v**2))

# Main function
#


def main():

    # Parse command line arguments
    #
    parser = argparse.ArgumentParser(description="GNSS data plotter")

    # Input file
    #
    parser.add_argument('inpFileName',
                        help='Input files (wildcard allowed)')
    parser.add_argument("-p", "--plotFileName", default=None,
                        help="Plot file name")
    parser.add_argument("-f", "--plotFileFormat", help="Plot file format",
                        default='eps')

    # Retrieve all command line arguments
    #
    args = parser.parse_args()

    data = {}

    # Get solution files
    #
    inpFileNames = glob(expanduser(args.inpFileName))
    for inpFileName in sorted(inpFileNames):

        # Print header line
        #
        print("Processing {}".format(inpFileName))

        # Intialize data structures for results
        #
        t0 = None
        t = []
        enu = []
        smode = []
        nep = 0

        # Open input files
        #
        with open(inpFileName, 'r') as f:

            # Parse file line by line
            #
            for line in f.read().splitlines():

                if 'mode' not in line:
                    continue

                fields = line.replace(',', '').split()

                epoch = toEpoch(fields[0]+' '+fields[1])
                if t0 is None:
                    t0 = epoch

                err_e = float(fields[6])
                err_n = float(fields[7])
                err_u = float(fields[8])
                err_2d = float(fields[10])
                mode = int(fields[12])

                t.append((epoch-t0).seconds/86400)
                enu.append([err_e, err_n, err_u, err_2d])
                smode.append(mode)
                nep += 1
                #print(epoch, err_e, err_n, err_u, err_2d, mode)

        # Convert to numpy arrays for indexing with np.where()
        #
        t = np.array(t).reshape((nep, 1))
        enu = np.array(enu).reshape((nep, 4))

        idx4 = np.where(np.array(smode, dtype=int) == 4)[0]
        idx5 = np.where(np.array(smode, dtype=int) == 5)[0]
        #idx0 = np.where(np.array(smode, dtype=int) == 0)[0]
        idx45 = np.where(np.array(smode, dtype=int) != 0)[0]

        # Statistics
        #
        rms2Dfloat = np.sqrt(np.mean(enu[idx5, 3]**2))
        if len(idx4) > 0:
            rms2Dfixed = np.sqrt(np.mean(enu[idx4, 3]**2))
            rmsUpfixed = np.sqrt(np.mean(enu[idx4, 2]**2))

        print()
        print("Solution statistics:")
        print()

        print("RMS float 2D solution {:4.1f} cm".format(rms2Dfloat*1e2))
        if len(idx4) > 0:
            print("RMS fixed 2D solution {:4.1f} cm".format(rms2Dfixed*1e2))
            print("RMS fixed up solution {:4.1f} cm".format(rmsUpfixed*1e2))
        print()

        # Find last index with solution over conversion limit
        #
        limConv = 0.1
        epoConv = np.where(enu[idx5, 3] > limConv)[0][-1]

        print("Converged solution")
        print()
        print("Time until convergence {:4.1f} min ({:3d} epochs)"
              .format(t[epoConv][0]*24*60, epoConv))
        print()
        for i, s in enumerate(("East", "North", "Up", "2D")):
            print("{:5s} {:4.1f}+/-{:2.1f} cm (RMS {:4.1f} cm)"
                  .format(s,
                          np.mean(enu[epoConv:, i])*1e2,
                          np.std(enu[epoConv:, i])*1e2,
                          rms(enu[epoConv:, i])*1e2))
        print()

        # Convergence time statistics
        #
        site = inpFileName.split('_')[-1][0:4]
        data.update({site: t[epoConv][0]*24*60})

        # Plot
        #
        fig = plt.figure(figsize=[9, 6])
        fig.set_rasterized(True)

        fmt = '%H:%M'

        ylim = 0.20
        tmaj = 0.05
        tmin = 0.01
        legend = ['East', 'North', 'Up']

        plt.subplot(2, 1, 1)
        plt.plot_date(t[idx45], enu[idx45, 0:3], markersize=4)

        plt.legend(legend)
        plt.ylabel("Position Error [m]")
        plt.grid()
        plt.ylim([-ylim, ylim])
        plt.gca().xaxis.set_major_formatter(md.DateFormatter(fmt))
        plt.gca().yaxis.set_major_locator(tck.MultipleLocator(base=tmaj))
        plt.gca().yaxis.set_minor_locator(tck.MultipleLocator(base=tmin))

        plt.subplot(2, 1, 2)
        #plt.plot_date(t[idx0], enu[idx0, 3], 'r.', markersize=4, label='None')
        if len(t[idx5]) > 0:
            plt.plot_date(t[idx5], enu[idx5, 3], 'y.',
                          markersize=4, label='Float')
        if len(t[idx4]) > 0:
            plt.plot_date(t[idx4], enu[idx4, 3], 'g.',
                          markersize=4, label='Fixed')
        if len(t[idx4]) > 0 and len(t[idx5]) > 0:
            plt.legend()

        #ylim = ylim/2
        #tmaj = tmaj/2
        #tmin = tmin/2

        plt.ylabel("Position Error 2D [m]")
        plt.grid()
        plt.ylim([0, ylim])
        plt.gca().xaxis.set_major_formatter(md.DateFormatter(fmt))
        plt.gca().yaxis.set_major_locator(tck.MultipleLocator(base=tmaj))
        plt.gca().yaxis.set_minor_locator(tck.MultipleLocator(base=tmin))

        plt.xlabel('Time [HH:MM]')

        plotFileFormat = 'eps'
        plotFileName = inpFileName.replace('.log', '_2D')
        plotFileName = '.'.join((plotFileName, plotFileFormat))

        plt.savefig(plotFileName, format=plotFileFormat,
                    bbox_inches='tight', dpi=300)
        # plt.show()

    # Overall statistics
    #
    print()
    print("Convergence time statistics:")
    print()

    values = list(data.values())
    mean = np.mean(values)
    sigma = np.std(values)
    p95 = np.percentile(values, 95)

    def sort_by_time(item):
        return item[1]

    sorted_data = dict(sorted(data.items(), key=sort_by_time))

    for site, dt in sorted_data.items():
        print("{:4s}  {:5.2f} min".format(site, dt))
    print()
    print("mean  {:5.2f} +/- {:5.1f} min (95% {:5.2f} min)"
          .format(mean, sigma, p95))


# Main program
#
if __name__ == "__main__":
    main()
