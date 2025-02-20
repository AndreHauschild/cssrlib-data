# CSSRlib-data

Sample scripts and dataset for CSSRlib
Pre-installation of CSSRlib is required.

## Prerequisites:

Additional python packages are required as prerequisites and can be installed via the following command

```
pip install -r requirements.txt
```

It is recommended to use a virtual environment to avoid version conflicts with the system-wide python installation.

## Installation in Virtual Environment

For this installation guide, it is assumed that the CSSRlib base repository has been cloned to a folder `cssrlib` and the samples repository has been cloned to a folder `cssrlib-data` in the current directory.

Create a virtual environment `cssrlib-venv` in the current directory with

```bash
python3 -m venv cssrlib-venv

```

The command `ls | grep cssrlib` should then return

```bash
cssrlib
cssrlib-data
cssrlib-venv
```

Active the environment and install the dependencies from the `requirements.txt` files of `cssrlib` and `cssrlib-data`:

```bash
source cssrlib-venv/bin/activate
pip install -r cssrlib/requirements.txt -r cssrlib-data/requirements.txt
```

Make sure the path to the `src` folder of the `cssrlib` base repository appears in the python path

```bash
echo $PYTHONPATH
```

If not, add it with the following export command, where `<path-to-cssrlib>` must be replaced with the full path to CSSRlib base repository.

```bash
export PYTHONPATH="$PYTHONPATH:<path-to-cssrlib>/src"
```

Check your installation

```bash
echo $PYTHONPATH
python -c "import cssrlib; print(cssrlib.__version__); print(cssrlib.__file__)"
```

## Handling Warnings from the PySolid module

The module `pysolid` is used for the computation of solid Earth tides. It contains a hard-coded leap second table with an expiration date, which is set to the next possible injection date of a leap second at the time of the last update. The table is frequently updated by the package maintainers. The following warning is issued when the expiration date is exceeded:

> Mild Warning -- time crossed leap second table boundaries.  Boundary edge value used instead

If you encounter this warning when executing CSSRlib scripts, it can most likely be removed by updating `pysolid` to the most recent version using

```bash
pip install --upgrade pysolid
```

If the warnings persist, installing the latest version from the repository can help

```bash
pip install -U git+https://github.com/insarlab/PySolid.git@main
```

## Ephemeris: RINEX/TLE

- test_eph.py reading/plotting ephemeris from RINEX 3
- test_tlesim.py reading/plotting TLE orbit

## Model

- test_trop.py tropospheric delay model
- cacode.py GPS/QZSS C/A code simulation

## Positioning: Static

- test_pntpos.py Standalone positioning
- test_rtk.py RTK positioning
- test_ppprtk.py PPP-RTK positioning

## Positioning: Kinematic

- test_pntpos2.py Standalone positioning
- test_rtk2.py RTK positioning
- test_ppprtk2.py PPP-RTK positioning

## ToDo-List

- [ ] Implement pole tide displacements
- [ ] Check and improve observation noise settings
- [ ] Add residual output
- [ ] Add check for observations, first two observations must be on different frequencies
- [ ] Number of frequencies `nav.nf` should be set automatically depending on specified signals
