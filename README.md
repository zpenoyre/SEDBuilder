# README #

`spectracles` is an MIT-licensed Python code that collects VizieR-tabulated photometric data points for any object with a 2MASS ID and uses them to make an SED.

I'm immediately sorry for the name.

# Installation #

To install spectracles, download this directory, navigate to it, and run:

`python setup.py install`

# Dependencies #
* numpy
* pathlib
* astropy
* astroquery

# Examples #

If you'd like to use the package, have a read through "introduction.ipynb". All the SED data published with Penoyre+2018 can also be found in "spectraData", and read with astropy tables or with the getSpectraFromFile() command.
