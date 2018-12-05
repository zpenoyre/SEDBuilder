# README #

`SEDBuilder` is an MIT-licensed Python code that collects VizieR and IRSA-tabulated photometric data points for any object with a 2MASS ID and uses them to make an SED. It contains SEDs for hundreds of example stars and YSOs in Taurus, Ophiuchus and Chamaeleon. It allows users to quickly access new photometric data, in the form of easy to use astropy tables, as well as to store and supplement that data.

# Installation #

To install SEDBuilder, download this directory, navigate to it, and run:

`python setup.py install`

# Dependencies #
* numpy
* pathlib
* astropy
* astroquery
(* matplotlib for the plotting function)

# Examples #

If you'd like to use the package, have a read through "introduction.ipynb". All the SED data published with Penoyre+2018 can also be found in "spectraData", the code used to generate it can be found in "dataCollation.ipynb". All saved data can be read with astropy tables or with the getSpectraFromFile() command.

# Contact #

We hope this will be a useful tool amongst a variety of acadmic communities. There is plenty of sope for it to grow and develop and if anyone should feel to make (or suggest) changes. You can contact me with any questions at zephyrpenoyre*at*gmail*dot*com
