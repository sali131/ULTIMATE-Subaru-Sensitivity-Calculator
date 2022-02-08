WFI Sensitivity Calculator
__________________________

-->This file describes the input parameters of the sensitivity calculator. 

The following options will be presented in order when running the code:

1) Filters in which the sensitivity is calculated can be chosen from a set of pre-existing filters, or custom filters made by the user. Choose 'moircs' or 'swims' to use existing moircs or swims filters, or choose 'custom' to specify your own custom filters.

2) If you chose 'moircs', you can select between 'bb' and 'nb' for broad-band or narrow-band filters. For 'swims', you can choose between 'bb', 'mb' or 'nb' - same as moircs but also with medium-band filters as an option.

3) If you chose 'custom', you can make your own custom filters by providing four parameters (separated by tab or comma) in a file named custom.asc that the code will read. The parameters are filter name, central wavelength, filter width and filter transmission in order. Please see the included custom.asc file for formatting - you can edit this file or make your own with the same name and place it in the working directory. Wavelengths should be given in microns (ideally between 0.9 - 2.5 um) and the transmission as a fraction (i.e. 80% is 0.8).

4) You can turn AO correction on by selecting 'glao' or have no AO correction by selecting 'noao'.

5) Exposure time should be given in hours (or fraction of hours). The default value is 5 hours if nothing is inputted.

6) The desired signal-to-noise ratio can also be given. The default value is 5 if nothing is inputted.

7) The code will create a plot of Limiting magnitude vs Wavelength for every chosen filter under three different seeing conditions (good, moderate and bad). A table will also be generated with the aforementioned values, with the sensitivities for both glao and noao included.