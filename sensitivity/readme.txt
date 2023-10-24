WFI Sensitivity Calculator
__________________________

-->This file describes the input parameters of the sensitivity calculator. 

The following options will be presented in order when running the code - it will use the default values if nothing is inputted:

1) Filters in which the sensitivity is calculated can be chosen from two instruments, or custom filters made by the user. Choose 'moircs' or 'swims' to use existing moircs or swims filters, or choose 'custom' to specify your own custom filters.

2) For either instrument, you can select between 'bb', 'mb' and 'nb' for broad, medium or narrow-band filters.

3) If you chose 'custom' for the first option, you can make your own custom filters by providing four parameters (separated by tab or comma) in a file named 'custom.asc'. The parameters are filter name, central wavelength, filter width and filter transmission in order. Please see the included custom.asc file for formatting - you can edit this file or make your own with the same name and place it in the working directory. Wavelengths should be given in microns (ideally between 0.9 - 2.5 um) and the transmission as a fraction (i.e. 80% is 0.8).

4) AO correction can be turned on by selecting 'glao' or have no AO correction by selecting 'noao'.

5) Pixel scale in arcsec/pixel. Default value is 0.1 (planned pixel scale for ULTIMATE).

6) Read noise in e-/rms. Default value is 16 (expected read noise for the detector).

7) Exposure time should be given in hours (or fraction of hours). The default value is 1 hour.

8) The desired S/N ratio. The default value is 5.

9) Combined instrument+filter throughput (between 0 to 1). Default value is 0.45.

10) The water vapor and airmass to select the background model. Default values are 1.0 mm and 1.0 respectively.

The code will create a plot of Limiting magnitude vs Wavelength for every chosen filter under three different seeing conditions (good, moderate and bad). A table will also be generated with the aforementioned values, with the sensitivities for both glao and noao included.
