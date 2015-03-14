## Stellar Instrumental Magnitudes + Filter Functions -> Calibration ##

#### <font color='blue'> Email big.macs.calibrate@gmail.com for updates/notifications </font> ####

An implementation of algorithms described in Kelly et al. 2014 (http://arxiv.org/abs/1208.0602)
Updated 8/25/12

#### Author List ####
**P. L. Kelly**, A. von der Linden, D. Applegate, M. Allen, S. W. Allen, P. R. Burchat, D. L. Burke, H. Ebeling, P. Capak, O. Czoske, D.  Donovan, A. Mantz, and R. G. Morris

#### 1-2% Color Calibration ####

Big MACS is a Python program that estimates an accurate photometric calibration from only an input catalog of stellar
magnitudes and filter transmission functions. The user does not have to measure color terms which can be difficult to characterize.

  * Supplied with filter transmission functions, Big MACS synthesizes an expected stellar locus for your data and then simultaneously solves for all unknown zeropoints when fitting to the instrumental locus.

  * Uses a spectroscopic model for the SDSS stellar locus in color-color space and filter functions to compute expected locus.

  * The stellar locus model is corrected for Milky Way reddening.

If SDSS or 2MASS photometry is available for stars in field, code can yield a highly accurate absolute calibration.

#### Magnitude calibration with 2MASS ####
Automatically match against 2MASS catalog (useful if 2MASS stars are not saturated in science exposures)

#### Incorporate SDSS photometry ####
Automatically match against SDSS catalog

#### Excellent for photometric redshifts ####
See Figure below and accompanying paper

#### Advantages ####
  * Calculates color transformation from SDSS to instrumental magnitudes automatically
  * Simultaneously solves for all zeropoint simultaneously using new, simplified algorithm
  * SDSS stellar locus is dereddened for Milky Way dust extinction
  * Improved zeropoint accuracy


#### Checkout with SVN or Download (by going to downloads section at top) ####
```
svn checkout http://big-macs-calibrate.googlecode.com/svn/ big-macs-calibrate
```

#### Why Big MACS? ####
Many of the galaxy clusters in our `Weighing the Giants' study
are members of the MACS sample and they are extremely massive.


![http://farm9.staticflickr.com/8283/7613582508_cbccece26e.jpg](http://farm9.staticflickr.com/8283/7613582508_cbccece26e.jpg)

Color calibration of photometry through nine Subaru SuprimeCam and CFHT MegaPrime filter bands of the RXJ1347-11 galaxy cluster field. We simultaneously vary eight zeropoints and hold one constant during the fit to maximize agreement between the instrumental stellar locus and the model stellar locus.

![http://farm8.staticflickr.com/7275/7723411818_8340373f95.jpg](http://farm8.staticflickr.com/7275/7723411818_8340373f95.jpg)

Photometric redshifts computed after calibration with Big MACS software.

#### Notes ####

Make sure to exclude stars with saturated or non-linear measurements.

The spectroscopic model for the locus is only trustworthy to approximately
10500 Angstroms (the red limit of the SDSS z' filter), so synthesized
magnitudes from near-IR filter functions are not likely to be correct.
Therefore, near-IR observations need to be calibrated separately,
although some work is ongoing to extend to templates to redder wavelengths.

The magnitudes, if any, that you hold "FIXED" need to be corrected for MW reddening
along the line of sight.