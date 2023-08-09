# Build a HAWKI ETC

## Exercise:

*Kieran wants to study the initial mass function of an open cluster in the Large Magellanic Cloud.
He plans to use the HAWK-I instrument on the VLT and would like to know what the faintest stars are that will be detectable (S/N > 5) with 1 hour of observing time in the K filter.
He assumes average (50th percentile) observing conditions will apply for the observation.*

## Deliverable:
The goal of this exercise is to build an exposure time calculator (ETC) module in Python for the [HAWK-I imager instrument on the VLT](https://www.eso.org/sci/facilities/paranal/instruments/hawki.html).

## Implementation:

1. Ensure you have Python >=3.7 installed.
2. Install the required packages:

```bash
pip install HowManyPhotons numpy astropy matplotlib
```

3. To run the ETC module and obtain the detectable star brightness:
```bash
python etc_module_hmbp.py
```

This will print the brightness of the faintest detectable star (with S/N > 5).

4. To run the test suite:
```bash
python test_etc.py
```

This will run the following tests:
- `test_sky_noise`: Tests the sky noise calculation for the Ks filter.
- `test_dark_current_noise`: Tests the dark current noise calculation.
- `test_read_noise`: Tests the read noise calculation.
- `test_compute_sn`: Tests the signal-to-noise ratio computation for a given star brightness.
- `test_detectable_star_brightness`: Tests the detectable star brightness calculation for a desired SNR.
- `test_plot_function`: Tests the plotting function for SNR as a function of star brightness.


