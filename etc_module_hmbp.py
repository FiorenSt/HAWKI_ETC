import numpy as np
from astropy import units as u
import hmbp
import matplotlib.pyplot as plt
import io
import sys


class HAWKI_ETC:
    """
    Exposure Time Calculator (ETC) for the HAWKI imager instrument on the VLT.
    This class provides methods to compute various noise components, signal-to-noise ratio,
    and plots the SNR as a function of brightness.
    """
    def __init__(self, exposure_time=3600 * u.s, airmass=2.0, pwv=5.0, seeing=0.8 * u.arcsec,
                 obstruction_factor=1, reflectance_factor=1):
        """
        Initializes the HAWKI_ETC with given parameters or defaults.
        :param exposure_time: Exposure time in seconds.
        :param airmass: Airmass value.
        :param pwv: Precipitable Water Vapor.
        :param seeing: Seeing in arcseconds.
        """
        self.exposure_time = exposure_time
        self.airmass = airmass
        self.pwv = pwv
        self.seeing = seeing
        # Telescope parameters
        self.telescope_size = np.pi * (4 * u.m) ** 2
        self.obstruction_factor = obstruction_factor # set to 1 due to no better feeling on it for HAWK
        self.reflectance_factor = reflectance_factor # set to 1
        self.efficiency_factor = self.obstruction_factor * self.reflectance_factor
        self.collection_area = self.telescope_size * self.efficiency_factor
        self.quantum_efficiency = 0.9 * (u.electron / u.ph) # approximation based on plots in the user manual
        self.pixel_scale = 0.106 * u.arcsec / u.pixel
        self.aper = np.pi * self.seeing ** 2.0 # simplified, just based on seeing (0.8'' as given by the calculator)

    def sky_noise(self, filter_name):
        """
        Calculates the sky noise for a given filter based on background, quantum efficiency, exposure time,
        collection area, and aperture.

        :param filter_name: Name of the filter.
        :return: Sky noise value.
        """
        # Redirect stdout and stderr to suppress unwanted warnings from the `hmbp.in_skycalc_background` function.
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = io.StringIO()
        sys.stderr = io.StringIO()

        try:
            sky = hmbp.in_skycalc_background(filter_name, airmass=self.airmass, pwv=self.pwv) / (u.arcsec ** 2)
        finally:
            # Reset stdout and stderr to their original values after capturing the output
            sys.stdout = old_stdout
            sys.stderr = old_stderr

        sky_aper = sky * self.quantum_efficiency * self.exposure_time * self.collection_area * self.aper
        return sky_aper

    def dark_current_noise(self):
        """
        Computes the dark current noise based on the aperture, pixel scale, and exposure time.
        :return: Dark current noise value.
        """
        n_pix = self.aper / self.pixel_scale ** 2.0
        dark_current = 0.0125 * u.electron / u.second / u.pixel ** 2 # between 0.10 and 0.15 in the manual
        dark_aper = dark_current * n_pix * self.exposure_time
        return dark_aper

    def read_noise(self):
        """
        Calculates the squared read noise over the aperture.
        :return: Squared read noise value in electron^2.
        """
        n_pix = self.aper / self.pixel_scale ** 2.0
        read_noise_per_pixel = 8.5 * u.electron / u.pixel  # between 5 and 12 in the manual
        squared_read_noise_aper = (read_noise_per_pixel ** 2) * n_pix
        return squared_read_noise_aper

    def compute_sn(self, star_flux, filter_name="Ks"):
        """
        Computes the signal-to-noise ratio (SNR) for a given star flux and filter name.
        :param star_flux: Star flux value.
        :param filter_name: Name of the filter.
        :return: SNR value.
        """
        photons = hmbp.for_flux_in_filter(filter_name, star_flux, instrument="HAWKI", observatory="Paranal")
        signal = photons * self.collection_area * self.quantum_efficiency * self.exposure_time
        total_noise = np.sqrt(
            signal.value
            + self.sky_noise(filter_name).value
            + self.dark_current_noise().value
            + self.read_noise().value
        )
        snr_aper = signal.value / total_noise
        return snr_aper

    def flux_to_mag(self, flux, zero_point_flux_photon):
        """
        Convert flux to magnitude using the zero-point flux.
        :param flux: The flux value.
        :param zero_point_flux_photon: The zero-point flux for the filter in photon units.
        :return: magnitude.
        """
        F0_electron = zero_point_flux_photon * self.collection_area * self.quantum_efficiency * self.exposure_time
        magnitude = -2.5 * np.log10(
            flux / F0_electron.value)
        return magnitude

    def detectable_star_brightness(self, desired_sn=5, filter_name="Ks"):
        """
        Determines the detectable star brightness for a desired SNR and filter name.
        :param desired_sn: Desired signal-to-noise ratio.
        :param filter_name: Name of the filter.
        :return: Detectable star brightness in magnitude.
        """
        # Calculate the total noise
        total_noise = np.sqrt(
            self.sky_noise(filter_name).value
            + self.dark_current_noise().value
            + self.read_noise().value
        )

        # Calculate the flux required for the desired SNR
        required_flux = desired_sn * total_noise * u.electron

        # Get the zero-point for the filter
        zero_point_flux_photon = hmbp.in_zero_vega_mags(filter_name, "HAWKI", "Paranal")

        # Convert this flux to a magnitude using the zero-point magnitude
        mag = self.flux_to_mag(required_flux.value, zero_point_flux_photon)

        return mag

    def plot(self, snr=None, savefig=False, filename="etc.png"):
        """
        Plots the SNR as a function of star brightness.

        :param snr: Signal-to-noise ratio threshold.
        :param savefig: Boolean indicating if the plot should be saved.
        :param filename: Filename to save the plot.
        :return: Plot figure.
        """
        fig, ax = plt.subplots()
        brightness_list = np.arange(26, 8, -0.1)
        snr_aper_list = [self.compute_sn(brightness * u.mag) for brightness in brightness_list]
        ax.scatter(brightness_list, snr_aper_list, s=5, color="black")
        if isinstance(snr, (float, int, list, np.ndarray)):
            ax.axhline(snr, color='red', linestyle='--')
        ax.set_xlim(min(brightness_list), max(brightness_list))
        ax.set_xlabel("Brightness (Magnitude)")
        ax.set_ylabel("SNR")
        ax.grid(True, which="major", ls="--", linewidth=0.5)  # <-- Changed this line
        ax.set_yscale("log")
        ax.set_title("SNR as a function of Star Brightness")
        if savefig:
            if isinstance(filename, str):
                plt.savefig(filename)
            else:
                raise TypeError(
                    f"filename has to be a string, {type(filename)} is given."
                )
        plt.show()
        return fig


if __name__ == "__main__":
    etc = HAWKI_ETC()
    # etc.print_units()
    print("Detectable star brightness:", etc.detectable_star_brightness(5, "Ks"))
