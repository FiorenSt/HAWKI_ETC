from etc_module_hmbp import HAWKI_ETC
from astropy import units as u

def test_sky_noise():
    etc = HAWKI_ETC()
    noise = etc.sky_noise("Ks")
    assert noise.unit == u.electron, f"Expected unit of electron, but got {noise.unit}"

def test_dark_current_noise():
    etc = HAWKI_ETC()
    noise = etc.dark_current_noise()
    assert noise.unit == u.electron, f"Expected unit of electron, but got {noise.unit}"

def test_read_noise():
    etc = HAWKI_ETC()
    noise = etc.read_noise()
    assert noise.unit == u.electron**2, f"Expected unit of electron, but got {noise.unit}"

def test_compute_sn():
    etc = HAWKI_ETC()
    snr = etc.compute_sn(23.3 * u.mag, "Ks")
    assert 4 < snr < 6, f"Expected SNR around 5, but got {snr}"

def test_detectable_star_brightness():
    etc = HAWKI_ETC()
    brightness = etc.detectable_star_brightness(5, "Ks")
    assert 23 < brightness < 24, f"Expected brightness around 23.15, but got {brightness}"

def test_plot_function():
    etc = HAWKI_ETC()
    fig = etc.plot(snr=5)
    assert fig, "Expected a plot figure, but got None."

if __name__ == "__main__":
    test_sky_noise()
    test_dark_current_noise()
    test_read_noise()
    test_compute_sn()
    test_detectable_star_brightness()
    test_plot_function()

    print("All tests passed!")
