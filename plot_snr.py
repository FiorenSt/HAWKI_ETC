from etc_module_hmbp import HAWKI_ETC

def main():
    """
    Main function to generate and display the SNR plot.
    """

    # Create an instance of the HAWKI_ETC class and plot
    etc = HAWKI_ETC()
    etc.plot(snr=5)

if __name__ == "__main__":
    main()
