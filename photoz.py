import numpy as np
import sys
from astropy.table import Table
import os
import yaml
import json
import dnf



def show_usage():
    """
    Show usage message when there is no arguments.
    """
    usage_message = """
Usage:
1. Edit the configuration file config.json or config.yaml
2. Run the command:
   python photoz.py train.fits valid.fits output.fits
"""
    print(usage_message)

def load_config():
    """
    load config from config.yaml o config.json.

    Returns:
        dict: Dictionary with the configuration
    """
    config_file = None
    config = None

    #look for config.yaml file or config.json
    if os.path.exists("config.yaml"):
        config_file = "config.yaml"
        with open(config_file, "r") as file:
            config = yaml.safe_load(file)
    elif os.path.exists("config.json"):
        config_file = "config.json"
        with open(config_file, "r") as file:
            config = json.load(file)
        
    if config_file:
        print(f"Configuration loaded from {config_file}")
        return config
    else:
        print("Error: No configuration file found (config.yaml or config.json).")
        sys.exit(1)




def process_galaxy_data(file_path, filters, filters_err):
    """Read and process galaxy file."""
    galaxy = Table.read(file_path)
    ngalaxies = len(galaxy)
    print(f'Ngalaxies={ngalaxies}')

    # Data to matrix assignation
    G = np.array([galaxy[f] for f in filters]).T
    Gerr = np.array([galaxy[f_err] for f_err in filters_err]).T
    
    return galaxy, G, Gerr

def main():
    """
    main code
    """
    if len(sys.argv) == 1:
        # There is no arguments, show usage
        show_usage()
        sys.exit(1)
    elif len(sys.argv) == 4:
        # load configuration and process files
        config = load_config()
        # Access specific parameters
        filters = config["filters"]
        filters_err = config["filters_err"]
        zspec= config["zspec"]
        bins = config["bins"]
        start, stop, nbins = bins['start'], bins['stop'], bins['nbins']

        fit = config["fit"]
        metric= config["metric"]
        pdf= config["pdf"]
        nfilters = len(filters)

        print("Filters:", filters)
        print("Bins:", bins)
        print("Fit option:", fit)
        print("pdf option:", pdf)
        print("Metric:", metric)

        # Reading training data
        TRAIN, T, Terr = process_galaxy_data(sys.argv[1], filters, filters_err)
        Ntrain = len(TRAIN)

        # Reading  data data
        VALID, V, Verr = process_galaxy_data(sys.argv[2], filters, filters_err)
        Nvalid = len(VALID)

        # pdz_grid generation
        pdz_grid = np.linspace(start, stop, nbins)

        # Train and Estimator instances
        Train = dnf.dnfTrain(T, Terr, TRAIN[zspec])
        Estimator = dnf.dnfEstimator(Train, pdz_grid, metric=metric, fit=fit, pdf=pdf)

        # Photometric redshift prediction
        z_photo, zerr_e, photozerr_param, photozerr_fit, z1, nneighbors, de1, d1, id1, C, Vpdf = Estimator.photometric_redshift(V, Verr)
        print(f"mean Nneighbors={np.mean(nneighbors)}")

        # Saving results to FITS file
        d = VALID
        d['DNF_Z'] = z_photo
        d['DNF_ZN'] = z1
        d['DNF_ZSIGMA'] = zerr_e
        d['DNF_ZERR_PARAM'] = photozerr_param
        d['DNF_ZERR_FIT'] = photozerr_fit
        d['DNF_D1'] = d1
        d['DNF_NNEIGHBORS'] = nneighbors
        d['DNF_DE1'] = de1
        d['DNF_ID1'] = id1
        if fit:
            d['C'] = C
        if pdf:
            d['pdf'] = Vpdf

        d.write(sys.argv[3], format='fits', overwrite=True)

    else:
        print("Error: Incorrect number of arguments.")
        show_usage()

if __name__ == "__main__":
    main()

    
