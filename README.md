# DNF: Directional Neighbourhood Fitting

DNF is a nearest-neighbor approach for photometric redshift estimation. DNF is an strategy that improves the kNN approaches. DNF computes the photo-z hyperplane that best fits the directional neighbourhood of a photometric galaxy in the training sample. A detailed description of DNF is available [here](https://arxiv.org/abs/1511.07623).

If you have any questions or suggestions, please don't hesitate to contact us at laura.toribio@ciemat.es.

# DNF User Manual
##  How to Get Started with DNF

**1. Create a New Directory:** 
You can name your directory as you prefer. For example, you can use the following command to create a directory called `photozDNF`:

```bash
mkdir photozDNF
```

**2. Download DNF and the Data Folder:**

- **Download the DNF, Photoz and config codes:**  
  Visit the official [DNF GitHub space](https://github.com/ltoribiosc/DNF_photoz) to download the latest version of:
  - `dnf.py`: This is the main code that calculates the photo-z. You don’t need to modify anything in this file.  
  - `photoz.py`: This code handles loading the configuration and reading the training data (`train.fits`) as well as the validation data or the data for which you want to calculate the photo-z (e.g., `valid.fits`). You don’t need to modify anything in this file. 
  - `config.json`: This file contains the selected configuration. You will need to modify this file, but don’t worry, we’ll explain how to do it later.  


- **Download the Data Folder:**  
  Download the data folder required to run the example. Make sure this folder contains all the necessary files. In the downloaded data folder, you will find two important files:
  - `train.fits`: This file is used for training DNF.
  - `valid.fits`: Use this file to test the quality of photo-z results.

You are now ready to begin using DNF with the provided data!

- **Run DNF with Sample Data:**
To run our example, you need to execute the following command:

```bash
python3 photoz.py data/train.fits data/valid.fits
```
PREGUNTAR A JUAN ESTO PORQUE HA CAMBIADO

In this command:
  - `data/train.fits` is the training data file, which contains the photometric information and spectroscopic redshift data for a sample of galaxies.
  - `data/valid.fits` is the validation data file, which contains the galaxies for which you want to calculate the photo-zs.

You will need to modify these paths when using your own data.

You have just run DNF. Congratulations on a successful run!

If everything has gone well, you should now have a new file named `results.fits`. This file contains all the galaxy information from the valid file, enriched with the new data provided by DNF.

If you encounter any errors or issues, please don't hesitate to contact us at laura.toribio@ciemat.es.

## Running DNF with Your Own Data
If you're interested in running DNF with your own data (that's probably why you're here.), follow these steps:

**1. Prepare Your Data: DNF requires two input files:**
  - `train.fits`: This file contains photometric information and spectroscopic redshift data for a sample of galaxies.
  - `valid.fits`: You will also need the file for which you want to calculate photoz. In this case, it should be a sample of galaxies with photometric information, exactly matching the filters used in the train file.

In the example, these files were named `train.fits` and `valid.fits`, but you can use whatever names you prefer. We recommend that you save your input files in the `data` folder, but this is not mandatory. You can create your own directory, just remember to include the correct path when running DNF.

**2. Modify the config.json file if is necessary:**
The `config.json` file contains all the necessary information for DNF to run with the required customization for each type of data and analysis. In the example, `config.json` was configured for the data in `train.py`, but now we will explain everything it contains so you can learn how to modify it.

  - **"filters"**: Contains the names of the filters used by DNF to calculate the photo-zs, which are available in the `train.py` file. In the example, we have four filters named:
    - `"BDF_MAG_G_CORRECTED"`
    - `"BDF_MAG_R_CORRECTED"`
    - `"BDF_MAG_I_CORRECTED"`
    - `"BDF_MAG_Z_CORRECTED"`
  You should modify these names to match those in your data.

  - **"filters_err"**: These are the error values corresponding to your data. You should modify these names to match those in your data.

  - **"bins"**: Contains the information for calculating the DNF, including the start, end, and number of bins. Unless you need to change this for a specific reason, you can leave the default configuration.

  - **"fit"**: Choose `True` or `False` depending on whether you want *** (please clarify what *** refers to).

  - **"pdf"**: Select `True` or `False` depending on whether you want to calculate the PDFs.

  - **"metric"**: DNF has three types of metrics for determining the DNF. Choose the one that best fits your needs: `ENF`, `ANF`, or `DNF`. DNF provides three different metrics for calculating photoz:
    - ENF: Euclidean neighbourhood. it's a common distance metric used in kNN (k-Nearest Neighbors) for photometric redshift prediction.
    - ANF: uses normalized inner product for more accurate photo-z predictions.
    - DNF: combines Euclidean and angular metrics, improving accuracy, especially for larger neighborhoods, and maintaining proportionality in observable content.
  You can learn more about these three metrics by referring to the following resource: [here](https://arxiv.org/abs/1511.07623).

  - **"output_file"**: Specify the path and name of the output file where the photo-z results will be saved.

**3. Launching DNF with Your Data**

Now, you can run DNF with your data by entering the following command in your console:

```bash
python3 photoz.py your_train_data.fits your_data.fits
```
