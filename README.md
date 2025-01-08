# DNF: Directional Neighbourhood Fitting

DNF is a nearest-neighbor approach for photometric redshift estimation developed at the CIEMAT (Centro de Investigaciones Energéticas Medioambientales y Tecnológicas). DNF is an strategy that improves the kNN approaches. DNF computes the photo-z hyperplane that best fits the directional neighbourhood of a photometric galaxy in the training sample. A detailed description of DNF is available [here](https://arxiv.org/abs/1511.07623).

If you have any questions or suggestions, please don't hesitate to contact us at laura.toribio@ciemat.es and/or juan.vicente@ciemat.es.

### Citing This Work

If you use this code in your research, we kindly ask you to cite the following article:  

De Vicente, J., Sánchez, E., Sevilla-Noarbe, I., 2016, Directional Neighbourhood Fitting (DNF): photometric redshifts using nearest neighbours, MNRAS, 459, 3078.  
[DOI: 10.1093/mnras/stw824](https://doi.org/10.1093/mnras/stw824)  
[arXiv:1511.07623](https://arxiv.org/abs/1511.07623)  

Here is the BibTeX entry:  

```bibtex
@article{DeVicente2016,
    author = {De Vicente, J. and Sánchez, E. and Sevilla-Noarbe, I.},
    title = {Directional Neighbourhood Fitting (DNF): photometric redshifts using nearest neighbours},
    journal = {Monthly Notices of the Royal Astronomical Society},
    volume = {459},
    issue = {3},
    pages = {3078-3091},
    year = {2016},
    doi = {10.1093/mnras/stw824},
    eprint = {1511.07623},
    archivePrefix = {arXiv}
}
```


# DNF User Manual
##  How to Get Started with DNF

**1. Create a New Directory:** 
You can name your directory as you prefer. For example, you can use the following command to create a directory called `photozDNF`:

```bash
mkdir DNF_photoz
```

**2. Download DNF and the Data Folder:**
- **Clone directly the repository:**

If you are familiar with using GitHub, you can clone the repository where you will find everything you need. Type git `clone` and paste the URL:

```
git clone https://github.com/ltoribiosc/DNF_photoz.git
```
If you are not famliar with GitHub:
- **Download the DNF, Photoz and config codes:**  
  Visit the official [DNF GitHub space](https://github.com/ltoribiosc/DNF_photoz) to download the latest version:
  - `dnf.py`: This is the main code that calculates the photo-z. You don’t need to modify anything in this file.  
  - `photoz.py`: This code handles loading the configuration and reading the training data (`train.fits`) as well as the validation data or the data for which you want to calculate the photo-z (e.g., `valid.fits`). You don’t need to modify anything in this file. 
  - `config.json`: This file contains the selected configuration. You will need to modify this file, but don’t worry, we’ll explain how to do it later.  


- **Download the Data Folder:**  
  Download the data folder required to run the example. Make sure this folder contains all the necessary files. In the downloaded data folder, you will find two important files:
  - `train.fits`: This file is used for training DNF.
  - `valid.fits`: Use this file to test the quality of photo-z results.



You are now ready to begin using DNF with the provided data!

**3. Run DNF with Sample Data:**
To run our example, you need to execute the following command:

```bash
python3 photoz.py data/train.fits data/valid.fits output_results.fits
```

You have just run DNF. 
Congratulations on a successful run!

If everything has gone well, you should now have a new file named `output_results.fits`. This file contains all the galaxy information from the `valid.fits` file,  with the photozs and other values provided by DNF.

If you encounter any errors or issues, please don't hesitate to contact us at laura.toribio@ciemat.es and/or juan.vicente@ciemat.es.

In this command:
  - `data/train.fits` is the training data file, which contains the photometric information and spectroscopic redshift data for a sample of galaxies.
  - `data/valid.fits` is the validation data file, which contains the galaxies for which you want to calculate the photo-zs.
  - `output_results.fits` is the results data file, which contins the galaxies of `data/valid.fits` files with the photo-zs and other parameters.

You will need to modify these paths when using your own data. If you want to save the output files in a specific folder or to have the inputs data in other folder, you must also indicate this in the path. In this example, `train.fits` and `valid.fits` are in the `data folder` while `output_results.fits` will be save in the `DNF_photoz folder`. 

## Running DNF with Your Own Data
If you're interested in running DNF with your own data (that's probably why you're here.), follow these steps:

**1. Prepare Your Data: DNF requires two input files:**
  - `train.fits`: This file contains photometric information and spectroscopic redshift data for a sample of galaxies.
  - `valid.fits`: You will also need the file for which you want to calculate photoz. In this case, it should be a sample of galaxies with photometric information, exactly matching the filters used in the train file.

In the example, these files were named `train.fits` and `valid.fits`, but you can use whatever names you prefer. We recommend that you save your input files in the `data` folder, but this is not mandatory. You can create your own directory, just remember to include the correct path when running DNF.

**2. Modify the config.json file if is necessary:**
The `config.json` file contains all the necessary information for DNF to run with the required customization for each type of data and analysis. In the example, `config.json` was configured for the data in `train.fits`, but now we will explain everything it contains so you can learn how to modify it.

  - **"filters"**: Contains the names of the filters used by DNF to calculate the photo-zs, which are available in the `train.fits` file. In the example, we have four filters named:
    - `"MAG_G"`
    - `"MAG_R"`
    - `"MAG_I"`
    - `"MAG_Z"`
      
  You should modify these names to match those in your data. You can include as many filters as you want in your analysis. You only have to take care to put the correct name that appears in the `train.fits` file.

  - **"filters_err"**: These are the error values corresponding to your data. You should modify these names to match those in your data.

  - **"bins"**: Contains the information for calculating the PDFs, including the start, end, and number of bins. Unless you need to change this for a specific reason, you can leave the default configuration.

  - **"fit"**: Choose `True` or `False` depending on whether you want use the fit fuction. The `fit` function performs a better adjustment for the photo-z calculation, but it increases the execution time.

  - **"pdf"**: Select `True` or `False` depending on whether you want to calculate the PDFs.

  - **"metric"**: DNF has three types of metrics for determining the DNF. Choose the one that best fits your needs: `ENF`, `ANF`, or `DNF`. DNF provides three different metrics for calculating photoz:
    - ENF: Euclidean neighbourhood. it's a common distance metric used in kNN (k-Nearest Neighbors) for photometric redshift prediction.
    - ANF: uses normalized inner product for more accurate photo-z predictions.
    - DNF: combines Euclidean and angular metrics, improving accuracy, especially for larger neighborhoods, and maintaining proportionality in observable content.
      
For more than 4 filters, **we recommend using ANF**. You can learn more about these three metrics by referring to the following resource: [here](https://arxiv.org/abs/1511.07623).

  - **"z_spec"**: Contains the name of the spectroscopic redshift column in the `train.fits` data.

**3. Launching DNF with Your Data**

Now, you can run DNF with your data by entering the following command in your console:

```bash
python3 photoz.py your_train_data.fits your_data.fits your_output_data_name.fits
```
