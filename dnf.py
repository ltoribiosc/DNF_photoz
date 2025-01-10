"""
This software is distributed without any warranty
Please, include the reference to the paper below in your publications if you use this software:
De Vicente, J.; Sanchez, E.; Sevilla-Noarbe, I.,"DNF - Galaxy photometric redshift by Directional Neighbourhood Fitting", MNRAS, Vol. 459, Issue 3, pag 3078-3088, 2016.
"""


#!/usr/bin/python
__author__ = "Juan de Vicente"
__copyright__ = "Copyright 2015, Juan de Vicente"
__version__ = "5.0"
__email__= "juan.vicente@ciemat.es"

"""
DNF 5.0 (2024)
Authors: Juan de Vicente, Laura Toribio
Update: class implementation
contains three alternative functions to compute the photometric redshift of a sample of galaxies:
enf: Euclidean Neighborhood Fit 
dnf: Directional Neighborhood Fit
anf: Angular Neighborhood Fit
"""


import math
import numpy as np
from sklearn import neighbors
import os
from astropy.io import fits
from sklearn.neighbors import KNeighborsRegressor


class dnfTrain:
    """
    Computes euclidean neighbors as a base for other metrics.

    Parameters:
    ----------
    T : ndarray, shape (Ntrain, nfilters)
        Magnitudes or fluxes of the training sample.
    Terr : ndarray, shape (Ntrain, nfilters)
        Photometric errors corresponding to the training sample.

    Attributes:
    ----------

    clf: neighbors regression
    Tnorm : ndarray
        Precomputed norms of the training set magnitudes.
    
    Methods:
    ----------
    def __init__(self, T, Terr, z)

    """
     
    def __init__(self, T, Terr, z):  
        
        '''
        Load training data and initialize variables
        '''

        self.T=T
        self.Terr=Terr
        self.z=z
        
        #Training euclidean metric
        self.clf=neighbors.KNeighborsRegressor()
        self.clf.fit(self.T, self.z)
        
        # Training variables
        self.Tnorm=np.linalg.norm(self.T,axis=1)
        

class dnfEstimator:
    """
    Computes photometric redshifts (photo-z) using the Directional Neighborhood Fit (DNF) method.

    DNF predicts the redshift of galaxies based on training data and photometric samples, using 
    various distance metrics and fitting strategies.

    Parameters:
    ----------
    
    V : ndarray, shape (Nvalid, nfilters)
        Magnitudes or fluxes of the photometric sample (validation/science set).
    Verr : ndarray, shape (Nvalid, nfilters)
        Photometric errors corresponding to the photometric sample.
    zgrid : ndarray
        1-dimensional array representing redshift bins for photo-z PDF estimation.
    metric : str, optional
        Distance metric to use: 'ENF' (Euclidean), 'ANF' (Angular), or 'DNF' (Directional). Default is 'ANF'.
    fit : bool, optional
        If True, computes the coefficients of the hyperplane fit for each galaxy. Default is False.
    pdf : bool, optional
        If True, computes the redshift Probability Density Function (PDF) for each galaxy. Default is False.
    presel : int
        Preselection threshold for neighbors based on Euclidean distance. Default is 80.
    Nneighbors : int
        Number of neighbors to use for the hyperplane predictor. Default is 500.

    Attributes:
    ----------
   
    fitIterations : int
        Number of iterations to refine the fitting by removing outliers.
    badtag : int or float
        Value used to tag galaxies without reliable photo-z estimates.
    nfilters : int
        Number of filters or features in the dataset.
    Ntrain : int
        Number of galaxies in the training set.
    Ts, Tsnorm : ndarray
        Preselected training data and their norms (initialized later).
    photoz, photozmean, photozfit : ndarray
        Estimated photo-z values for the photometric sample.
    photozerr, photozerr_param, photozerr_fit : ndarray
        Uncertainty estimates for the photo-z predictions.
    nneighbors : ndarray
        Number of neighbors used for the photo-z estimation of each galaxy.
    NEIGHBORS : structured ndarray
        Stores neighbor indices, distances, and redshifts for each galaxy.
    Vpdf : ndarray
        Redshift Probability Density Functions (PDFs) for the photometric sample.
    zpdf, wpdf : ndarray
        Redshift and weight arrays of neighbors used for PDF estimation.

    Methods:
    ----------
 ******************
    def __init__(self, dnfTrain, zgrid, metric='ANF', fit=False, pdf=False,Nneighbors=80,presel=500):
    def validate_columns(self,V):
    def photometric_redshift(self, V, Verr):
        return (photoz, photozerr, photozerr_param, photozerr_fit, z1, nneighbors, de1, d1, id1, C, Vpdf,)
    def preselection(self, V, Verr):
        return (NEIGHBORS, Ts, Tsnorm, de1, )
    def metric_computation(self,V, NEIGHBORS)
        return NEIGHBORS
    def compute_euclidean_distance(self, V, Ts)
        return d
    def compute_angular_distance(self, V, Ts, Tsnorm)
        return alpha
    def compute_directional_distance(self, V, Ts, Tsnorm)
        return d
    def compute_photoz_mean_routliers(self, NEIGHBORS, Verr):
        return (photoz, photozerr, photozerr_param, photozerr_fit, z1, nneighborgs, Vpdf,  NEIGHBORS, )
    def compute_photoz_fit(self, V, Verr, NEIGHBORS): 
        return (photoz, photozerr, photozerr_param, photozerr_fit, nneighborgs, C.)                  
    def compute_pdfs(self,zpdf, wpdf):
        return Vpdf    
    def calculate_pdf_vectorized(self, zpdf, wpdf):
        return Vpdf

 ***********

    """
    def __init__(self, dnfTrain, zgrid, metric='ANF', fit=False, pdf=False,Nneighbors=80,presel=500):
        
        '''
        Load target data and initialize variables. In addition load dnfTrain instance.
        '''

        # Configuration parameters
        self.Nneighbors=Nneighbors
        self.presel=presel
        self.fitIterations = 4
        self.badtag = 99.0
        self.fit = fit
        self.zgrid = zgrid
        self.pdf = pdf
        self.metric = metric
        self.presel= presel
        

        # Training variables
        self.T=dnfTrain.T
        self.Terr=dnfTrain.Terr
        self.z=dnfTrain.z
        self.clf=dnfTrain.clf
        self.Tnorm=dnfTrain.Tnorm
        self.Ntrain, self.nfilters = self.T.shape
        self.Nneighbors_presel = self.presel if self.Ntrain > self.presel else self.Ntrain
        self.Te = np.hstack([self.T, np.ones((self.Ntrain, 1), dtype='double')])


    def validate_columns(self,V):
        """
        Validates that the columns of T and V have the same names.
        
        Parameters
        ----------
        T : np.ndarray
            Training data.
        V : np.ndarray
            Validation data.
        
        Raises
        ------
        ValueError
            If the column names of T and V do not match.
        """
        if not np.array_equal(self.T.dtype.names, V.dtype.names):
            raise ValueError("The columns of T and V do not match. Please ensure that both T and V have the same features.")


    def photometric_redshift(self, V, Verr):
        """
        Compute the photometric redshifts for the validation or science sample.
    
        Returns
        -------
            - photoz: Estimated photometric redshift.
            - photozerr: Error on the photometric redshift.
            - photozerr_param: Redshift error due to parameters.
            - photozerr_fit:Redshift error due to fit.
            - z1: Closest redshift estimate.
            - nneighbors: Number of neighbors considered.
            - de1: Distances Euclidea to the closest neighbor.
            - d1: Distances to the closest neighbor.
            - id1: Index of the closest neighbor.
            - C: Additional computed parameters.
            - zpdf: Matrix containing the redshifts of neighboring galaxies.
            - wpdf: Matrix of weights corresponding to the neighboring redshifts.
            - Vpdf: Probability Density Functions (PDFs) for the photometric redshifts of the validation set.
        """
        
        C=0 #coefficients by default

        # Step 1: Preselection
        NEIGHBORS, Ts, Tsnorm, de1=self.preselection(V, Verr)

        # Step 2: Metric computation
        NEIGHBORS=self.metric_computation(V, NEIGHBORS, Ts, Tsnorm) 
        

        # Step 3: Compute mean redshift removing outliers
        photoz, photozerr, photozerr_param, photozerr_fit, z1, d1, id1, nneighbors, Vpdf, NEIGHBORS=self.compute_photoz_mean_routliers(NEIGHBORS, Verr)

        # Step 4: Optional fitting of redshifts
        if self.fit:
            photoz, photozerr, photozerr_param, photozerr_fit, nneighbors, C=self.compute_photoz_fit( V, Verr, NEIGHBORS, photoz, photozerr, photozerr_param, photozerr_fit) 
        
        # Return
        return (
            photoz,
            photozerr,
            photozerr_param,
            photozerr_fit,
            z1,
            nneighbors,
            de1,
            d1,
            id1,
            C,
            Vpdf,
         )


    def preselection(self, V, Verr):
        """
        Perform the preselection process for photometric redshift estimation.
        """
        # Ensure V is set
        if V is None:
            raise ValueError("Validation data 'V' is not set. Ensure it is initialized.")

        # Validation variables
        self.validate_columns(V)
        
        
        Nvalid=V.shape[0]
        

        # Compute distances and indices for preselected neighbors
        Ndistances, Nindices = self.clf.kneighbors(V, n_neighbors=self.Nneighbors_presel)

        # Handle cases where distance is zero
        mask = Ndistances[:, 0] == 0
        Ndistances[mask] = np.roll(Ndistances[mask], -1, axis=1)
        Nindices[mask] = np.roll(Nindices[mask], -1, axis=1)
        Ndistances[mask, -1] = Ndistances[mask, 0]
        Nindices[mask, -1] = Nindices[mask, 0]

        # Store distances and indices
        de1 = Ndistances[:, 0]
        

        # Initialize NEIGHBORS array to store indices, distances, and redshifts
        NEIGHBORS = np.zeros(
            (Nvalid, self.Nneighbors_presel),
            dtype=[('index', 'i4'), ('distance', 'f8'), ('z', 'f8')]
        )
        NEIGHBORS['index'] = Nindices
        NEIGHBORS['distance'] = Ndistances
        NEIGHBORS['z'] = self.z[Nindices]

        # Extract training preselection data
        Ts = self.T[Nindices]
        Tsnorm = self.Tnorm[Nindices]

        return NEIGHBORS, Ts, Tsnorm, de1

    def metric_computation(self,V, NEIGHBORS, Ts, Tsnorm):
        """
        Compute distances based on the selected metric, sort neighbors, and store the closest neighbors.
        """              
        # Compute distances based on the chosen metric
        if self.metric == 'ENF':
            NEIGHBORS['distance'] = self.compute_euclidean_distance(V, Ts)
        elif self.metric == 'ANF':
            NEIGHBORS['distance'] = self.compute_angular_distance(V, Ts, Tsnorm)
        elif self.metric == 'DNF': 
            NEIGHBORS['distance'] = self.compute_directional_distance(V, Ts, Tsnorm)

        # Sort the top Nneighbors
        NEIGHBORS = np.sort(NEIGHBORS, order='distance', axis=1)
       
        # Select redshift of closer neighbors until n numbers of neighbors
        NEIGHBORS=NEIGHBORS[:,:self.Nneighbors]
        
        return NEIGHBORS 


    def compute_euclidean_distance(self, V, Ts):
        """
        Compute distances based on Euclidean metric.
        """      
         
        D = V[:, np.newaxis,:] - Ts
        Dsquare = D[:] * D[:]
        D2 = np.sum(Dsquare[:], axis=2)
        d = np.sqrt(D2)
        return d


    def compute_angular_distance(self, V, Ts, Tsnorm):
        """
        Compute distances based on angular (ANF) metric.
        """    
         
        Vnorm = np.linalg.norm(V, axis=1) #[:,np.newaxis])
        pescalar = np.sum(V[:, np.newaxis,:] * Ts, axis=2) 
        normalization = Vnorm[:, np.newaxis] * Tsnorm 
        NIP = pescalar / normalization
        alpha = np.sqrt(1.0 - NIP**2)
        return alpha
         
         
    def compute_directional_distance(self, V, Ts, Tsnorm):
        """
        Compute distances based on directional (DNF) metric.
        """ 
         
        d1 = self.compute_euclidean_distance(V, Ts) 
        d2 = self.compute_angular_distance(V, Ts, Tsnorm)
        d = d1 * d2
        return d

    
    def compute_photoz_mean_routliers(self, NEIGHBORS, Verr):
        """
        Compute the mean photometric redshift removing outliers
        """
        
        # Extract distances and redshifts from neighbors 
        distances = NEIGHBORS['distance']
        zmatrix = NEIGHBORS['z']
        indices = NEIGHBORS['index']

        # Store nearest redshift, distance and index
        z1= zmatrix[:,0]
        d1 = distances[:, 0]
        id1 = indices[:,0]
        
        # --- Outlier Detection and Weighting ---
        # Calculate mean distance for each sample
        mean_absolute_deviation =distances.mean(axis=1)
        # Define the threshold for outlier detection
        threshold = mean_absolute_deviation #Adjust multiplier if needed (e.g., *2)
       
        # Create a mask for non-outlier distances
        outliers_weights = distances < threshold[:, None]
       
        # Update the number of valid neighbors per sample
        nneighbors = np.sum(outliers_weights, axis=1)
        cutNneighbors = np.max(nneighbors) # Maximum number of valid neighbors
            
        # --- Distance Weighting ---
        # Compute inverse distances for weighting
        inverse_distances = 1.0 / distances
         
        # Apply the outlier weigh to inverse distances and distances
        inverse_distances = inverse_distances * outliers_weights
        distances = distances * outliers_weights

        # Normalize weights by the sum of inverse distances per sample
        row_sum = inverse_distances.sum(axis=1, keepdims=True)
        wmatrix = inverse_distances / row_sum
        wmatrix = np.nan_to_num(wmatrix) # Handle potential NaN values from division by zero
            
        # --- Photometric Redshift Computation ---
        # Compute the weighted mean redshift for each sample
        photozmean = np.sum(zmatrix * wmatrix, axis=1)
        photoz = photozmean
         
        # --- Redshift Error Computation ---
        # Compute the standard deviation of redshifts (fit error)         
        zt = photozmean[:, np.newaxis] #column vector
        zerrmatrix = (zmatrix - zt) ** 2

        # Compute the error based in the fit and the parameters
        Verrnorm = np.linalg.norm(Verr, axis=1)
        photozerr_fit = np.sqrt(np.sum(zerrmatrix * wmatrix, axis=1))
        photozerr_param = np.abs(0.2 * Verrnorm * (1.0 + photoz))
         
        # Combine errors to calculate the total redshift error
        photozerr = np.sqrt(photozerr_param**2 + photozerr_fit**2)

        # --- PDF Computation Setup ---
        # Select the top Nneighbors redshifts and weights for PDF computation
        zpdf = zmatrix[:, :cutNneighbors]
        wpdf = wmatrix[:, :cutNneighbors]

        # Update NEIGHBORS array to include only the top Nneighbors
        NEIGHBORS = NEIGHBORS[:, :cutNneighbors]

        if self.pdf:
            Vpdf=self.compute_pdfs(zpdf, wpdf)

        return photoz, photozerr, photozerr_param, photozerr_fit, z1, d1, id1, nneighbors, Vpdf,  NEIGHBORS

 
    def compute_photoz_fit(self, V, Verr, NEIGHBORS, photoz, photozerr, photozerr_param, photozerr_fit): 
        """
        Compute the photometric redshift fit by iteratively removing outliers.
        """   
        # Initialize output parameters
        Nvalid=V.shape[0]
        # Initialize output variables with default values        
        '''
        photoz=self.badtag * np.ones(Nvalid, dtype='double')
        photozerr=self.badtag * np.ones(Nvalid, dtype='double')
        photozerr_param=self.badtag * np.ones(Nvalid, dtype='double')
        photozerr_fit=self.badtag * np.ones(Nvalid, dtype='double')
        '''
        photoz=photoz
        photozerr=photozerr
        photozerr_param=photozerr_param
        photozerr_fit=photozerr_fit
        print('NEIGHBORS',NEIGHBORS)
        
        #photozfit=self.badtag *np.zeros(Nvalid, dtype='double')
        photozfit=99.0 *np.zeros(Nvalid, dtype='double')
        nneighbors=np.zeros(Nvalid,dtype='double')
        
        if self.fit:
            C = np.zeros((Nvalid, self.nfilters + 1), dtype='double')
        else:
            C = 0
        # Increase dimensionality of validation and training data for offsets in fit
        Ve = np.hstack([V, np.ones((Nvalid, 1), dtype='double')])
        

        # Loop over all validation samples
        for i in range(0, Nvalid):
           NEIGHBORSs = NEIGHBORS[i] # Get neighbors for the current sample
           nneighbors[i]=len(NEIGHBORSs)
           # Perform iterative fitting
           for h in range(0, self.fitIterations):
               # Build the design matrix (A) and target vector (B) for the neighbors
               A = self.Te[NEIGHBORSs['index']]  
               B = self.z[NEIGHBORSs['index']]
                
               # Solve the least squares problem
               X = np.linalg.lstsq(A, B)
               residuals = B - np.dot(A, X[0]) # Compute residuals
                
               # Identify outliers using a 3-sigma threshold
               abs_residuals = np.abs(residuals)      
               sigma3 = 3.0 * np.mean(abs_residuals)
               selection = (abs_residuals < sigma3)

               # Update the number of selected neighbors
               nsel = np.sum(selection)
               nneighbors[i] = nsel
               
               # If enough neighbors remain, update NEIGHBORSs; otherwise, stop iteration
               if nsel > 10:
                   NEIGHBORSs = NEIGHBORSs[selection]
               else:
                   break
                    
               # Save the solution vector
           C[i] = X[0]
                
           # Compute the photometric redshift fit for the current sample
           photozfit[i] = np.inner(X[0], Ve[i])
           nneighbors[i]= NEIGHBORSs.shape[0] 
           # Calculate error metrics
           if X[1].size != 0:      # Check if residuals from the least squares solution exist
              # Compute the error based in parameters
              param_error = np.abs(X[0][:-1]) * Verr[i]
              photozerr_param[i] = np.sqrt(np.sum(param_error**2)) #,axis=1))
                  
              # Compute the error based in fit
              photozerr_fit[i] = np.sqrt(X[1] / nneighbors[i])
                  
              # Combine errors to get total redshfit error
              photozerr[i] = np.sqrt(photozerr_param[i]**2 + photozerr_fit[i]**2)
                  
        # Assign the computed redshift fits to the output variable
        photoz = photozfit
                  
        return photoz, photozerr, photozerr_param, photozerr_fit, nneighbors, C       
                  
    def compute_pdfs(self, zpdf, wpdf):
       """
       Compute the pdfs from neighbor redshifts and weights. 
       """
       
       Nvalid=zpdf.shape[0]

       # Initialize PDF-related variables if required
       if self.pdf:
            Vpdf = np.zeros((Nvalid, len(self.zgrid)), dtype='double')
       else:
            Vpdf = 0

       bin_indices = np.digitize(zpdf,self.zgrid)  # Indices de los bins para cada elemento de zpdf (l x m)
       
       # Mask matrix for each bin (l x m x len(bins)) 
       masks = (bin_indices[..., np.newaxis] == np.arange(len(self.zgrid))) 
           
       # mask per weights and sum by columns
       histograms = np.einsum('ij,ijk->ik', wpdf, masks)
       Vpdf=histograms

       return Vpdf
    

       
    def calculate_pdf_vectorized(self, zpdf, wpdf):
        """
        Compute the probability density function from zpdf and wpdf 
        """
        # Validacion de entrada
        if zpdf.shape != wpdf.shape:
            raise ValueError("zpdf y wpdf should have the same dimensions.")
        
        # Normalization constant
        normalizer = 1 / np.sqrt(2 * np.pi)
        
        # Expand dimensions for computing differences with broadcasting
        # x_grid has the shape (1, 1, nbins) for operating with A in vectorized way
        diffs = self.zgrid[None, None, :] - zpdf[:, :, None]

        # Compute gaussians over zgrid: exp(-0.5 * diffs^2) with normalization
        gaussians = normalizer * np.exp(-0.5 * diffs**2)
        
        # Adding weights self.wpdf expanding to operate with gaussians
        # B[:, :, None] has shape (m, n, 1) for multiplying with gaussians
        weighted_gaussians = gaussians * wpdf[:, :, None]
        
        # Sum the contributions of each point in each file
        Vpdf = np.sum(weighted_gaussians, axis=1)

        return Vpdf




        

