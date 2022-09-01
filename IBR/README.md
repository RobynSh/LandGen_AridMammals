# Isolation by resistance 

The R code *IBR.R* documents the process of:
* Calculating mean pairwise Euclidean genetic distance within a 5km area (to match the pixel size of the rasters used in this analysis)
* Running single surface optimisation in ResistanceGA to determine the best performing raster within correlated variable sets
* Running multi-surface optimisation in ResistanceGA across the final raster set, including optimising raster transformations across a maximum of four combined layers, outputting final MLPE model selection, and running bootstrap analysis


### Input files:
* Genlight objects (*SNP_Filtering/SampleClean_outputs/RelClean.gl.[Spp code].rdata*) after SNP filtering.
* Metadata for final cleaned samples (*Data/SampleMetaData_Cleaned.csv*)
* Rasters aggregrated to a 5 km pixel size in the *Rasters_Shapefiles* folder in the home directory. Note that the final set of uncorrelated rasters used in the multisurface optimisation must be placed in a separate folder when running.
* Coordinate (*Data/Coords_[Spp code].csv*) and genetic metric (*Data/EucDist_[marker type/Spp code].csv*) csv files generated in the first section of this R code (read back in when running ResistanceGA)


### Output files:
This code generates the following files (output to the *IBR_outputs* folder):
* Model selection and bootstrap results (results aggregated/averaged over 1000 bootstrap iterations) for all raster combinations included in "all_comb" run (*[Spp code.marker type].ModelSelection.csv* and *[Spp code.marker type].Bootstrap.csv*)
* Folders containing all relevant information for the best-performing model/optimised resistance surfaces for all species and marker types:
 - *Nt.MSAT.TopModel_CF*
 - *Nt.SNP.TopModel_SOMO29.FOR.VRM*
 - *Pc.MSAT.TopModel_CF*
 - *Pc.SNP.TopModel_SOMO29.WII*
 - *Ph.MSAT.TopModel_Clay*
 - *Ph.SNP.TopModel_VRM*
 
 The following outputs can be found within each of the best-performing model subfolders:
      - *[Specific layer name/s here]_commuteDistance_distMat.csv* - Pairwise distance matrix with commute distances (i.e. distance between points calculated using the optimised resistance surface)
      - *[Specific layer name/s here]_DiagnosticPlots.tif* - Diagnostic plots for MLPE model
      - *[Specific layer name/s here]_full.rds* - R data
      - *[Specific layer name/s here].asc* - Optimised resistance surface
      - *[Specific layer name/s here].rds* - R data
      - *MLPE_coeff_Table.csv* - MLPE model summary based on optimised resistance surface
      - *Multisurface_Optim_Summary.txt* - Parameters and details for optimisation
      - *Percent_Contribution.csv* - Percent each layer contributes to the final optimised resistance surface


&nbsp;

&nbsp;

<div align="center">
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
</div>
