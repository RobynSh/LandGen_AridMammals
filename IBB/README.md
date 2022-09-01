# Isolation by barrier analyses

The R code *IBB.R* documents the process of:
* Running Tess3R for each species
* Creating STRUCTURE-style barcharts, pie charts with admixture proportions and cross entropy plots.


### Input files:
* Genlight objects (*SNP_Filtering/SampleClean_outputs/RelClean.gl.[Spp code].rdata*) after SNP filtering
* Shape files of the study region (*PilbaraIBRA.shp*) and of Western Australia (*IBRA_regions_WA.shp*) in the *Rasters_Shapefiles* folder in the home directory.
* Structure format files for microsatellite datasets (*[Spp code]_LevySubset_Microsats.str*) for each species (note that the excel microsatellite data are also provided in this folder, from which structure files were generated)


### Output files:
This code generates the following files (output to the *IBB_outputs* folder):
* *CrossVal.Tess.Multi.pdf* - Figure: Tess3R Cross entropy plot (for choosing best value for K) across all species and marker types
* *Tess.obj.FullDat.best100reps.K1_7.[Spp code]_[marker type].R* - Tess3R full analysis saved as an R object for each species/marker type
* *TessPieBar_K2.Nt.pdf* - Figure: Tess3r Structure-like barchart and ancestry proportions displayed as pie charts on a map of the study region   


&nbsp;

&nbsp;

<div align="center">
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
</div>
