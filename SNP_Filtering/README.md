# SNP Filtering

The code *DArTSNPFilt.R* documents the filtering pipeline for each species (Nt = *Ningaui timealeyi*, Pc = *Pseudomys chapmani*, Ph = *P. hermannsburgensis*). 

### Input files:  
* Supplied by DArT: 1-row SNP csv - *Data/[Spp code]_DArT_1Row.csv*
* Supplied by DArT: 2-row SNP csv - *Data/[Spp code]_DArT_2Row.csv* 
* Supplied by DArT: Read count csv - *Data/[Spp code]_DArT_Counts.csv* 
* Sample meta-data csv (including ID, pop, lon/lat, other info) - *Data/SampleMetaData.csv*
* Mus musculus chromosome info to match BLAST data in rodents - *Data/Mus.musculus_ChromInfo*  
* Shape file of study region in *Rasters_Shapefiles* folder in the home directory - *PilbaraIBRA.shp* 

### Output files:  
Individual plots and the filtered data are saved in the *Filtering_outputs* folder:
- *[Spp code]_Avg.Ind.ReadCount.jpg* - Combined smear plot, individual call rate, and map displaying samples coloured by read count
- *[Spp code]_RawHistograms.Mapped-UnMapped.jpg* - Histograms of important metrics across SNPs mapped to a genome vs. unmapped SNPs
- *[Spp code]_Chrom.jpg* - SNP positions in the genome
- *[Spp code]_Corr.FinalFilt.jpg* - Correlations between important metrics in the final filtered data set
- *[Spp code]_Corr.QualFilt.jpg* - Correlations between important metrics in the data set after quality control filters
- *[Spp code]_Corr.Raw.jpg* - Correlations between important metrics in the raw data set
- *[Spp code]_Final.PCoA.jpg* - PCoA of the final filtered dataset
- *[Spp code]_Ind.ReadCount.jpg* - Combined plot showing how individual read counts effect the called genotype, allele balance and important population level metrics
- *[Spp code]_Priv.Al.jpg* - Private alleles in each population before and after MAF filter
- *[Spp code]_SummStat.FinalFilt.jpg* - Histograms and scatter plots for important summary statistics on the final filtered data set
- *[Spp code]_SummStat.QualFilt.jpg* - Histograms and scatter plots for important summary statistics after quality control filtering
- *[Spp code]_SummStat.Raw.jpg* - Histograms and scatter plots for important summary statistics before filtering
- *gl.[Spp code]_FinalFilt.rdata* - Genlight for the final filtered SNP data set, saved as an R object
- *[Spp code]_Mapped_gl.CR.RC.Rp.MAF.lD.gds* - gds file for use in SNPRelate (LD pruning - mapped loci)
- *[Spp code]_UnMapped_gl.CR.RC.Rp.MAF.lD.gds* - gds file for use in SNPRelate (LD pruning - unmapped loci)

See R code for more details.


# Sample Cleaning

Following SNP filtering, we carried out sample cleaning using the R code *SampleCleaning.R*, by removing highly related individuals. We calculate Wang's (2002) pairwise relatedness estimate and remove one individual from the pair if relatedness is greater >=0.24 to avoid biasing population genetic analyses. Finally, we also generate sample buffers of 5km to group samples in similar locations for later use.

### Input files:
* Sample meta-data file - *Data/SampleMetaData.csv*
* The final filtered genlight generateds (SNP filtering) - *Filtering_outputs/gl.[Spp code]_FinalFilt.rdata*

### Output files:
* *RelClean.gl.[Spp code].rdata* - Final cleaned genlight with highly related individuals removed, saved as an R object
* *[Spp code]_related.txt* - Data preparation for Wang's 2002 pairwise relatedness for each species
* *[Spp code]_wang0.24.csv* - Results for Wang's 2002 pairwise relatedness for each species. Note that only pairs with relatedness greater than 0.24 are output.
* *SamplesToDelete_Related.csv* - List of samples to remove that are highly related.
See R code for more details.


&nbsp;

&nbsp;
<div align="center">
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
</div>
