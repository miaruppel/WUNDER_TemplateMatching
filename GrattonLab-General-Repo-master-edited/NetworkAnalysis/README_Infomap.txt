Updated 4/2022 - Running vertex-wise Infomap on Quest

Started from WashU versions of Infomap scripts

Vertex-wise Infomap can be run using the batch script 'Run_Infomap.sh' as follows:
$ sbatch Run_Infomap.sh /path/to/Infomap_params.m

Where Infomap_params.m (see the template file) sets the following variables:
  1. subject (i.e., 'INET002')
  2. datafile (i.e., path to the same type of .xlsx file used by FDcalc, FCprocess, etc.)
  3. ciftiDir (folder containing the processed CIFTI data; ex. '/projects/b1081/iNetworks/Nifti/derivatives/postFCproc_CIFTI_20.2.0')
  	NOTE: based on the runs included in the datalist, the script will grab those corresponding CIFTIs:
	ciftiDir/sub-{subject}/ses-{session}/cifti_timeseries_normalwall/sub-{subject}_ses-{session}_task-{task}_run-{run}_LR_surf_subcort_222_32k_fsLR_smooth2.55.dtseries.nii
  4. outDir (will create if not already made; makes a subject-named folder inside)
  5. dmatname (path to a distance matrix to be used for distance exclusion)
  	NOTE: the code is set up to run Infomap on the surface (59,412 cortical surface vertices). The dmat in the template params file 'Cifti_geo_distances_xhemisphere_large.mat' is a 66697x66697 matrix; if the dmat is larger than the data correlation matrix (59412x59412), the script will assume that it should crop the distance matrix to fit the data.
  6. other necessary Infomap parameters:
  	xdistance (exclusion distance in mm; Infomap doesn't consider correlations of nearby vertices. Typical value for vertexwise (surface) data is 30)
 	thresholdarray (keeping this at [.003 .004 .005:.005:.05] for now; this seems to cut the time down to ~12 hours and produces a similar result)
	makebinary (0)
	numpools (currently set at 8 parallel pools that MATLAB will use to run Infomap across thresholds)
	structure_indices (kept for now as an empty matrix, but if you want to include subcortical structures, would use this variable to label the nodes)
	cortexOnly (1; currently set up this way)

OUTPUTS
  1. a time-stamped .mat file in the subject's output folder with a log of input parameters
  2. thresholds.txt - text file containing the thresholds 
  3. rawassn.txt - text file containing a community assignment for each vertex across all thresholds (e.g., if 12 thresholds, size would be 59412x12)

