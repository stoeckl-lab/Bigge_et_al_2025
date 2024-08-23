# Bigge_et_al_2024
Analysis files for the publication Bigge et al. (2024) Integration of parallel pathways for flight control in a hawkmoth reflects prevalence and relevance of natural visual cues [LINK]
The raw data for this publication is deposited at https://figshare.com/s/e680da3be83fe172a5e4.

The scripts and functions in this repository build two analysis pipelines. One to extract and process flight tracks from a tunnel setup for hawkmoths, and one to analyse videos of natural visual scenes.  

# Flight track extraction and processing
1) autoTrack_tunnel_dorsal reads the raw videos that filmed the hawkmoth flights through the tunnel, and a) automatically extracts the flight coordinates and area the hawkmoth subtends in the video using an image-difference-based method, or b) re-routs of DLTdv6 for manual tracking https://biomech.web.unc.edu/dltdv/

2) analyse_tunnel_tracks_dorsal extracts multiple parameters from the flight tracks, such as speed, in-flight variation of position, position in the tunnel, and tracking area of the hawkmoth, for one experimental condition, and generates a comp_measures.mat file, which contains these parameters, as well as the original tracks. This file is deposited in the figshare repository for each condition, and is used for subsequent analysis across conditions.

3) tunnelStats_across_conditions_dorsoventral loads the comp_measures.mat file upon custom-selection (hard-coded lists at the top of the script, or manual folder selection after prompt) and plots the different flight track parameters across conditions. It also performs a statistical analysis with either ANOVA or Kruskal-Wallis test, depending on the normality of the ANOVA-test residuals (assessed with "lillietest") and saves the results as stats.txt files, indicating which test was performed.

4) plot_all_tracks plots all selected flight tracks within the tunnel dimensions


# Analysis of natural visual scenes

1) ImageDifferenceAnalysis computes the angle and magnitude of translational optic flow from videos of natural visual scenes (using the Methods described here: 10.1016/j.cub.2021.02.022). This script generates two .mat files for each video, one with the optic flow angles for each video (anglestack) and one with the mangnitudes corresponding to the angles (magnistack)

2) final_Analysis_new extracts the magnitude of translational optic flow for each quadrant of the visual field (dorsal, ventral, lateral left and right) for each video, and generates a summary file (filename.mat)
3) final_Analysis_contrast magnitude of contrast edges for each quadrant of the visual field (dorsal, ventral, lateral left and right) for each video, and generates a summary file (filename_contrast.mat)

4) summaryAnalysis generates sumamry heat maps for optic flow and contrast, respectively, from all selected .mat files generated in the previous step. They also plot boxplots of the resulting optic flow and contrast values across all habitat types and quadrants of the visual field. 

5) summaryAnalysis_Tunnel generates summary heat maps for optic flow and contrast, respectively, from all selected .mat files generated in the previous step for data filmed in the flight tunnel. It generates stacked bar plots of the resulting optic flow and contrast values in the quadrants of the visual field across all tunnel conditions.

6) GLMM_Mean_Test_DV_Imaging.R statistical analysis of the optic flow and contrast magnitudes across different habitat types in the quadrants of the visual field, using a linear mixed effects model


# Dependencies 
1) euclid_dist computes the Euclidian distance (D) and the angle relative to the horizontal (alpha) between points in a 2 column vector (if one input) or sets of points a and b if two inputs rows are individual observations, columns are x and then y value of the point
2) shadedErrorBar_anna, which is a slightly altered version of shadedErrorBar by Rob Campbell https://de.mathworks.com/matlabcentral/fileexchange/26311-raacampbell-shadederrorbar (which would also work)
5) Violin Plot Version 1.7.0.0 by Holger Hoffmann https://de.mathworks.com/matlabcentral/fileexchange/45134-violin-plot

