This is an instruction of generating Tables 3-4, Figures 1-10 in the paper "Comparing Bayesian Models with Spatially Varying Coefficients When Count Data Are Over-dispersed and Spatially Discontinuous" submitted to IJGIS. 

1. For real data application:     
The code is "real_data_code_upload.R", data used are "chicago_boundary_count_model_2023.csv" and "same_side_matrix.csv".     
To generate plots of Chicago, additional files used are "chicago_community_areas_current.geojson" and "chicago_police_districts_current.geojson".    

Step 1. install and library R packages     
Setp 2. run SECTION 1: load data     
Setp 3. run SECTION 2-SECTION 6: model preparation     
Step 4. run SECTION 7: MODEL1-Bayesian MGWNBR    
Step 5. run SECTION 8: MODEL 2-SVC-NBR     
Step 6. run SECTION 9: outputs     
---- generate Table 4 for model comparison     
---- generate Figure 7 (a)(d): estimates for \beta_0(s)     
---- generate Figure 7 (b)(e): estimates for \beta_1(s)     
---- generate Figure 7 (c)(f): estimates for \beta_2(s)     
---- generate Figure 8: Absolute difference between MGWNBR and SVC-NBR posterior mean intercept estimates     
---- generate Figure 6: Spatial distribution of crime counts     
---- generate Figure 10: Boundary-distance profiles of absolute coefficient discrepancies     
---- generate Figure 9: Relationship between absolute intercept discrepancy j∆β0(s)j and distance     
     
     
2. For simulation study:     
The code is "simulation_upload.R".     
     
Step 1. run Section 1: Simulation Setpip     
---- generate simulated datasets for different scenarios     
Step 2. run Section 2: run the simulations via HPC     
---- generate the simulated results with R=200 replications     
Step 3. run Section 3: collect simulation results and generate plots     
---- collect the 200 simulation results files, summarise and then obtain plots     
---- generate Figure 1: True coefficient surfaces under different simulation settings     
---- generate Figure 2: Global coefficient recovery measured by the mean RMSE ratio     
---- generate Table 3: Median information-criterion differences     
---- generate Figure 3: Information-criterion comparison between MGWNBR and SVC-NBR     
---- generate Figure 4: Boundary-distance RMSE ratio under hard-discontinuity scenarios     
---- generate Figure 5: Boundary-distance win probability for MGWNBR     


