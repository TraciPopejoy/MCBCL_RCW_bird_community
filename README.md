# RCW as an umbrella species in Camp LeJeune

## Scripts
1. **U_load_data.R**
Loads in the two excel sheets that contains point count results from 2009 and 2010. Renames some of the column names for ease of use in the future. Filters the raw point count results to exclude rare birds (<1% observations) and birds > 300 m away. Creates the function sp_filter() that pulls one species from the point count results (DCERP_filt) and formats it to run within a ds model (distance sampling model from package Distance). It also creates a dataframe of all the sampled points and unique covariates at those points. Our table of habitat guilds is also generated here. 

2. **U_run_models.R** 
Code to run the distance sampling models in Distance. This creates one file for each species that contains all the model information for all distance models run. Can use this information to rerun the ds function if needed. 

3. **U_dist_mod_analysis.R** 
Evaluates all the distance models run and identifies the best model. Best models are defined as those that match reality (GOF p value > 0.05) and have the lowest AIC. If no models match reality, we take the lowest AIC value from the five models with the highest GOF p value. 
- Table S3: description of Detection Probability models from distance sampling
- Figure S2: graph of count of model parameters from the top models
- model.appendix.all: table of top distance sampling models

4. **U_model_abundance.R**
Using the top distance sampling model for each species to estimate detection probability, we then completed non-parametric bootstraps to estimate abundance at each point and then relate that to RCW habitat quality. This creates two files for each day and species: one that contains model coefficients and another that contains raw abundance estimates for each point from each model.

5. **U_simp_occ.R**
Runs occupancy models using the unmarked package.
- Figure 3: RCW occupancy and abund
- Figure 4: BACS occupancy and abund
- Table S2: Occupancy detection portions
- Figure 2: taxa richness
- Figure S1: checker plot to quickly summarize findings

6. **U_quantify_bootstraps** 
Evaluates the output of nonparametrically bootstrapping the distance and abundance models. ! Takes a while to run and a lot of memory ! 
- abundance.good: raw abundance used
- abundance_summary: 95% confidence at each habitat category for each bird
- Figure 6: Coefficients between RCW habitat and occupancy and density
- Table 3: sites and points used
- Table 4: detection, RCW habitat affect on density, density averages  

7. **U_evaluate_abundance** 
Builds a pdf that depicts the median abundance of each species within each habitat group. Creates some tables to aid in writing the results. Runs the fuzzy analysis. 
- Figure 7: fit of fuzzy ord
- Figure 8: bird community across sites
- Figure 5: total bird abundance
- Table S4 & Figure S2: detection diffs