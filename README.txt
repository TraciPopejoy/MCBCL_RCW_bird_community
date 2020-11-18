RCW as an umbrella species

Scripts:
U_load_data.R - loads in the two excel sheets that contains point count results from 2009 and 2010. Renames some of the column names for ease of use in the future. Filters the raw point count results to exclude rare birds (<1% observations) and birds > 300 m away. Creates the function sp_filter() that pulls one species from the point count results (DCERP_filt) and formats it to run within a ds model (distance sampling model from package Distance). It also creates a dataframe of all the sampled points and unique covariates at those points. Our table of habitat guilds is also generated here. 

U_run_models.R - code to run the distance sampling models in Distance. Output: one file for each species with all the run models information within it. Can use this information to rerun the ds function if needed. 

U_data_analysis.R - evaluates all the distance models run and identifies the best model. Best models are defined as those that match reality (GOF p value > 0.05) and have the lowest AIC. If no models match reality, we take the lowest AIC value from the five models with the highest GOF p value. Outputs: Table S3, Figure S1, and model.appendix.all.

U_quantify_bootstraps - evaluates the output of nonparametrically bootstrapping the distance and abundance models. ! Takes a while to run and a lot of memory ! Outputs: abundance.good (raw abundance used), abundance_summary (95% confidence at each habitat category for each bird), Figure 6, Table 3 information, Table 4  

U_evaluate_abundance - builds a pdf that depicts the median abundance of each species within each habitat group. Creates some tables to aid in writing the results. Runs the fuzzy analysis. Outputs: Figure 7 (fit of fuzzy ord), Figure 8 (bird community across sites), Figure 5 (total bird abundance), Table S4 & Figure S2 (detection diffs)

U_occupancy_models - runs occupancy models using the unmarked package. Outputs: Figure 3, Figure 4, Table S2, part of Table 4, Figure 2 for taxa richness