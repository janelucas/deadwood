# deadwood
Suspended sections within downed dead wood are drier, have altered decomposer communities, and slower decomposition

# Title of Dataset
ReadMe file for data submitted as part of Barrera et al. 2023 publication entitled 'Suspended sections within downed dead wood are drier, have altered decomposer 
2 communities, and slower decomposition' in Ecosystems Journal

## Description of the data and file structure
Experiment	Files	Columns	Detailed descriptions
Wood decomposition survey	ws_50ha_susp_decomp_1to5.csv	percentage_suspended	0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100. What percent of total dead wood volume is not  touching the soil, which may include a piece of trunk that is visible standing (i.e., 100). 
Wood decomposition survey	ws_50ha_susp_decomp_1to5.csv	volume_remaining	Percentage of total volume remaining (i.e., 0-100%)
Wood decomposition survey	tri2tu_all_remaining.csv	percentage_suspended	0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100. What percent of total dead wood volume is not  touching the soil, which may include a tree trunk that is standing (i.e., 100% suspended). 
Wood decomposition survey	tri2tu_all_remaining.csv	source	Categorical term delineating whether the tree was marked dead in a 5-year census (i.e., plot.census) interval or annual census interval (i.e., ws_50ha_nona_elev1_tri2tu_trimdendro)
Wood decomposition survey	tri2tu_all_remaining.csv	years_since_death	Number of years since tree died
Wood decomposition survey	tri2tu_all_remaining.csv	Volume loss (%)	Percentage of total volume lost 
Wood decomposition survey	tri2tu_by_years_since_death.csv	source	Categorical term delineating whether the tree was marked dead in a 5-year census (i.e., plot.census) interval or annual census interval (i.e., ws_50ha_nona_elev1_tri2tu_trimdendro)
Wood decomposition survey	tri2tu_by_years_since_death.csv	years_since_death	Number of years since tree died
Wood decomposition survey	tri2tu_by_years_since_death.csv	mean_percentage_suspended	Mean of the percent of the total dead wood volume that is not touching the soil, which may include a tree trunk that is standing (i.e., 100%). 
Wood decomposition survey	tri2tu_by_years_since_death.csv	sample_size	Number of Trichilia tuberculata trees that died that year (years_since_death)
Wood decomposition survey	tri2tu_by_years_since_death.csv	Volume loss (%)	Percentage of volume lost 
Wood respiration and macroorganism measurements	respiration_model_results.csv	wood_ID	Number corresponding to the tree in order of collection (1,2,3...)
Wood respiration and macroorganism measurements	respiration_model_results.csv	treatment	Part of the trunk from which the wood sample was collected. gu (ground up), gd (ground down), eu (elevated up), ed (elevated down)
Wood respiration and macroorganism measurements	respiration_model_results.csv	decay_status	Wood decay scale(1-5) according to Larjavaara&Muller-Landau (2010)
Wood respiration and macroorganism measurements	respiration_model_results.csv	moisture_content	Moisture content of wood (%) resulted from sustracting wet weight (initial_weight) minus dry weight (final_weight), divided by wet weight (initial_weight) and multiplied by 100
Wood respiration and macroorganism measurements	respiration_model_results.csv	wood_density	Ratio between the dry weight of wood sample (final_weight) divided by the volume of the wood sample (wood_volume) (g/cm^3)
Wood respiration and macroorganism measurements	respiration_model_results.csv	slope	Respiration rate (ppm CO2 s−1)
Wood respiration and macroorganism measurements	respiration_model_results.csv	visible_algae_lichen_num	Presence of photosynthetic growth in the wood sample collected (yes=1, no=0)
Wood respiration and macroorganism measurements	respiration_model_results.csv	visible_hyphae_num	Presence of  visible hyphae (e.g., hyphal strands, fans, rhizomorphs) in the wood sample collected (yes=1, no=0)
Wood respiration and macroorganism measurements	respiration_model_results.csv	termites_num	Presence of termites or termite galleries in the wood sample collected (yes=1, no=0)
Wood respiration and macroorganism measurements	respiration_model_results.csv	moist_cont_dry	Moisture content of wood (%) resulted from sustracting wet weight (initial_weight) minus dry weight (final_weight), divided by dry weight (final_weight) and multiplied by 100
Experimental manipulations of fine wood connectedness and ground contact	mass_loss.csv	group	Experimental subset to which the block belonged
Experimental manipulations of fine wood connectedness and ground contact	mass_loss.csv	block	Main block of wood from which that block was obtained
Experimental manipulations of fine wood connectedness and ground contact	mass_loss.csv	mass.remaining	Amount of mass in grams that remained after the experiment, after 16 months in the field
Experimental manipulations of fine wood connectedness and ground contact	mass_loss.csv	updown	Wood suspension (downed or suspended)
Experimental manipulations of fine wood connectedness and ground contact	mass_loss.csv	connection	Wood connectivity (suspended sections either connected to or separated from a downed section of the same wood piece) 
Experimental manipulations of fine wood connectedness and ground contact	mass_loss.csv	mass.loss	Amount of mass in grams that was decomposed after the experiment, after 16 months in the field
Experimental manipulations of fine wood connectedness and ground contact	mass_loss.csv	logit.mass.loss	Logit transformed mass loss data (`mass.loss`)![image](https://github.com/janelucas/deadwood/assets/50459345/506e5534-0e1a-48ee-9ad6-ec40fdb0a3bb)

## Sharing/Access information
NA
## Code/Software
Included and described above
