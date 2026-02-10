# 14C Antarctic Production Rate Calculator

Antarctic production rate calibration calculator v1.0

This code is structured as more of a notebook file than an adaptable code however the same methods may be applied in other applications using the same outline that this method applies.

This code was designed specifically to calibrate the production rate of in-situ 14C in Antarctica using the CRONUS-A intercomparison and reference sample as calibration data. It relies on the many measurements that have been made, and based on the fact that CRONUS-A is saturated with respect to in-situ 14C.

The calibration uses code developed by Greg Balco to calculate the production rate of the sample using Stone Scaling techniques and code written by Allie Koester and Nat Lifton, which calculates the LSDn scaled production rate.

CRONUS-A data is compiled from all published measured concentrations and used as a compilation to determine the most likely true concentration. 


## Within the main branch there are 4 directories: calibration, expage_v202403, age_offset, and figures.




#### Calibration performs the bulk of the work.

main14C_prodcalib.m is the main script. Additional files are needed for the production rate calibration but are not meant to be run. The script will run all the way through and produce a number of figures throughout that are used in the analysis and some of which are included in the paper.

Input files:
	all_Antarctic_14C.txt (All Antarctic 14C Data)
	CRONUSA.mat (All compiled CRONUS-A data)
	LSDn_inputs.txt (CRONUS-A data to input for the production rate)
Output files:
	At multiple locations in the code lines can be uncommented to save figures to your local machine, they are all commented in the current version to prevent many files from being written.
Outputs:
	The main production rate calibration is stored in the LSD and St structures which are written during the code run
	Additional outputs are stored in multiple variables with the other most likely referenced output is the CRONUS-A concentration which is stored in the CA structure. 






### expage_v202403 is the calculator which performs the age recalcualtion that we employ based on the expage calcualtor written by Jakob Heyman.
The expage calcualtor code was left unchanged except for the file which writes the constants including the reference production rate. We modify this specific file to alter the reference production rate.

Input files:
	input_lassitercoast.txt: Data for the lassiter coast samples
	input_TAMs: Data for the Northern Transantarctic Mountains
	input_ASE.txt: Data for the Amundsen sea nearby Pope glacier
	input_DML.txt: Data for Dronning Maud Land collected frm blue ice moraines there
	input_CRONUSA.txt: Calculates the age of CRONUS_A which we assume to be saturated, this is just a check
Outputs





### age_offset is used to explore the difference in the exposure ages that results from using different calcualtors and different production rates.
There is only a single code which is a matlab live script. Each code block is a different figure or multiple similar figures. Descriptions for each are included as markdown.







### Figures includes figures output from the code for reference

Additionally, the entire CRONUS-A compilation, and all of the production rate calibration data that are used in the paper are included as files in the repository for reference. The full production rate calibration data is taken directly from the supplement of Koester and Lifton, 2023 and is updated to include our CRONUS-A production rate calibration.



Written By Jason Drebber, 2026, Colorado School of Mines
jason_drebber@mines.edu
