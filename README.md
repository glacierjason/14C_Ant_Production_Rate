# 14C Antarctic Production Rate Calculator

Antarctic production rate calibration calculator v1.0

This code is structured as more of a notebook file than an adaptable code however the same methods may be applied in other applications using the same outline proposed here.

This code was designed specifically to calibrate the production rate of in-situ 14C in Antarctica using the CRONUS-A intercomparison and reference sample as calibration data.

The calibration uses code developed by Greg Balco to calculate the production rate of the sample using Stone Scaling techniques and made by Allie Koester and Nat Lifton, which calculates the LSDn scaled production rate.

CRONUS-A data is compiled from all published measured concentrations and used as a compilation to determine the most likely true concentration. 


main14C_prodcalib.m is the main script. 

Input files:
	all_Antarctic_14C.txt (All Antarctic 14C Data)
	CRONUSA.mat (All compiled CRONUS-A data)
	LSDn_inputs.txt (CRONUS-A data to input for the production rate)


Written By Jason Drebber, 2025, Colorado School of Mines
jason_drebber@mines.edu
