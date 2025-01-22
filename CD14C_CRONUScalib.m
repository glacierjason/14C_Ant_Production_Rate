function out = CD14C_CRONUScalib(sampledata, flag, PR)
% wrapper script to calculate the compositionally dependent in situ 14C production rate for a sample
%
% syntax = CD14C('sampledata.txt')
% output is a binary .mat file

% Written by Allie Koester and Nat Lifton 2022, Purdue University
% koestea@purdue.edu

% Modified by Jason Drebber specifically for the purpose of calculating
% production rates for quartz while at the Colorado School of Mines, 2025
% jason_drebber@mines.edu

% Based on code by Greg Balco (Balco et al., 2008)
% April, 2007
% Part of the CRONUS-Earth online calculators: 
%      http://hess.ess.washington.edu/math
% Copyright 2001-2007, University of Washington
% All rights reserved
% Developed in part with funding from the National Science Foundation.

% Also based on code by Nat Lifton (Lifton et al., 2014; Lifton, 2016);
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 2,
% as published by the Free Software Foundation (www.fsf.org).

version = '1.0 - 05/2022';

% If no production rate is passed as an argument the code will assume that
% a production rate calibration is being performed for in-situ 14C in
% quartz and will assign a linear range of possible production rate values
% to the production rate. 
if nargin == 2
    P14.min = 10; %minimum estimate
    P14.max = 20; %maximum estimate
    P14.intervals = (P14.max-P14.min)*1000; %number of interval values to check
    P14.range = linspace(P14.min, P14.max, P14.intervals);
elseif nargin == 3
    P14.calib = PR;
end

%% Import data

%Input data format is Sample_Name; Latitude (deg N); Longitude (deg E); Elevation (m); 
% SiO2, TiO2, Al2O3, Fe2O3, FeO, MnO, MgO, CaO, Na2O, K2O, and P2O5 all in
% weight percent
if flag == 1
    FID = fopen(sampledata);
    data = textscan(FID,'%s %n %n %n %n %n %n %n %n %n %n %n %n %n %n');
    fclose(FID);
    dstring='';

    %Import the sample data vectors
    all_sample_name = data{1}; %sample name
    all_lat = data{2}; %latitude
    all_long = data{3}; %longitude
    all_elv = data{4}; %elevation
    all_xrf = cell2mat(data(5:end)); %matrix of all the sample major element data

    %calculate elemental number density for each sample 
    all_ND = numdCD(all_xrf);

else
    data = sampledata;
    all_sample_name = data(:,1);
    all_lat = data(:,2); %latitude
    all_long = data(:,3); %longitude
    all_elv = data(:,4); %elevation
    all_xrf = cell2mat(data(:,5:end)); %matrix of all the sample major element data


    %calculate elemental number density for each sample 
    all_ND = numdCD(all_xrf);
end

%% Load the constants file
load consts_CD14C.mat;

%% Scaling and Production calculations 
if flag == 1
    sample.sample_name = all_sample_name;
    sample.lat = all_lat;
    sample.long = all_long;
    sample.elv = all_elv;
    sample.ND = all_ND;
    sample.pressure = ERA40atm(sample.lat,sample.long,sample.elv); %Pressure in hPa

    scaling14 = ScalingLSD_CD(sample,consts);


    for i=1:length(P14.range)

        %Calculate the spallation production rate using gridded Rc values
        P14_CD(i,:) = scaling14.SF_LS_CD.*P14.range(i);  % Equation 3
                                    %LSD calibrated qtz spallation PR
    end
else
    num_samples = length(all_lat);
    for a = 1:num_samples
        sample.sample_name = all_sample_name(a);
        sample.lat = cell2mat(all_lat(a));
        sample.long = cell2mat(all_long(a));
        sample.elv = cell2mat(all_elv(a));
        sample.ND = all_ND(a,:);

        sample.pressure = ERA40atm(sample.lat,sample.long,sample.elv);


        
        scaling14 = ScalingLSD_CD(sample,consts);
    
        %Calculate the spallation production rate using gridded Rc values
        P14_CD(a,:) = scaling14.SF_LS_CD.*P14.calib;  % Equation 3
                                    %LSD calibrated qtz spallation PR
     

    end
end


%% output the results
%save output data to workspace
out.ID = all_sample_name';
out.tv = scaling14.tv; %time vector, relative to 2010 (t = 0)
out.P14_CD = P14_CD; %time- and compositionally dependent site production 
% rate relative to that geologically calibrated for quartz
%out.P14gd_CD = P14gd_CD; %time- and compositionallydependent site production 
% rate relative to that geologically calibrated for quartz, for geocentric dipole


%save(replace(sampledata,".txt","results.mat"),'out');
%disp(['Constants version ' consts.version]);
%disp('Saved'); 

end