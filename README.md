# shetran-automatic-calibration-groundwater

## Features
- Executable to run automatic calibration of Shetran

## Usage
1) in command prompt change folder to Shetran folder 
e.g cd C:\Users\steve\Documents\openclim\shetran-automaticcalibration-groundwater\example\shetran
2) type shetran-automaticcalibration-groundwater.exe 33029

This reads ../33029/33029_LibraryFile.xml and ../33029/optimise.csv

In optimise.csv there are the following 20 lines:

optimise parameters
NRFA_DailyFlows_33029_19800101-20101231.txt
CalibrationTimesFlow,3654,7305
ValidationTimesFlow,7306,10958
GreatThornsFarm.csv
CalibrationTimesGroundwater,113,230
ValidationTimesGroundwater,231,355
output_WaterTable_Element133.txt
deep_soil_conductivity,0.1,200.0
shallow_soil_conductivity,10.0,200.0
shallow_soil_depth,0.5,5
AePe_ratio,1.0,1.5
UrbanPrecFraction,0.1,0.5
NSE,40
PBias,0
KGE,0
LogNSE,0
RMSE,20
Pearson,40
Spearman,0

"NRFA_DailyFlows_33029_19800101-20101231.txt" is the mneasured flow data

"CalibrationTimesFlow" and "ValidationTimesFlow" are the row numbers in the measured flow data for calibration and validation 

"GreatThornsFarm.csv" is the measured water level data. The formatting needs tobe the same as this one

"CalibrationTimesGroundwater" and "ValidationTimesGroundwater" are the row numbers in the measured water level data for calibration and validation

"output_WaterTable_Element133.txt" is the simulated water level data from the simulation start time. The element number 133 must correspond with the location of the measured data

Lines 9-13 are the minimum and maximum values for the 5 paramters that are calibrated.

Line 14-20 are the objective function amount. 14-17 are for flow data. 18-20 are for groundwater data


A two soil cateogry file "soil-2types.asc" and single land-use file "landcover-2types.asc" is produced and used. Category 2 is urban everything else is category 1. Urban areas have very low conductivities and high runoff Strickler coefficient. The urban category is half the urban fraction in the original land cover map (e.g. 38012_LandCover.asc) as the urban areas in this maps contains lots of gradens parks, etc. Also prooduced and used is "Urban_PET.csv" and "Urban_Precip.csv". For non-urban these are the same as before for urban areas the PET is set to zero and the precipitation depends on the "Urban_seperate_sewer_fraction".

3) shetran-automaticcalibration.exe uses the SCE-UA global optimization method Duan (1994).

! NoP = number of partitions
! NoN = number of points in each complex
! NoM = Number of members in each complex
! NoQ = Number of points in each sub-complex
! NoBeta = number of steps taken by each complex before the complexes are shuffled
! NoS = Number of sample points (= NoP * NoM)
! NoIter = Number of iterations of algorithm
! there are NoS initial runs. Then for each iteration NoP*Nobeta runs which are repeated NoIter times

!number of simulations
! if NoP=2,NoN=5,NoM=2*NoN+1,NoQ=NoN+1,NoBeta=2*NoN+1,NoS=NoP*NoM,noptvalues=2,NoIter=20
! NoS random = NoS = 2*11 = 22
! plus NoBeta* NoP * NoIter = 2 * 11 * 20 = 440
! total  = 462


For each simulation a new library file is produced.

LibraryFile1.xml

4) Each simulation runs:
shetran-prepare-snow.exe
and
shetran.exe

As usual shetran-prepare-snow.exe produces the input*** files and shetran.exe the output****

5) at the end of the simulation shetran-automaticcalibration.exe calculates the NSE and then appends

results.csv  - a line with the results and the parameters addded

results-complex.csv - details of the results for all the complex members added

6) when all the simulations are finished the best is run again.
