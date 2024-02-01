1) in command prompt change folder to Shetran folder 
e.g cd C:\Users\steve\Documents\openclim\shetran-automaticcalibration\example\shetran
2) type shetran-automatic-calibration.exe 38017

This reads ../38017/38017_LibraryFile.xml and ../38017/optimise.csv

In optimise.csv there are the following 9 lines:

optimise parameters
NRFA_DailyFlows_38012_19800101-20101231.txt
Calibration_start_and_end_times,3654,7305
Validation_start_and_end_times,7306,10958
deep_soil_conductivity,0.0001,1.0
shallow_soil_conductivity,1,100
shallow_soil_depth,0.5,4
AePe_ratio,0.5,2.0
Urban_seperate_sewer_fraction,0.1,0.5


The "Calibration_start_and_end_times" are the times used in the caculation of the calibration NSE. In this case from 1/1/1990 to 31/12/1999  

The "Validation_start_and_end_times" are the times used in the caculation of the calibration NSE. In this case from 1/1/2000 to 31/12/2009  


Lines 5-9 are the minimum and maximum values for the 5 paramters that are calibrated. 

A two soil cateogry file "soil-2types.asc" and single land-use file "landcover-2types.asc" is produced and used. Category 2 is urban everything else is category 1. Urban areas have very low conductivities and high runoff Strickler coefficient. The urban category is half the urban fraction in the original land cover map (e.g. 38012_LandCover.asc) as the urban areas in this maps contains lots of gradens parks, etc. Also prooduced and used is "Urban_PET.csv" and "Urban_Precip.csv". For non-urban these are the same as before for urban areas the PET is set to zero and the precipitation depends on the "Urban_seperate_sewer_fraction". 

3) shetran-automaticcalibration.exe uses the SCE-UA global optimization method Duan (1994).

currenlty it carries out 546 simulations

! NoP = number of partitions
! NoN = number of points in each complex
! NoM = Number of members in each complex
! NoQ = Number of points in each sub-complex
! NoBeta = number of steps taken by each complex before the complexes are shuffled
! NoS = Number of sample points (= NoP * NoM)
! NoIter = Number of iterations of algorithm
! there are NoS initial runs. Then for each iteration NoP*Nobeta runs which are repeated NoIter times

!standard catchment
integer ,parameter       :: NoP=2,NoN=5,NoM=2*NoN+1,NoQ=NoN+1,NoBeta=2*NoN+1,NoS=NoP*NoM,noptvalues=2,NoIter=20
! big catchments with fewer runs
!integer ,parameter       :: NoP=1,NoN=5,NoM=2*NoN+1,NoQ=NoN+1,NoBeta=2*NoN+1,NoS=NoP*NoM,noptvalues=2,NoIter=10

!number of simulations
! if NoP=2,NoN=5,NoM=2*NoN+1,NoQ=NoN+1,NoBeta=2*NoN+1,NoS=NoP*NoM,noptvalues=2,NoIter=20
! NoS random = NoS = 2*11 = 22
! plus NoBeta* NoP * NoIter = 2 * 11 * 20 = 440
! total  = 462
!if  NoP=1,NoN=5,NoM=2*NoN+1,NoQ=NoN+1,NoBeta=2*NoN+1,NoS=NoP*NoM,noptvalues=2,NoIter=10
! NoS random = NoS = 1*11 = 11
! plus NoBeta* NoP * NoIter = 1 * 11 * 10 = 110
! total  = 121


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

 