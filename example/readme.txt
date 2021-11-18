1) in command prompt change folder to Shetran folder 
e.g cd C:\Users\steve\Documents\openclim\shetran-automaticcalibration\example\shetran
2) type shetran-automaticcalibration.exe 38017

This reads ../38017/38017_LibraryFile.xml and ../38017/optimise.csv

In optimise.csv there are the following 9 lines:

optimise parameters
FlowData1990.01.01-2001.12.31.csv
Number_of_spin_up_values,730     
deep_soil_conductivity,0.0001,100
shallow_soil_conductivity,1,100
shallow_soil_depth,0.5,4
AePe_ratio,0.5,2
Strickler,0.2,5
initial_water_table,0,10

The "Number_of_spin_up_values" igonres the first 730 discharge data values in the caculation of the NSE

Lines 4-9 are the minimum and maximum values for the 6 paramters that are calibrated. 

A sinlge soil cateogry file "soil-1type.asc" and single land-use file "landcover-1type.asc" is produced and used.

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
integer ,parameter       :: NoP=2,NoN=6,NoM=2*NoN+1,NoQ=NoN+1,NoBeta=2*NoN+1,NoS=NoP*NoM,noptvalues=2,NoIter=20


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

 