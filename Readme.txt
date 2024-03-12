080324
******

Note: set up for 50m of aquifer depth so need correct version of shetran-prepare. Also need correct version to change end of the rundata file to include addtional water-table depths.

Note: use latest version of Shetran.exe for outputting extra water level data

optimise.csv now includes waterlevel file and times for calibration and validation plus simulated file name


optimise parameters
NRFA_DailyFlows_33029_19800101-20101231.txt
CalibrationTimesFlow,3654,7305
ValidationTimesFlow,7306,10958
GreatThornsFarm.csv
CalibrationTimesGroundwater,113,230
ValidationTimesGroundwater,231,355
output_33029_psl-element133.txt
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




061023
******

for catchment 33029
integer:: calstart=113,calfin=230,valstart=231,valfinish=355


measuredWaterTable = '../'//trim(libraryfile)//'/GreatThornsFarm.csv'
simulatedWaterTable =  '../'//trim(libraryfile)//'/output_33029_psl-element133.txt'


   do z=1,NoDischargeValues
       if (MeasuredDischarge(z).GT.0) then
           logMeasured(z)=log10(MeasuredDischarge(z))
       else
           logMeasured(z)=-2
       endif
  enddo
   do z=1,NoDischargeValues
       if (SimulatedDischarge(z).GT.0) then
           logSimulated(z)=log10(SimulatedDischarge(z))
       else
           logSimulated(z)=-2
       endif

In measured discharge replace -999 with sensible values so that LogNSE works.

Change PBias so the aim is to achieve 0 not a high value
    Objfn(i)=NSETemp*NSE(i) + BiasTemp*(1/abs(Bias(i)))+ KGETemp*KGE(i) + LogNSETemp*LogNSE(i) + RMSETemp*(1/RMSE(i)) + PeaTemp*Pea(i) + SpeTemp*Spe(i)
Objfn=NSETemp*NSE + BiasTemp*(1/abs(Bias)) + KGETemp*KGE + LogNSETemp*LogNSE + RMSETemp*(1/RMSE) + PeaTemp*Pea + SpeTemp*Spe



090522
******
bug in code when there are lots of missing values. Calculation of mean values wrong

integer        :: nvalues



   MeanMeasured=0
   MeanSimulated=0
   nvalues=0
   do z=1, mindischargevalues-SpinUpValues
     if ((MeasuredDisShort(z)).GE.(0.0)) then
       MeanMeasured = MeanMeasured + MeasuredDisShort(z)
       MeanSimulated = MeanSimulated + SimulatedDisShort(z)
       nvalues=nvalues+1
     endif
   enddo
   MeanMeasured=MeanMeasured/nvalues
   MeanSimulated=MeanSimulated/nvalues





Set values to zero to avoid error 

Need to change linux version

  do z=1,mindischargevalues-SpinUpValues
     if ((MeasuredDisShort(z)).GE.(0.0)) then
       MeasuredSquare(z)=(MeasuredDisShort(z)-MeanMeasured)**2
       SimuMeasSquare(z)=(SimulatedDisShort(z)-MeasuredDisShort(z))**2
     else
        MeasuredSquare(z)=0
        SimuMeasSquare(z)=0
     endif
