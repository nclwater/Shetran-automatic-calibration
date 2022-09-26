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
