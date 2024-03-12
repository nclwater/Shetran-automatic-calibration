!  shetran-automaticcalibration.f90 
!
!  FUNCTIONS:
!  Shetranoptimise - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Shetranoptimise
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

program ShetranAutoCalibration


implicit none

! Variables

integer        :: n1,i,res,namelength,j,k,l,m,n,x,y,z,simnumber=0,simbest
integer        :: io,filelength,vegend,soilpropend,soildetend,ncols,nrows
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

integer        :: NoDischargeValues,itemp
integer        :: SpinUpValues

!integer(1), dimension(:),allocatable            :: arrayof1s
integer, dimension(:,:),allocatable            :: vegtype,vegtypenew
character(300) :: optimisefull,optimisefilename,librarytext,resultsfile,libraryfilepath,dischargefilepath,resultsfile2
character(300) :: GroundwaterLevelFile,GroundLevelFilePath,SimulatedLevelFile,SimulatedLevelFilepath
character(300) :: title,libraryfile,discharge,text1,catchmentname,dischargefilename,vegmap,soilmap,vegfilename,soilfilename
character(300) :: vegfilename1type,soilfilename1type,textasc(6),rundatafilename,text2,text3,text4
character(300), dimension(:),allocatable  :: text
character(300)  :: simnumberstring,libraryfileoutput



real            :: random,temp,temparray(NoS)
real            :: DeepSoilCond(noptvalues),ShallowSoilCond(noptvalues)
real            :: ShallowSoilDepth(noptvalues),AePeRatio(noptvalues)
real            :: Urban(noptvalues)
real            :: Strickler(noptvalues),InitWaterTable(noptvalues)
real            :: DeepSoilCond1(NoS),ShallowSoilCond1(NoS)
real            :: ShallowSoilDepth1(NoS),AePeRatio1(NoS)
real            :: Urban1(NoS)
real            :: prob(NoS),rtemp1,rtemp2,ProbTimesRandom(NoS)
real, dimension(:),allocatable            :: MeasuredDischarge
real                                      :: BiasCal(NoS),NSECal(NoS),KGECal(NoS),LogNSECal(NoS)
real                                      :: BiasVal(NoS),NSEVal(NoS),KGEVal(NoS),LogNSEval(NoS)
real                                      :: RMSECal(NoS),PeaCal(NoS),SpeCal(NoS)
real                                      :: RMSEVal(NoS),PeaVal(NoS),SpeVal(NoS)
real                                      :: ObjFn(NoS),ObjFnTemp,NSETempArray(NoS),ObjFnBest=-9999
real                                      :: BiasCalTemp,NSECalTemp,KGECalTemp,LogNSECalTemp,BiasValTemp,NSEvalTemp,KGEValTemp,LogNSEValTemp
real                                      :: RMSECalTemp,peaCalTemp,SpeCalTemp,RMSEValTemp,PeavalTemp,SpeValTemp
real                                      ::DeepSoilCondTemp,ShallowSoilCondTemp,ShallowSoilDepthTemp,AePeRatioTemp,UrbanTemp
real                                      ::DeepSoilCondTemp1,ShallowSoilCondTemp1,ShallowSoilDepthTemp1,AePeRatioTemp1,UrbanTemp1
real                                      ::DeepSoilCondTemp2,ShallowSoilCondTemp2,ShallowSoilDepthTemp2,AePeRatioTemp2,UrbanTemp2
real                                      ::NSECalComplex(NoM,NoP),BiasCalComplex(NoM,NoP),KGECalComplex(NoM,NoP),LogNSECalComplex(NoM,NoP)
real                                      ::NSEValComplex(NoM,NoP),BiasValComplex(NoM,NoP),KGEValComplex(NoM,NoP),LogNSEValComplex(NoM,NoP)
real                                      ::RMSECalComplex(NoM,NoP),PeaCalComplex(NoM,NoP),SpeCalComplex(NoM,NoP)
real                                      ::RMSEValComplex(NoM,NoP),PeaValComplex(NoM,NoP),SpeValComplex(NoM,NoP)
real                                      ::DeepSoilCondComplex(NoM,NoP),ShallowSoilCondComplex(NoM,NoP)
real                                      ::ShallowSoilDepthComplex(NoM,NoP),AePeRatioComplex(NoM,NoP),UrbanComplex(NoM,NoP)
real                                      ::NSEFrac,BiasFrac,KGEFrac,LogNSEFrac,RMSEFrac,PeaFrac,SpeFrac
logical        ::outside1(NoN),outside2(NoN)
logical        :: first
integer        :: calibrationTimes(2),validationTimes(2)
integer        :: calstart=113,calfin=230,valstart=231,valfinish=355








! Body of Shetranoptim
print*
print*, 'Shetran automatic calibration'
print*, '*****************************'
print* 
      
! read n1 from the command line- catchmentnumber
n1=1
CALL GETARG(n1,libraryfile)

!libraryfile='33029'
!libraryfile='UCauca5000_easy'

!do not change the optimise.csv file here. But will need to change 
!Number_of_spin_up_values,730
!set calibration from 1/1/1990 to 31/12/1999 (3654 to 7305)
!validation from 1/1/2000 to 31/12/2009 (7306 to 10958 )
! Strickler,0.2,1 take value in libraryfile.xml
! initial_water_table,0,10, hard code to 4.0
! urban rain fraction hard code to 0.1-0.5



!optimise.csv is the list of parameter values
optimisefull = '../'//trim(libraryfile)//'/optimise.csv'
!LibraryFile.xml is the original libbrary file
libraryfilepath = '../'//trim(libraryfile)//'/'//trim(libraryfile)//'_LibraryFile.xml'
open(10,FILE=optimisefull,err=9999,status='old')
open(11,FILE=trim(libraryfilepath),err=9998,status='old')



!output files
!***********
!results.csv is the model results from all the models
resultsfile= '../'//trim(libraryfile)//'/results.csv'
resultsfile2= '../'//trim(libraryfile)//'/results-complex.csv'
open(15,FILE=resultsfile)
open(19,FILE=resultsfile2)
write(15,'(A16)') 'Shetran-Optimise'
write(15,'(A16)') '****************'
write(15,'(A216)') 'number,deep_soil_conductivity,shallow_soil_conductivity,shallow_soil_depth,AePe_ratio,urban-prec-fraction,NSECal,BiasCal,KGECal,LogNSECal,RMSECal,PeaCal,SpeCal,NSEVal,BiasVal,KGEVal,LogNSEVal,RMSEVal,PeaVal,Speval'
write(19,'(A231)') 'iteration,complex,beta,deep_soil_conductivity,shallow_soil_conductivity,shallow_soil_depth,AePe_ratio,urban-prec-fraction,NSECal,BiasCal,KGECal,LogNSECal,RMSECal,PeaCal,SpeCal,NSEVal,BiasVal,KGEVal,LogNSEVal,RMSEVal,PeaVal,Speval'


read (10,*) title
!read (10,*) libraryfile

!Read measured flow data file and calibration and validation row numbers
read (10,*) discharge
dischargefilepath = '../'//trim(libraryfile)//'/'//discharge
open(12,FILE=trim(dischargefilepath),err=9997,status='old')
!print*, optimisefull,resultsfile,libraryfilepath,dischargefilepath
read(10,*) text1,calibrationTimes(1),calibrationTimes(2)
read(10,*) text1,validationTimes(1),validationTimes(2)

!Read Groundwater level file and calibration and validation row numbers
read (10,*) GroundwaterLevelFile
GroundLevelFilePath = '../'//trim(libraryfile)//'/'//GroundwaterLevelFile
read(10,*) text1,calstart,calfin
read(10,*) text1,valstart,valfinish

!Read Simulated Flow file. This needs the correct version of SHetran.exe to produce this fike
read (10,*) SimulatedLevelFile
SimulatedLevelFilepath = '../'//trim(libraryfile)//'/'//SimulatedLevelFile


!optimsed parameters
read(10,*) text1,(DeepSoilCond(i),i=1,noptvalues)
read(10,*) text1,(ShallowSoilCond(i),i=1,noptvalues)
read(10,*) text1,(ShallowSoilDepth(i),i=1,noptvalues)
read(10,*) text1,(AePeRatio(i),i=1,noptvalues)
read(10,*) text1,(Urban(i),i=1,noptvalues)

!optimsed objective functions
read(10,*) text1,NSEFrac
read(10,*) text1,BiasFrac
read(10,*) text1,KGEFrac
read(10,*) text1,LogNSEFrac
read(10,*) text1,RMSEFrac
read(10,*) text1,peaFrac
read(10,*) text1,SpeFrac



close(10)

!number of values in discharge file
NoDischargeValues = 0
read(12,*)
do 
  Read(12,'(A11,f10.4)',IOSTAT=io) Text1
  if (io < 0) then
    EXIT
  else
    NoDischargeValues= NoDischargeValues + 1
  endif
enddo
!define array sizes for discharge data
allocate (MeasuredDischarge(NoDischargeValues))
rewind(12)

!read measured discharge data
read(12,*)
do i=1,NoDischargeValues
  Read(12,'(A11,f10.0)') Text1,MeasuredDischarge(i) 
enddo
close(12)

!find length of xml library file
filelength=0
do 
read(11,'(A)',IOSTAT=io) text1
  if (io < 0) then
    EXIT
  else
    filelength= filelength + 1
  endif
enddo
allocate (text(filelength))
rewind(11)

!read xml library file and row location of key points
do i=1,filelength
   read(11,'(A)') text(i)
   if (trim(text(i)).eq.'</VegetationDetails>') then
      vegend=i
   endif
   if (trim(text(i)).eq.'</SoilProperties>') then
      soilpropend=i
   endif
   if (trim(text(i)).eq.'</SoilDetails>') then
      soildetend=i
   endif
enddo
close(11)
!print*,filelength,NoDischargeValues,vegend,soilpropend,soildetend

!shetran files names
text1=text(3)
namelength=len_trim(text1)
catchmentname=text1(16:namelength-16)
dischargefilename='../'//trim(libraryfile)//'/output_'//trim(catchmentname)//'_discharge_sim_regulartimestep.txt'
rundatafilename = '../'//trim(libraryfile)//'/rundata_'//trim(catchmentname)//'.txt'







!produce file with 2 land-use type. type 2 urban, type 1 everything else
! only every second urban
text1=text(7)
namelength=len_trim(text1)
vegmap=text1(9:namelength-9)
vegfilename='../'//trim(libraryfile)//'/'//trim(vegmap)
vegfilename1type='../'//trim(libraryfile)//'/'//'landcover-2types.asc'
open(16,FILE=trim(vegfilename),err=9995,status='old')
open(18,FILE=trim(vegfilename1type))
read (16,*) text1,ncols
write (18,'(A5,i4)') trim(text1),ncols
read (16,*) text1,nrows
write (18,'(A5,i4)') trim(text1),nrows
do z=1,4
    read (16,*) textasc(z)
    write (18,'(A)') textasc(z)
enddo
allocate (vegtype(nrows,ncols))
allocate (vegtypenew(nrows,ncols))
first=.true.
do i=1,nrows
   read (16,*) (vegtype(i,z),z=1,ncols)
   do j=1,ncols
       if (Vegtype(i,j).eq.1) then
           Vegtypenew(i,j)=1
       elseif (Vegtype(i,j).eq.2) then
           Vegtypenew(i,j)=1
       elseif (Vegtype(i,j).eq.3) then
           Vegtypenew(i,j)=1
       elseif (Vegtype(i,j).eq.4) then
           Vegtypenew(i,j)=1
       elseif (Vegtype(i,j).eq.5) then
           Vegtypenew(i,j)=1
       elseif (Vegtype(i,j).eq.6) then
           Vegtypenew(i,j)=1
       elseif (Vegtype(i,j).eq.7) then
           if (first) then
               Vegtypenew(i,j)=2
              first=.false.
           else
              Vegtypenew(i,j)=1
              first=.true.
          endif
               
       else
           Vegtypenew(i,j)=-9999           
       endif
    enddo
   write (18,'(*(I5,1X))') (Vegtypenew(i,z),z=1,ncols)
enddo


close(16)
close(18)

text1=text(8)
namelength=len_trim(text1)
soilfilename1type='../'//trim(libraryfile)//'/'//'soil-2types.asc'
open(16,FILE=trim(vegfilename),err=9995,status='old')
open(18,FILE=trim(soilfilename1type))
read (16,*) text1,ncols
write (18,'(A5,i4)') trim(text1),ncols
read (16,*) text1,nrows
write (18,'(A5,i4)') trim(text1),nrows
do z=1,4
    read (16,*) textasc(z)
    write (18,'(A)') textasc(z)
enddo
do i=1,nrows
   write (18,'(*(I5,1X))') (Vegtypenew(i,z),z=1,ncols)
enddo


close(16)
close(18)






!goto 754

!Run random points
!*****************
! Duan step (1)
!NoS random ponits
simnumber=0
do i=1,NoS
!    OrigSim(i)=i
    simnumber=simnumber+1
    call random_number(random)
    ! put deep soil conductivity over a log scale
!    temp=(log10(DeepSoilCond(2))-log10(DeepSoilCond(1)))*random+log10(DeepSoilCond(1))
!    print*,temp
!    DeepSoilCond1(i)=10.0**temp
!    print*,random(1),temp,DeepSoilCond1(i)
    DeepSoilCond1(i)=(DeepSoilCond(2)-DeepSoilCond(1))*random+DeepSoilCond(1)
    call random_number(random)
    ShallowSoilCond1(i)=(ShallowSoilCond(2)-ShallowSoilCond(1))*random+ShallowSoilCond(1)
    call random_number(random)
    ShallowSoilDepth1(i)=(ShallowSoilDepth(2)-ShallowSoilDepth(1))*random+ShallowSoilDepth(1)
    call random_number(random)
    AePeRatio1(i)=(AePeRatio(2)-AePeRatio(1))*random+AePeRatio(1)

!     call random_number(random)
!   Strickler1(i)=(Strickler(2)-Strickler(1))*random+Strickler(1)

    call random_number(random)
   Urban1(i)=(Urban(2)-Urban(1))*random+Urban(1)

!write XML file, run shetran and work out NSE
    Call WriteXmlFile(dischargefilename,rundatafilename,libraryfile,filelength,text,simnumber,vegend,soilpropend,soildetend,NoDischargeValues,calibrationTimes,validationTimes,NoS,MeasuredDischarge, &
      AePeRatio1(i),ShallowSoilCond1(i),DeepSoilCond1(i),ShallowSoilDepth1(i),Urban1(i),BiasCal(i),NSECal(i),KGECal(i),LogNSECal(i),RMSECal(i),PeaCal(i),SpeCal(i), &
        BiasVal(i),NSEval(i),KGEVal(i),LogNSEVal(i),RMSEVal(i),PeaVal(i),Speval(i),nrows,ncols,Vegtypenew,calstart,calfin,valstart,valfinish,GroundLevelFilePath,SimulatedLevelFilepath)

    
    
       write(15,'(i4,A1,18(f10.4,A1),f10.4)') Simnumber,',',DeepSoilCond1(i),',',ShallowSoilCond1(i),',',ShallowSoilDepth1(i),',',AePeRatio1(i),',',Urban1(i),',',NSECal(i),',',BiasCal(i),',',KGECal(i),',',LogNSECal(i),',',RMSECal(i),',',PeaCal(i),',',SpeCal(i),',',NSEVal(i),',',BiasVal(i),',',KGEVal(i),',',LogNSEVal(i),',',RMSEVal(i),',',PeaVal(i),',',SpeVal(i)

       

    
enddo

!temp
!754 DeepSoilCond1= (/ 0.0001,0.103,0.0003,28.13,0.0003 /)
!ShallowSoilCond1= (/ 3.52,91.62,88.94,10.67,91.27/)
!ShallowSoilDepth1= (/ 1.73,3.29,2.95,0.64,0.82/)
!AePeRatio1= (/ 1.5,1.75,1.6,0.63,1.46/)
!Strickler1= (/ 4.82,1.86,1.64,2.88,4.29/)
!Urban1= (/ 8.38,8.71,0.5,9.26,1.21/)
!NSE= (/-3.68,-8.69,-10.61,-16.35,-16.67/)
!Bias= (/-30.71,177.43,207.53,347.9,211.08/)
!OrigSim = (/1,2,3,4,5/)
!write(642,'(A124)') 'deep_soil_conductivity,shallow_soil_conductivity,shallow_soil_depth,AePe_ratio,Strickler,initial_water_table,NSE,Bias'
!write(643,'(A124)') 'deep_soil_conductivity,shallow_soil_conductivity,shallow_soil_depth,AePe_ratio,Strickler,initial_water_table,NSE,Bias'
!write(644,'(A124)') 'deep_soil_conductivity,shallow_soil_conductivity,shallow_soil_depth,AePe_ratio,Strickler,initial_water_table,NSE,Bias'
!write(645,'(A124)') 'deep_soil_conductivity,shallow_soil_conductivity,shallow_soil_depth,AePe_ratio,Strickler,initial_water_table,NSE,Bias'
!write(646,'(A124)') 'deep_soil_conductivity,shallow_soil_conductivity,shallow_soil_depth,AePe_ratio,Strickler,initial_water_table,NSE,Bias'

!do i=1,NoS
!  write(642,'(7(f10.4,A1),f10.4)') DeepSoilCond1(i),',',ShallowSoilCond1(i),',',ShallowSoilDepth1(i),',',AePeRatio1(i),',',Strickler1(i),',',Urban1(i),',',NSE(i),',',Bias(i)
!enddo
!  write(642,*)
!end temp

!number of iterations of each 
!do x=1,0 
do x=1,NoIter  
  
!*****************************
!Duan step (2) rank points  
    
Call objectivefunction(NoS,NSECal,BiasCal,KGECal,logNSECal,RMSECal,PeaCal,SpeCal,NSEFrac,BiasFrac,KGEFrac,LogNSEFrac,RMSEFrac,PeaFrac,SpeFrac,Objfn)
    
    
!ObjFn=NSECal
!NSEtempArray=NSE
!order by NSE
Call orderdata(NoS,ObjFn,DeepSoilCond1,ShallowSoilCond1,ShallowSoilDepth1,AePeRatio1,Urban1,NSECal,BiasCal,KGEcal,LogNSECal,RMSECal,PeaCal,SpeCal,NSEVal,BiasVal,KGEVal,LogNSEVal,RMSEVal,PeaVal,Speval)
 
!*****************************
!Duan step (3) partition into p complexes each containg m points
do i=1,NoM
    do j=1,NoP
          itemp=NoP*(i-1) + j
          DeepSoilCondComplex(i,j)=DeepSoilCond1(itemp)
          ShallowSoilCondComplex(i,j)=ShallowSoilCond1(itemp)
          ShallowSoilDepthComplex(i,j)=ShallowSoilDepth1(itemp)
          AePeRatioComplex(i,j)=AePeRatio1(itemp)
!          StricklerComplex(i,j)=Strickler1(itemp)
          UrbanComplex(i,j)=Urban1(itemp)
          NSECalComplex(i,j)=NSECal(itemp)
          BiasCalComplex(i,j)=BiasCal(itemp)
          KGECalComplex(i,j)=KGECal(itemp)
          LogNSECalComplex(i,j)=LogNSECal(itemp)
          RMSECalComplex(i,j)=RMSECal(itemp)
          PeaCalComplex(i,j)=PeaCal(itemp)
          SpeCalComplex(i,j)=SpeCal(itemp)
          NSEValComplex(i,j)=NSEVal(itemp)
          BiasValComplex(i,j)=BiasVal(itemp)
          KGEValComplex(i,j)=KGEVal(itemp)
          LogNSEValComplex(i,j)=LogNSEVal(itemp)
          RMSEValComplex(i,j)=RMSEVal(itemp)
          PeaValComplex(i,j)=PeaVal(itemp)
          SpeValComplex(i,j)=SpeVal(itemp)
    enddo
enddo


!Duan step (4) evolve each partitions p ussing CCE algorithm
!start of loop over p partitions
do y=1,NoP
  do i=1,NoM
       DeepSoilCond1(i)=DeepSoilCondComplex(i,y)  
       ShallowSoilCond1(i)=ShallowSoilCondComplex(i,y)  
       ShallowSoilDepth1(i)=ShallowSoilDepthComplex(i,y)  
       AePeRatio1(i)=AePeRatioComplex(i,y)  
!       Strickler1(i)=StricklerComplex(i,y)  
       Urban1(i)=UrbanComplex(i,y)  
       NSECal(i)=NSECalComplex(i,y)  
       BiasCal(i)=BiasCalComplex(i,y)  
       KGECal(i)=KGECalComplex(i,y)  
       LogNSECal(i)=logNSECalComplex(i,y)  
       RMSECal(i)=RMSECalComplex(i,y)  
       PeaCal(i)=PeaCalComplex(i,y)  
       SpeCal(i)=SpeCalComplex(i,y)  
       NSEVal(i)=NSEValComplex(i,y)  
       BiasVal(i)=BiasValComplex(i,y)  
       KGEVal(i)=KGEValComplex(i,y)  
       LogNSEVal(i)=logNSEValComplex(i,y)  
       RMSEVal(i)=RMSEValComplex(i,y)  
       SpeVal(i)=SpeValComplex(i,y)  
       SpeVal(i)=SpeValComplex(i,y)  
  enddo

!start of loop over beta values
  do j=1,NoBeta
    simnumber=simnumber+1
    outside1=.false.
    outside2=.false.
!    do i=1,NoM
!     write(643,'(7(f10.4,A1),f10.4)') DeepSoilCond1(i),',',ShallowSoilCond1(i),',',ShallowSoilDepth1(i),',',AePeRatio1(i),',',Strickler1(i),',',Urban1(i),',',NSE(i),',',Bias(i)
!    enddo
!      write(643,*)

!   NSEtempArray=NSE
!order by NSE
Call objectivefunction(NoS,NSECal,BiasCal,KGECal,logNSECal,RMSECal,PeaCal,SpeCal,NSEFrac,BiasFrac,KGEFrac,LogNSEFrac,RMSEFrac,PeaFrac,SpeFrac,Objfn)
!    ObjFn=NSECal
Call orderdata(NoS,ObjFn,DeepSoilCond1,ShallowSoilCond1,ShallowSoilDepth1,AePeRatio1,Urban1,NSECal,BiasCal,KGEcal,LogNSECal,RMSECal,PeaCal,SpeCal,NSEVal,BiasVal,KGEVal,LogNSEVal,RMSEVal,PeaVal,Speval)

  do i=1,NoM
    write(19,'(3(I4,A1),18(f10.4,A1),f10.4)') x,',',y,',',j,',',DeepSoilCond1(i),',',ShallowSoilCond1(i),',',ShallowSoilDepth1(i),',',AePeRatio1(i),',',Urban1(i),',',NSECal(i),',',BiasCal(i),',',KGECal(i),',',LogNSECal(i),',',RMSECal(i),',',PeaCal(i),',',SpeCal(i),',',NSEVal(i),',',BiasVal(i),',',KGEVal(i),',',LogNSEVal(i),',',RMSEVal(i),',',PeaVal(i),',',SpeVal(i)
  enddo
      write(19,*)
!   do i=1,NoM
!      write(644,'(7(f10.4,A1),f10.4)') DeepSoilCond1(i),',',ShallowSoilCond1(i),',',ShallowSoilDepth1(i),',',AePeRatio1(i),',',Strickler1(i),',',Urban1(i),',',NSE(i),',',Bias(i)
!   enddo
!      write(644,*)

!CCE Step 1 (Duan-jh-1994)
    do i=1,NoM
    rtemp1=2*(NoM+1-i)
    rtemp2=NoM*(NoM+1)
    prob(i)=rtemp1/rtemp2
    call random_number(random)
!    print*,random
    ProbTimesRandom(i)=prob(i)*random
   enddo

!order by ProbTimesRandom
  Call orderdata(NoM,ProbTimesRandom,DeepSoilCond1,ShallowSoilCond1,ShallowSoilDepth1,AePeRatio1,Urban1,NSECal,BiasCal,KGEcal,LogNSECal,RMSECal,PeaCal,SpeCal,NSEVal,BiasVal,KGEVal,LogNSEVal,RMSEVal,PeaVal,Speval)

!  do i=1,NoM
!    write(645,'(7(f10.4,A1),f10.4)') DeepSoilCond1(i),',',ShallowSoilCond1(i),',',ShallowSoilDepth1(i),',',AePeRatio1(i),',',Strickler1(i),',',Urban1(i),',',NSE(i),',',Bias(i)
!  enddo
!    write(645,*)


!  do i=1,NoS
!    print*,NSE(i),DeepSoilCond1(i),ProbTimesRandom(i)
!  enddo

!CCE Step ii and iii (Duan-jh-1994)
!order q point using NSE, q is now worst point
!  NSEtempArray=NSE
  Call objectivefunction(NoQ,NSECal,BiasCal,KGECal,logNSECal,RMSECal,PeaCal,SpeCal,NSEFrac,BiasFrac,KGEFrac,LogNSEFrac,RMSEFrac,PeaFrac,SpeFrac,Objfn)
!    ObjFn=NSECal
  Call orderdata(NoQ,ObjFn(1:NoQ),DeepSoilCond1(1:NoQ),ShallowSoilCond1(1:NoQ),ShallowSoilDepth1(1:NoQ),AePeRatio1(1:NoQ),Urban1(1:NoQ),NSECal(1:NoQ),BiasCal(1:NoQ),KGEcal(1:NoQ),LogNSECal(1:NoQ),RMSECal(1:NoQ),PeaCal(1:NoQ),SpeCal(1:NoQ),NSEVal(1:NoQ),BiasVal(1:NoQ),KGEVal(1:NoQ),LogNSEVal(1:NoQ),RMSEVal(1:NoQ),PeaVal(1:NoQ),Speval(1:NoQ))

!  do i=1,NoQ
!     write(646,'(7(f10.4,A1),f10.4)') DeepSoilCond1(i),',',ShallowSoilCond1(i),',',ShallowSoilDepth1(i),',',AePeRatio1(i),',',Strickler1(i),',',Urban1(i),',',NSE(i),',',Bias(i)
!   enddo
!      write(646,*)


! put deep soil conductivity over a log scale
! reflection step
!    temp=sum(log10(DeepSoilCond1(1:NoQ-1)))/real(NoQ-1)
!    DeepSoilCondTemp=10.0**temp+(10.0**temp-DeepSoilCond1(NoQ))

! reflection 
    temp=sum(DeepSoilCond1(1:NoQ-1))/real(NoQ-1)
    DeepSoilCondTemp1=temp+(temp-DeepSoilCond1(NoQ))
    temp=sum(ShallowSoilCond1(1:NoQ-1))/real(NoQ-1)
    ShallowSoilCondTemp1=temp+(temp-ShallowSoilCond1(NoQ))
    temp=sum(ShallowSoilDepth1(1:NoQ-1))/real(NoQ-1)
    ShallowSoilDepthTemp1=temp+(temp-ShallowSoilDepth1(NoQ))
    temp=sum(AePeRatio1(1:NoQ-1))/real(NoQ-1)
    AePeRatioTemp1=temp+(temp-AePeRatio1(NoQ))
!    temp=sum(Strickler1(1:NoQ-1))/real(NoQ-1)
!    StricklerTemp1=temp+(temp-Strickler1(NoQ))
    temp=sum(Urban1(1:NoQ-1))/real(NoQ-1)
    UrbanTemp1=temp+(temp-Urban1(NoQ))

!if ((DeepSoilCondTemp.lt.DeepSoilCond(1)).or.(DeepSoilCondTemp.gt.DeepSoilCond(2))) outside(1)=.true.
   if ((DeepSoilCondTemp1.lt.0.0).or.(DeepSoilCondTemp1.gt.DeepSoilCond(2)*2.0)) outside1(1)=.true.
!if ((ShallowSoilCondTemp.lt.ShallowSoilCond(1)).or.(ShallowSoilCondTemp.gt.ShallowSoilCond(2))) outside(2)=.true.
  if ((ShallowSoilCondTemp1.lt.0.0).or.(ShallowSoilCondTemp1.gt.ShallowSoilCond(2)*2.0)) outside1(2)=.true.
!if ((ShallowSoilDepthTemp.lt.ShallowSoilDepth(1)).or.(ShallowSoilDepthTemp.gt.ShallowSoilDepth(2))) outside(3)=.true.
  if ((ShallowSoilDepthTemp1.lt.0.1).or.(ShallowSoilDepthTemp1.gt.ShallowSoilDepth(2)*2.0)) outside1(3)=.true.
!if ((AePeRatioTemp.lt.AePeRatio(1)).or.(AePeRatioTemp.gt.AePeRatio(2))) outside(4)=.true.
  if ((AePeRatioTemp1.lt.0.1).or.(AePeRatioTemp1.gt.AePeRatio(2)*2.0)) outside1(4)=.true.
!if ((StricklerTemp.lt.Strickler(1)).or.(StricklerTemp.gt.Strickler(2))) outside(5)=.true.
!  if ((StricklerTemp1.lt.0.1).or.(StricklerTemp1.gt.Strickler(2)*2.0)) outside1(5)=.true.
!if ((UrbanTemp.lt.Urban(1)).or.(UrbanTemp.gt.Urban(2))) outside(6)=.true.
  if ((UrbanTemp1.lt.0.0).or.(UrbanTemp1.gt.Urban(2)*2.0)) outside1(5)=.true.

! half reflection
    temp=sum(DeepSoilCond1(1:NoQ-1))/real(NoQ-1)
    DeepSoilCondTemp2=temp+0.5*(temp-DeepSoilCond1(NoQ))
    temp=sum(ShallowSoilCond1(1:NoQ-1))/real(NoQ-1)
    ShallowSoilCondTemp2=temp+0.5*(temp-ShallowSoilCond1(NoQ))
    temp=sum(ShallowSoilDepth1(1:NoQ-1))/real(NoQ-1)
    ShallowSoilDepthTemp2=temp+0.5*(temp-ShallowSoilDepth1(NoQ))
    temp=sum(AePeRatio1(1:NoQ-1))/real(NoQ-1)
    AePeRatioTemp2=temp+0.5*(temp-AePeRatio1(NoQ))
!    temp=sum(Strickler1(1:NoQ-1))/real(NoQ-1)
!    StricklerTemp2=temp+0.5*(temp-Strickler1(NoQ))
    temp=sum(Urban1(1:NoQ-1))/real(NoQ-1)
    UrbanTemp2=temp+0.5*(temp-Urban1(NoQ))

!if ((DeepSoilCondTemp.lt.DeepSoilCond(1)).or.(DeepSoilCondTemp.gt.DeepSoilCond(2))) outside(1)=.true.
   if ((DeepSoilCondTemp2.lt.0.0).or.(DeepSoilCondTemp2.gt.DeepSoilCond(2)*2.0)) outside2(1)=.true.
!if ((ShallowSoilCondTemp.lt.ShallowSoilCond(1)).or.(ShallowSoilCondTemp.gt.ShallowSoilCond(2))) outside(2)=.true.
  if ((ShallowSoilCondTemp2.lt.0.0).or.(ShallowSoilCondTemp2.gt.ShallowSoilCond(2)*2.0)) outside2(2)=.true.
!if ((ShallowSoilDepthTemp.lt.ShallowSoilDepth(1)).or.(ShallowSoilDepthTemp.gt.ShallowSoilDepth(2))) outside(3)=.true.
  if ((ShallowSoilDepthTemp2.lt.0.1).or.(ShallowSoilDepthTemp2.gt.ShallowSoilDepth(2)*2.0)) outside2(3)=.true.
!if ((AePeRatioTemp.lt.AePeRatio(1)).or.(AePeRatioTemp.gt.AePeRatio(2))) outside(4)=.true.
  if ((AePeRatioTemp2.lt.0.1).or.(AePeRatioTemp2.gt.AePeRatio(2)*2.0)) outside2(4)=.true.
!if ((StricklerTemp.lt.Strickler(1)).or.(StricklerTemp.gt.Strickler(2))) outside(5)=.true.
!  if ((StricklerTemp2.lt.0.1).or.(StricklerTemp2.gt.Strickler(2)*2.0)) outside2(5)=.true.
!if ((UrbanTemp.lt.Urban(1)).or.(UrbanTemp.gt.Urban(2))) outside(6)=.true.
  if ((UrbanTemp2.lt.0.0).or.(UrbanTemp2.gt.Urban(2)*2.0)) outside2(5)=.true.

  
    
  !if outside feasible space then random point
  If ((Any(outside1)).and.(Any(outside2))) then

      
      !Step VI random point
    call random_number(random)
!    print*,random
!    temp=(log10(DeepSoilCond(2))-log10(DeepSoilCond(1)))*random+log10(DeepSoilCond(1))
!    DeepSoilCond1(NoQ)=10.0**temp
    DeepSoilCond1(NoQ)=(DeepSoilCond(2)-DeepSoilCond(1))*random+DeepSoilCond(1)
    call random_number(random)
    ShallowSoilCond1(NoQ)=(ShallowSoilCond(2)-ShallowSoilCond(1))*random+ShallowSoilCond(1)
    call random_number(random)
    ShallowSoilDepth1(NoQ)=(ShallowSoilDepth(2)-ShallowSoilDepth(1))*random+ShallowSoilDepth(1)
    call random_number(random)
    AePeRatio1(NoQ)=(AePeRatio(2)-AePeRatio(1))*random+AePeRatio(1)
!    call random_number(random)
!    Strickler1(NoQ)=(Strickler(2)-Strickler(1))*random+Strickler(1)
    call random_number(random)
    Urban1(NoQ)=(Urban(2)-Urban(1))*random+Urban(1)
    Call WriteXmlFile(dischargefilename,rundatafilename,libraryfile,filelength,text,simnumber,vegend,soilpropend,soildetend,NoDischargeValues,calibrationTimes,validationTimes,NoS,MeasuredDischarge, &
      AePeRatio1(NoQ),ShallowSoilCond1(NoQ),DeepSoilCond1(NoQ),ShallowSoilDepth1(NoQ),Urban1(NoQ),BiasCal(NoQ),NSECal(NoQ),KGECal(NoQ),LogNSECal(NoQ),RMSECal(NoQ),PeaCal(NoQ),SpeCal(NoQ),BiasVal(NoQ),NSEval(NoQ),KGEVal(NoQ),LogNSEVal(NoQ),RMSEVal(NoQ),PeaVal(NoQ),Speval(NoQ),nrows,ncols,Vegtypenew, &
      calstart,calfin,valstart,valfinish,GroundLevelFilePath,SimulatedLevelFilepath)

  else
!inside feasible space
     if (Any(outside1)) then
!half reflection
        Call WriteXmlFile(dischargefilename,rundatafilename,libraryfile,filelength,text,simnumber,vegend,soilpropend,soildetend,NoDischargeValues,calibrationTimes,validationTimes,NoS,MeasuredDischarge, &
          AePeRatioTemp2,ShallowSoilCondTemp2,DeepSoilCondTemp2,ShallowSoilDepthTemp2,UrbanTemp2,BiasCalTemp,NSECalTemp,KGECalTemp,LogNSECalTemp,RMSECalTemp,PeaCalTemp,SpeCalTemp,BiasValTemp,NSEvalTemp,KGEValTemp,LogNSEValTemp,RMSEValTemp,PeaValTemp,SpevalTemp,nrows,ncols,Vegtypenew, &
          calstart,calfin,valstart,valfinish,GroundLevelFilePath,SimulatedLevelFilepath)
     else
!full reflection
         Call WriteXmlFile(dischargefilename,rundatafilename,libraryfile,filelength,text,simnumber,vegend,soilpropend,soildetend,NoDischargeValues,calibrationTimes,validationTimes,NoS,MeasuredDischarge, &
         AePeRatioTemp1,ShallowSoilCondTemp1,DeepSoilCondTemp1,ShallowSoilDepthTemp1,UrbanTemp1,BiasCalTemp,NSECalTemp,KGECalTemp,LogNSECalTemp,RMSECalTemp,PeaCalTemp,SpeCalTemp,BiasValTemp,NSEvalTemp,KGEValTemp,LogNSEValTemp,RMSEValTemp,PeaValTemp,SpevalTemp,nrows,ncols,Vegtypenew, &
         calstart,calfin,valstart,valfinish,GroundLevelFilePath,SimulatedLevelFilepath)
     endif
         

!    print*, DeepSoilCondTemp,ShallowSoilCondTemp,ShallowSoilDepthTemp,AePeRatioTemp,StricklerTemp, UrbanTemp

  endif

Call objectivefunction1(NSECalTemp,BiasCalTemp,KGECalTemp,logNSECalTemp,RMSECalTemp,PeaCalTemp,SpeCalTemp,NSEFrac,BiasFrac,KGEFrac,LogNSEFrac,RMSEFrac,PeaFrac,SpeFrac,ObjfnTemp)

!  ObjFnTemp=NSECalTemp 
!CCE Step iV and V (Duan-jh-1994)
  if ((Any(outside1)).and.(Any(outside2))) then
    write(678,*) 'outside feasible space'
  else if ((ObjFnTemp.GT.ObjFn(NoQ))) then
!reflection is good
    if (Any(outside1)) then
       write(678,*) 'reflection-half'
       DeepSoilCond1(NoQ)=DeepSoilCondTemp2
       ShallowSoilCond1(NoQ)=ShallowSoilCondTemp2
       ShallowSoilDepth1(NoQ)= ShallowSoilDepthTemp2
       AePeRatio1(NoQ)=AePeRatioTemp2
!       Strickler1(NoQ)=StricklerTemp2
       Urban1(NoQ)=UrbanTemp2
       BiasCal(NoQ)=BiasCalTemp
       NSECal(NoQ)=NSECalTemp
       KGECal(NoQ)=KGECalTemp
       logNSECal(NoQ)=LogNSECalTemp
       RMSECal(NoQ)=RMSECalTemp
       PeaCal(NoQ)=PeaCalTemp
       SpeCal(NoQ)=SpeCalTemp
       BiasVal(NoQ)=BiasValTemp
       NSEVal(NoQ)=NSEValTemp
       KGEVal(NoQ)=KGEValTemp
       logNSEVal(NoQ)=LogNSEValTemp
       RMSEVal(NoQ)=RMSEValTemp
       PeaVal(NoQ)=PeaValTemp
       SpeVal(NoQ)=SpeValTemp
    else
       write(678,*) 'reflection-full'
      DeepSoilCond1(NoQ)=DeepSoilCondTemp1
      ShallowSoilCond1(NoQ)=ShallowSoilCondTemp1
      ShallowSoilDepth1(NoQ)= ShallowSoilDepthTemp1
      AePeRatio1(NoQ)=AePeRatioTemp1
!      Strickler1(NoQ)=StricklerTemp1
      Urban1(NoQ)=UrbanTemp1
       BiasCal(NoQ)=BiasCalTemp
       NSECal(NoQ)=NSECalTemp
       KGECal(NoQ)=KGECalTemp
       logNSECal(NoQ)=LogNSECalTemp
       RMSECal(NoQ)=RMSECalTemp
       PeaCal(NoQ)=PeaCalTemp
       SpeCal(NoQ)=SpeCalTemp
       BiasVal(NoQ)=BiasValTemp
       NSEVal(NoQ)=NSEValTemp
       KGEVal(NoQ)=KGEValTemp
       logNSEVal(NoQ)=LogNSEValTemp
       RMSEVal(NoQ)=RMSEValTemp
       PeaVal(NoQ)=PeaValTemp
       SpeVal(NoQ)=SpeValTemp
    endif
  else
! contraction    
!    temp=sum(log10(DeepSoilCond1(1:NoQ-1)))/real(NoQ-1)
!    DeepSoilCondTemp=(10.0**temp+DeepSoilCond1(NoQ))/2.0
    temp=sum(DeepSoilCond1(1:NoQ-1))/real(NoQ-1)
    DeepSoilCondTemp=(temp+DeepSoilCond1(NoQ))/2.0
    temp=sum(ShallowSoilCond1(1:NoQ-1))/real(NoQ-1)
    ShallowSoilCondTemp=(temp+ShallowSoilCond1(NoQ))/2.0
    temp=sum(ShallowSoilDepth1(1:NoQ-1))/real(NoQ-1)
    ShallowSoilDepthTemp=(temp+ShallowSoilDepth1(NoQ))/2.0
    temp=sum(AePeRatio1(1:NoQ-1))/real(NoQ-1)
    AePeRatioTemp=(temp+AePeRatio1(NoQ))/2.0
!    temp=sum(Strickler1(1:NoQ-1))/real(NoQ-1)
!    StricklerTemp=(temp+Strickler1(NoQ))/2.0
    temp=sum(Urban1(1:NoQ-1))/real(NoQ-1)
    UrbanTemp=(temp+Urban1(NoQ))/2.0

    Call WriteXmlFile(dischargefilename,rundatafilename,libraryfile,filelength,text,simnumber,vegend,soilpropend,soildetend,NoDischargeValues,calibrationTimes,validationTimes,NoS,MeasuredDischarge, &
      AePeRatioTemp,ShallowSoilCondTemp,DeepSoilCondTemp,ShallowSoilDepthTemp,UrbanTemp,BiasCalTemp,NSECalTemp,KGECalTemp,LogNSECalTemp,RMSECalTemp,PeaCalTemp,SpeCalTemp,BiasValTemp,NSEvalTemp,KGEValTemp,LogNSEValTemp,RMSEValTemp,PeaValTemp,SpevalTemp,nrows,ncols,Vegtypenew, &
      calstart,calfin,valstart,valfinish,GroundLevelFilePath,SimulatedLevelFilepath)

!      write(651,'(7(f10.4,A1),f10.4)') DeepSoilCondTemp,',',ShallowSoilCondTemp,',',ShallowSoilCondTemp,',',AePeRatioTemp,',',StricklerTemp,',',UrbanTemp,',',NSETemp,',',BiasTemp
!       write(651,*) NoQ,NSE(NoQ)
 
    
Call objectivefunction1(NSECalTemp,BiasCalTemp,KGECalTemp,logNSECalTemp,RMSECalTemp,PeaCalTemp,SpeCalTemp,NSEFrac,BiasFrac,KGEFrac,LogNSEFrac,RMSEFrac,PeaFrac,SpeFrac,ObjfnTemp)
    
!     ObjFnTemp=NSECalTemp 
    
    if (ObjFnTemp.GT.ObjFn(NoQ)) then
! contraction is good
        write(678,*) 'contraction'
       DeepSoilCond1(NoQ)=DeepSoilCondTemp
       ShallowSoilCond1(NoQ)=ShallowSoilCondTemp
       ShallowSoilDepth1(NoQ)= ShallowSoilDepthTemp
       AePeRatio1(NoQ)=AePeRatioTemp
!       Strickler1(NoQ)=StricklerTemp
       Urban1(NoQ)=UrbanTemp
       BiasCal(NoQ)=BiasCalTemp
       NSECal(NoQ)=NSECalTemp
       KGECal(NoQ)=KGECalTemp
       logNSECal(NoQ)=LogNSECalTemp
       RMSECal(NoQ)=RMSECalTemp
       PeaCal(NoQ)=PeaCalTemp
       SpeCal(NoQ)=SpeCalTemp
       BiasVal(NoQ)=BiasValTemp
       NSEVal(NoQ)=NSEValTemp
       KGEVal(NoQ)=KGEValTemp
       logNSEVal(NoQ)=LogNSEValTemp
       RMSEVal(NoQ)=RMSEValTemp
       PeaVal(NoQ)=PeaValTemp
       SpeVal(NoQ)=SpeValTemp
      
    else
!Step VI random point
    call random_number(random)
    write(678,*) 'random'
!    temp=(log10(DeepSoilCond(2))-log10(DeepSoilCond(1)))*random+log10(DeepSoilCond(1))
!    DeepSoilCond1(NoQ)=10.0**temp
    DeepSoilCond1(NoQ)=(DeepSoilCond(2)-DeepSoilCond(1))*random+DeepSoilCond(1)
    call random_number(random)
    ShallowSoilCond1(NoQ)=(ShallowSoilCond(2)-ShallowSoilCond(1))*random+ShallowSoilCond(1)
    call random_number(random)
    ShallowSoilDepth1(NoQ)=(ShallowSoilDepth(2)-ShallowSoilDepth(1))*random+ShallowSoilDepth(1)
    call random_number(random)
    AePeRatio1(NoQ)=(AePeRatio(2)-AePeRatio(1))*random+AePeRatio(1)
!    call random_number(random)
!    Strickler1(NoQ)=(Strickler(2)-Strickler(1))*random+Strickler(1)
    call random_number(random)
    Urban1(NoQ)=(Urban(2)-Urban(1))*random+Urban(1)

!write XML file, run shetran and work out NSE
    Call WriteXmlFile(dischargefilename,rundatafilename,libraryfile,filelength,text,simnumber,vegend,soilpropend,soildetend,NoDischargeValues,calibrationTimes,validationTimes,NoS,MeasuredDischarge, &
      AePeRatio1(NoQ),ShallowSoilCond1(NoQ),DeepSoilCond1(NoQ),ShallowSoilDepth1(NoQ),Urban1(NoQ),BiasCal(NoQ),NSECal(NoQ),KGECal(NoQ),LogNSECal(NoQ),RMSECal(NoQ),PeaCal(NoQ),SpeCal(NoQ),BiasVal(NoQ),NSEval(NoQ),KGEVal(NoQ),LogNSEVal(NoQ),RMSEVal(NoQ),PeaVal(NoQ),Speval(NoQ),nrows,ncols,Vegtypenew, &
      calstart,calfin,valstart,valfinish,GroundLevelFilePath,SimulatedLevelFilepath)
        
    endif

  endif
!  write(15,*) 'complex',  y  ,'beta in complex',j
  write(15,'(i4,A1,18(f10.4,A1),f10.4)') Simnumber,',',DeepSoilCond1(NoQ),',',ShallowSoilCond1(NoQ),',',ShallowSoilDepth1(NoQ),',',AePeRatio1(NoQ),',',Urban1(NoQ),',',NSECal(NoQ),',',BiasCal(NoQ),',',KGECal(NoQ),',',LogNSECal(NoQ),',',RMSECal(noQ),',',PeaCal(noQ),',',SpeCal(NoQ),',',NSEVal(NoQ),',',BiasVal(NoQ),',',KGEVal(NoQ),',',LogNSEVal(NoQ),',',RMSEVal(NoQ),',',PeaVal(noQ),',',SpeVal(NoQ)
  if (ObjFn(NoQ).GT.ObjFnBest) then
      Simbest = simnumber
      ObjFnBest = ObjFn(NoQ)
  endif

!  do i=1,NoM
!    write(19,'(2I4,7(f10.4,A1),f10.4)') y,j,DeepSoilCond1(i),',',ShallowSoilCond1(i),',',ShallowSoilDepth1(i),',',AePeRatio1(i),',',Strickler1(i),',',Urban1(i),',',NSE(i),',',Bias(i)
!  enddo
!      write(19,*)

!end of loop over beta values
  enddo 

  
    do i=1,NoM
       DeepSoilCondComplex(i,y)=DeepSoilCond1(i)
       ShallowSoilCondComplex(i,y)=ShallowSoilCond1(i)  
       ShallowSoilDepthComplex(i,y)=ShallowSoilDepth1(i) 
       AePeRatioComplex(i,y)=AePeRatio1(i)  
!       StricklerComplex(i,y)= Strickler1(i) 
       UrbanComplex(i,y)=Urban1(i) 
       NSECalComplex(i,y)=NSECal(i)
       BiasCalComplex(i,y)=BiasCal(i)
       KGECalComplex(i,y)=KGECal(i)
       LogNSECalComplex(i,y)=LogNSECal(i)
       RMSECalComplex(i,y)=RMSECal(i)
       PeaCalComplex(i,y)=PeaCal(i)
       SpeCalComplex(i,y)=SpeCal(i)
       NSEValComplex(i,y)=NSEVal(i)
       BiasValComplex(i,y)=BiasVal(i)
       KGEValComplex(i,y)=KGEVal(i)
       LogNSEValComplex(i,y)=LogNSEVal(i)
       RMSEValComplex(i,y)=RMSEVal(i)
       PeaValComplex(i,y)=peaVal(i)
       SpeValComplex(i,y)=SpeVal(i)
    enddo


!end loop over partitions
  enddo

!do i=1,NoS
!    write(542,*) NSE(i)
!enddo
!    write(542,*) 
!*****************************
!Duan step (5) combine complexes 
  do i=1,NoM
    do j=1,NoP
          itemp=NoP*(i-1) + j
          DeepSoilCond1(itemp)=DeepSoilCondComplex(i,j)
          ShallowSoilCond1(itemp)=ShallowSoilCondComplex(i,j)
          ShallowSoilDepth1(itemp)=ShallowSoilDepthComplex(i,j)
          AePeRatio1(itemp)=AePeRatioComplex(i,j)
  !        Strickler1(itemp)=StricklerComplex(i,j)
          Urban1(itemp)=UrbanComplex(i,j)
          NSECal(itemp)=NSECalComplex(i,j)
          BiasCal(itemp)=BiasCalComplex(i,j)
          KGECal(itemp)=KGECalComplex(i,j)
          LogNSECal(itemp)=LogNSECalComplex(i,j)
          RMSECal(itemp)=RMSECalComplex(i,j)
          peaCal(itemp)=PeaCalComplex(i,j)
          SpeCal(itemp)=SpeCalComplex(i,j)
          NSEVal(itemp)=NSEValComplex(i,j)
          BiasVal(itemp)=BiasValComplex(i,j)
          KGEVal(itemp)=KGEValComplex(i,j)
          LogNSEVal(itemp)=LogNSEValComplex(i,j)
          RMSEVal(itemp)=RMSEValComplex(i,j)
          peaVal(itemp)=PeaValComplex(i,j)
          SpeVal(itemp)=SpeValComplex(i,j)
    enddo
  enddo

!*****************************
!Duan go  back to step 2 (new iteration) where the sample is sorted and split into complexes. In Duan this is done in Step 5.
enddo


print*, 'run best simulation again'
!Find the best simulation and run it
!NSEtempArray=NSE
Call orderdata(NoS,ObjFn,DeepSoilCond1,ShallowSoilCond1,ShallowSoilDepth1,AePeRatio1,Urban1,NSECal,BiasCal,KGEcal,LogNSECal,RMSECal,PeaCal,SpeCal,NSEVal,BiasVal,KGEVal,LogNSEVal,RMSEVal,PeaVal,Speval)

!use simbest
Call WriteXmlFile(dischargefilename,rundatafilename,libraryfile,filelength,text,simbest,vegend,soilpropend,soildetend,NoDischargeValues,calibrationTimes,validationTimes,NoS,MeasuredDischarge, &
      AePeRatio1(1),ShallowSoilCond1(1),DeepSoilCond1(1),ShallowSoilDepth1(1),Urban1(1),BiasCal(1),NSECal(1),KGECal(1),LogNSECal(1),RMSECal(1),PeaCal(1),SpeCal(1),BiasVal(1),NSEVal(1),KGEVal(1),LogNSEVal(1),RMSEVal(1),PeaVal(1),Speval(1),nrows,ncols,Vegtypenew, &
     calstart,calfin,valstart,valfinish,GroundLevelFilePath,SimulatedLevelFilepath)



write(15,*) 
write(15,'(A18,I4)') 'Best Simulation = ',simbest
write(15,'(A216)') 'number,deep_soil_conductivity,shallow_soil_conductivity,shallow_soil_depth,AePe_ratio,urban-prec-fraction,NSECal,BiasCal,KGECal,LogNSECal,RMSECal,PeaCal,SpeCal,NSEVal,BiasVal,KGEVal,LogNSEVal,RMSEVal,PeaVal,Speval'
write(15,'(i4,A1,20(f10.4,A1),f10.4)') Simnumber,',',DeepSoilCond1(1),',',ShallowSoilCond1(1),',',ShallowSoilDepth1(1),',',AePeRatio1(1),',',Urban1(1),',',NSECal(1),',',BiasCal(1),',',KGECal(1),',',LogNSECal(1),',',RMSECal(1),',',PeaCal(1),',',SpeCal(1),',',NSEVal(1),',',BiasVal(1),',',KGEVal(1),',',LogNSEVal(1),',',RMSEVal(1),',',PeaVal(1),',',SpeVal(1)
!correct reduced precipitation in urban areas. Seperate sewers  

!pause
 stop

9999 write (*,*) 'Error openinig file ',optimisefull
write(*,'(''paused, type [enter] to continue'')')
read (*,*)
stop 
9998 write (*,*) 'Error openinig file ',libraryfilepath
write(*,'(''paused, type [enter] to continue'')')
read (*,*)
stop 
9997 write (*,*) 'Error openinig file ',dischargefilepath
write(*,'(''paused, type [enter] to continue'')')
read (*,*)
stop 
9995 write (*,*) 'Error openinig file ',vegfilename
write(*,'(''paused, type [enter] to continue'')')
read (*,*)
stop 
9994 write (*,*) 'Error openinig file ',soilfilename
write(*,'(''paused, type [enter] to continue'')')
read (*,*)
stop 



end program ShetranAutoCalibration


subroutine WriteXmlFile(dischargefilename,rundatafilename,libraryfile,filelength,text,simnumber,vegend,soilpropend,soildetend,NoDischargeValues,calibrationTimes,validationTimes,NoS,MeasuredDischarge, &
      SubAePeRatio,SubShallowSoilCond,SubDeepSoilCond,SubShallowSoilDepth,SubUrban,subBiasCal,subNSECal,subKGECal,subLogNSECal,RMSECal,pearsonCal,spearmanCal,subBiasVal,subNSEVal,subKGEval,subLogNSEVal, &
      RMSEVal,pearsonVal,spearmanVal,nrows,ncols,Vegtypenew,calstart,calfin,valstart,valfinish,GroundLevelFilePath,SimulatedLevelFilepath)

implicit none


character(300),intent(in) :: dischargefilename
character(300),intent(in) :: rundatafilename,libraryfile
integer,intent(in)        :: filelength
character(300),intent(in) :: text(filelength)
integer,intent(in)        :: simnumber,vegend,soilpropend,soildetend,NoDischargeValues,calibrationTimes(2),validationTimes(2),NoS  
real,intent(in)           :: MeasuredDischarge(NoDischargeValues)
real,intent(in)           :: SubAePeRatio,SubShallowSoilCond,SubDeepSoilCond,SubShallowSoilDepth,SubUrban
real,intent(out)          :: subBiasCal,subNSECal,subKGECal,subLogNSECal,subBiasVal,subNSEVal,subKGEval,subLogNSEVal
real,intent(out)          :: RMSECal,RMSEVal,pearsonCal,spearmanCal,pearsonVal,spearmanVal
integer,intent(in)        :: nrows,ncols,Vegtypenew(nrows,ncols)
integer,intent(in)        :: calstart,calfin,valstart,valfinish
character(300),intent(in) :: GroundLevelFilePath,SimulatedLevelFilepath


!integer        :: nvalues
!integer        :: Measuredday,Measuredmonth,Measuredyear,StartDate,CurrentDate(1000)
!integer        :: julian_date

!character(300)  :: measuredWaterTable,simulatedWaterTable 

!local 
character(300)            :: simnumberstring,libraryfileoutput,header
integer                   :: z,nosimdischargevalues,mindischargevalues,io
real                      :: MeanMeasured,MeanSimulated
real                      :: MeanLogMeasured,MeanLogSimulated
real                      :: SimulatedDischarge(NoDischargeValues),logMeasured(NoDischargeValues),logSimulated(NoDischargeValues)
real                      :: MeasuredSquare(NoDischargeValues),SimuMeasSquare(NoDischargeValues)
real                      :: LogMeasuredSquare(NoDischargeValues),LogSimuMeasSquare(NoDischargeValues)
real                      :: SimulatedSquare(NoDischargeValues),SimMeas(NoDischargeValues)
real                      :: stdDevMeas,stdDevSim,SumMeasuredSquare,SumSimulatedSquare,sumSimMeas,correlation
real                      :: precip(NoDischargeValues)
real                      :: measuredPSL(1000),SimulatedPSLMeasuredTimes(1000)
real                      :: SumPslDiffSqur,PslDiffsqure(1000)
real                      :: OrderMeasuredPSL(1000),OrderSimulatedPSLMeasuredTimes(1000),locOrderSimulatedPSLMeasuredTimes(1000)
integer        :: Measuredday,Measuredmonth,Measuredyear,StartDate,CurrentDate(1000)
integer        :: NoValues,julian_date



real                      :: SimulatedPSL(20000)
character(300) :: text1,precipFile,precipFilename,precipfilenameUrban,text2
integer :: namelength,vegorder,i,j,nvalues
real :: vegurban(nrows*ncols)
real:: temp1,temp2
real :: MeanPslMeas,meanPSLSim,CovPSL,StDevMeas,StDevSim,sumRankDiff



    Write( simnumberstring, '(i4)' )simnumber
    simnumberstring=adjustl(simnumberstring)
    libraryfileoutput='../'//trim(libraryfile)//'/LibraryFile'//trim(simnumberstring)//'.xml'
    open(13,FILE=libraryfileoutput)
    do z=1,6
    write (13,'(A)') text(z)
    enddo
    write(13,'(A37)') '<VegMap>landcover-2types.asc</VegMap>'
    write(13,'(A34)') '<SoilMap>soil-2types.asc</SoilMap>'
    do z=9,13
    write (13,'(A)') text(z)
    enddo
    write(13,'(A39,F5.2,A23)') '<VegetationDetail>1, veg1, 3.0, 1, 1.0,',SubAePeRatio,',2.0</VegetationDetail>'
    write(13,'(A68)') '<VegetationDetail>2, urban, 0.3,0.3,0.5,1.0,12.0,</VegetationDetail>'
    do z=vegend,vegend+2
    write (13,'(A)') text(z)
    enddo
    write(13,'(A35,F9.4,A31)') '<SoilProperty>1,Soil1,0.430, 0.010,',SubShallowSoilCond,', 0.0083, 1.2539</SoilProperty>'
    write(13,'(A38,F9.4,A26)') '<SoilProperty>2,Aquifer1,0.200, 0.100,',SubDeepSoilCond,', 0.01, 5.0</SoilProperty>'
    write(13,'(A74)') '<SoilProperty>3,urban1,0.430, 0.010, 0.0001, 0.0083, 1.2539</SoilProperty>'
    do z=soilpropend,soilpropend+2
    write (13,'(A)') text(z)
    enddo
    write(13,'(A20,F6.3,A13)') '<SoilDetail>1, 1, 1,',SubShallowSoilDepth,'</SoilDetail>'
    write(13,'(A39)') '<SoilDetail>1, 2, 2, 50.0</SoilDetail>'
    write(13,'(A39)') '<SoilDetail>2, 1, 3, 0.10</SoilDetail>'
    do z=soildetend,soildetend+1
    write (13,'(A)') text(z)
    enddo

call reducedurbanrain(libraryfile,filelength,text,soildetend,NoDischargeValues,SubUrban,nrows,ncols,Vegtypenew)
    

!Zero PET in urban areas    
    text1=text(soildetend+4)
    namelength=len_trim(text1)
    precipFile=text1(28:namelength-28)
    precipFilename='../'//trim(libraryfile)//'/'//trim(precipFile)
    precipfilenameUrban='../'//trim(libraryfile)//'/'//'Urban_PET.csv'
    open(16,FILE=trim(precipFilename),err=9995,status='old')
    Read(16,*) header
    open(18,FILE=trim(precipfilenameUrban))
    write(18,'(A50)') 'Urban PET. Zero PET in urban areas'
    vegorder=0
     do i=1,nrows
        do j=1,ncols
            if (Vegtypenew(i,j).eq.1) then
               vegorder=vegorder+1
               vegurban(vegorder)=1.0
            elseif (Vegtypenew(i,j).eq.2) then
                !urban
               vegorder=vegorder+1
               vegurban(vegorder)=0.0
            endif
        enddo
    enddo
!!! is the Precip the same timestep as the discharge? *24 if hourly and discharge daily
   do i=1,NoDischargeValues
        Read(16,*) (precip(j),j=1,vegorder)
        write(18,'(*(f6.2, ", "))') (precip(j)*vegurban(j),j=1,vegorder)
    enddo
    close(16)
    close(18)

 
                
    
!    SubUrban
    write(13,'(A87)') '<PrecipitationTimeSeriesData>Urban_Precip.csv</PrecipitationTimeSeriesData> urban precip'
    do z=soildetend+3,soildetend+3
    write (13,'(A)') text(z)
    enddo
    write(13,'(A78)') '<EvaporationTimeSeriesData>Urban_PET.csv</EvaporationTimeSeriesData> urban evap'
    do z=soildetend+5,filelength
    write (13,'(A)') text(z)
    enddo
    close(13)

    CALL execute_command_line('shetran-prepare-snow.exe '//trim(libraryfileoutput))
    CALL execute_command_line('shetran.exe -f '//trim(rundatafilename))
    
!read simulated dischargew
    open(14,FILE=trim(dischargefilename),err=9996,status='old')
    read(14,*)
    nosimdischargevalues=0
    do z=1,NoDischargeValues
      Read(14,*,IOSTAT=io) SimulatedDischarge(z) 
      if (io < 0) then
           EXIT
      else
         nosimdischargevalues=nosimdischargevalues+1
      endif
    enddo
    close(14)
    
!**********************  
!WATER TABLE Caculations
!**********************  
 !water table files
!measuredWaterTable = '../'//trim(libraryfile)//'/GreatThornsFarm.csv'
!simulatedWaterTable =  '../'//trim(libraryfile)//'/output_33029_psl-element133.txt'


!open(20,FILE=measuredWaterTable,err=9989,status='old')
!open(21,FILE=simulatedWaterTable,err=9988,status='old')

open(22,FILE=trim(GroundLevelFilePath),err=9989,status='old')
open(21,FILE=trim(SimulatedLevelFilepath),err=9988,status='old')


    
!read simulated and measured PSL
NoValues=0
do 
  Read(21,'(A11,f10.4)',IOSTAT=io) Text1
  if (io < 0) then
    EXIT
  else
    NoValues= NoValues + 1
  endif
enddo
rewind(21)
Read(21,*) Text1
do i=1,NoValues-1
  Read(21,*) Text1,SimulatedPSL(i) 
enddo
close(21)


StartDate= julian_date(1980,1,1)
NoValues=0
do 
  Read(22,'(A11,f10.4)',IOSTAT=io) Text1
  if (io < 0) then
    EXIT
  else
    NoValues= NoValues + 1
  endif
enddo
rewind(22)
read(22,*)
do i=1,NoValues-1
  Read(22,'(A100)') Text1
!  n2=index(text1,',')
!  endline=len_trim(Text1)
  Text2 = Text1(1:2)
  read(Text2,*) MeasuredDay
  Text2 = Text1(4:5)
  read(Text2,*) MeasuredMonth
  Text2 = Text1(7:10)
  read(Text2,*) MeasuredYear
  CurrentDate(i)=julian_date(MeasuredYear,MeasuredMonth,Measuredday)-StartDate
  Text2 = Text1(12:20)
  read(Text2,*) MeasuredPSL(i)
enddo
close(22)

do i=1,NoValues-1
    SimulatedPSLMeasuredTimes(i)=SimulatedPSL(CurrentDate(i))
enddo






SumPslDiffSqur=0
meanPSLmeas=0
meanPSLSim=0
nvalues=0
covPSL=0
StDevmeas=0
StDevSim=0
OrderMeasuredPSL=MeasuredPSL
OrderSimulatedPSLMeasuredTimes=SimulatedPSLMeasuredTimes
SumRankDiff=0

do z=calstart,calfin
 PslDiffsqure(z) = (SimulatedPSLMeasuredTimes(z)-MeasuredPSL(z))**2
 SumPslDiffSqur = SumPslDiffSqur + PslDiffsqure(z)
 meanPSLmeas=meanPSLmeas + MeasuredPSL(z)
 meanPSLSim=meanPSLSim + SimulatedPSLMeasuredTimes(z)
 nvalues=nvalues+1
enddo
meanPSLmeas=meanPSLmeas/nvalues
meanPSLsim=meanPSLsim/nvalues
RMSECal=sqrt(SumPslDiffSqur/nvalues)


do z=calstart,calfin
    covPSL=covPSL+(MeasuredPSL(z)-meanPSLmeas)*(simulatedPSLMeasuredTimes(z)-meanPSLsim)
    StDevmeas = StDevmeas + (MeasuredPSL(z)-meanPSLmeas)**2
    StDevSim = StDevSim + (simulatedPSLMeasuredTimes(z)-meanPSLsim)**2
enddo
StDevmeas = sqrt(StDevmeas)
StDevSim = sqrt(StDevSim)
PearsonCal = covPSL/(StDevmeas*StDevSim)


!print*, meanPSLmeas, meanPSLsim, RMSECal, covPSL, StDevmeas, StDevSim, PearsonCal

do z=calstart,calfin
 do i=calstart,calfin-1
  if (OrderMeasuredPSL(i).le.OrderMeasuredPSL(i+1)) then
    temp1 = OrderMeasuredPSL(i)
    OrderMeasuredPSL(i) = OrderMeasuredPSL(i+1)
    OrderMeasuredPSL(i+1) = temp1
    temp1 = OrderSimulatedPSLMeasuredTimes(i)
    OrderSimulatedPSLMeasuredTimes(i) = OrderSimulatedPSLMeasuredTimes(i+1)
    OrderSimulatedPSLMeasuredTimes(i+1) = temp1
    endif
 enddo
enddo
do z=calstart,calfin
    locOrderSimulatedPSLMeasuredTimes(z)=z
enddo
do z=calstart,calfin
 do i=calstart,calfin-1
  if (OrderSimulatedPSLMeasuredTimes(i).le.OrderSimulatedPSLMeasuredTimes(i+1)) then
    temp1 = OrderSimulatedPSLMeasuredTimes(i)
    OrderSimulatedPSLMeasuredTimes(i) = OrderSimulatedPSLMeasuredTimes(i+1)
    OrderSimulatedPSLMeasuredTimes(i+1) = temp1
    temp1 = locOrderSimulatedPSLMeasuredTimes(i)
    locOrderSimulatedPSLMeasuredTimes(i) = locOrderSimulatedPSLMeasuredTimes(i+1)
    locOrderSimulatedPSLMeasuredTimes(i+1) = temp1
    endif
 enddo
enddo
do z=calstart,calfin
    SumRankDiff=SumRankDiff + (locOrderSimulatedPSLMeasuredTimes(z)-z)**2
enddo
SpearmanCal=1-6*SumRankDiff/(nvalues*(nvalues**2-1))

!***validation
SumPslDiffSqur=0
meanPSLmeas=0
meanPSLSim=0
nvalues=0
covPSL=0
StDevmeas=0
StDevSim=0
OrderMeasuredPSL=MeasuredPSL
OrderSimulatedPSLMeasuredTimes=SimulatedPSLMeasuredTimes
SumRankDiff=0

do z=valstart,valfinish
 PslDiffsqure(z) = (SimulatedPSLMeasuredTimes(z)-MeasuredPSL(z))**2
 SumPslDiffSqur = SumPslDiffSqur + PslDiffsqure(z)
 meanPSLmeas=meanPSLmeas + MeasuredPSL(z)
 meanPSLSim=meanPSLSim + SimulatedPSLMeasuredTimes(z)
 nvalues=nvalues+1
enddo

!print*, meanPSLmeas, meanPSLsim, nvalues, valstart,valfinish

meanPSLmeas=meanPSLmeas/nvalues
meanPSLsim=meanPSLsim/nvalues
RMSEval=sqrt(SumPslDiffSqur/nvalues)


do z=valstart,valfinish
    covPSL=covPSL+(MeasuredPSL(z)-meanPSLmeas)*(simulatedPSLMeasuredTimes(z)-meanPSLsim)
    StDevmeas = StDevmeas + (MeasuredPSL(z)-meanPSLmeas)**2
    StDevSim = StDevSim + (simulatedPSLMeasuredTimes(z)-meanPSLsim)**2
enddo
StDevmeas = sqrt(StDevmeas)
StDevSim = sqrt(StDevSim)
PearsonVal = covPSL/(StDevmeas*StDevSim)

!print*, meanPSLmeas, meanPSLsim, RMSEval, covPSL, StDevmeas, StDevSim, PearsonVal

do z=valstart,valfinish
 do i=valstart,valfinish-1
  if (OrderMeasuredPSL(i).le.OrderMeasuredPSL(i+1)) then
    temp1 = OrderMeasuredPSL(i)
    OrderMeasuredPSL(i) = OrderMeasuredPSL(i+1)
    OrderMeasuredPSL(i+1) = temp1
    temp1 = OrderSimulatedPSLMeasuredTimes(i)
    OrderSimulatedPSLMeasuredTimes(i) = OrderSimulatedPSLMeasuredTimes(i+1)
    OrderSimulatedPSLMeasuredTimes(i+1) = temp1
    endif
 enddo
enddo
do z=valstart,valfinish
    locOrderSimulatedPSLMeasuredTimes(z)=z
enddo
do z=valstart,valfinish
 do i=valstart,valfinish-1
  if (OrderSimulatedPSLMeasuredTimes(i).le.OrderSimulatedPSLMeasuredTimes(i+1)) then
    temp1 = OrderSimulatedPSLMeasuredTimes(i)
    OrderSimulatedPSLMeasuredTimes(i) = OrderSimulatedPSLMeasuredTimes(i+1)
    OrderSimulatedPSLMeasuredTimes(i+1) = temp1
    temp1 = locOrderSimulatedPSLMeasuredTimes(i)
    locOrderSimulatedPSLMeasuredTimes(i) = locOrderSimulatedPSLMeasuredTimes(i+1)
    locOrderSimulatedPSLMeasuredTimes(i+1) = temp1
    endif
 enddo
enddo
do z=valstart,valfinish
    SumRankDiff=SumRankDiff + (locOrderSimulatedPSLMeasuredTimes(z)-z)**2
enddo
SpearmanVal=1-6*SumRankDiff/(nvalues*(nvalues**2-1))



!print*, RMSECal,PearsonCal,SpearmanCal

!**********************  
!END OF WATER TABLE Caculation
!*****************************



!**********************  
!Discharge Caculations
!**********************  

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
  enddo
   MeanMeasured=0
   MeanSimulated=0
   MeanLogMeasured=0
   MeanLogSimulated=0
   nvalues=0
   do z=calibrationTimes(1), calibrationTimes(2)
     if ((MeasuredDischarge(z)).GE.(0.0)) then
       MeanMeasured = MeanMeasured + MeasuredDischarge(z)
       MeanSimulated = MeanSimulated + SimulatedDischarge(z)
       MeanLogMeasured = MeanLogMeasured + LogMeasured(z)
       MeanLogSimulated = MeanLogSimulated + LogSimulated(z)
       nvalues=nvalues+1
     endif
   enddo
   MeanMeasured=MeanMeasured/nvalues
   MeanSimulated=MeanSimulated/nvalues
   MeanLogMeasured=MeanLogMeasured/nvalues
   MeanLogSimulated=MeanLogSimulated/nvalues
!   print*,mindischargevalues,spinupvalues,meanmeasured,meansimulated,MeasuredDisShort(1),SimulatedDisShort(1)
   SubBiasCal = 100* (MeanSimulated-MeanMeasured)/MeanMeasured
   do z=calibrationTimes(1), calibrationTimes(2)
     if ((MeasuredDischarge(z)).GE.(0.0)) then
       MeasuredSquare(z)=(MeasuredDischarge(z)-MeanMeasured)**2
       SimulatedSquare(z)=(SimulatedDischarge(z)-MeanSimulated)**2
       SimuMeasSquare(z)=(SimulatedDischarge(z)-MeasuredDischarge(z))**2
       SimMeas(z)=(MeasuredDischarge(z)-MeanMeasured)*(SimulatedDischarge(z)-MeanSimulated)
       LogMeasuredSquare(z)=(LogMeasured(z)-MeanLogMeasured)**2
       LogSimuMeasSquare(z)=(LogSimulated(z)-LogMeasured(z))**2
     else
        MeasuredSquare(z)=0
        SimulatedSquare(z)=0
        SimuMeasSquare(z)=0
        SimMeas(z)=0
        LogMeasuredSquare(z)=0
        LogSimuMeasSquare(z)=0
     endif
!  if (z.gt.10200) then 
!       print*,z,MeasuredDisShort(z),SimulatedDisShort(z),MeasuredSquare(z),SimuMeasSquare(z)
!   endif
   enddo
     stdDevMeas=sqrt(sum(MeasuredSquare(calibrationTimes(1):calibrationTimes(2)))/nvalues)
     stdDevSim=sqrt(sum(SimulatedSquare(calibrationTimes(1):calibrationTimes(2)))/nvalues)
     SumMeasuredSquare=sum(MeasuredSquare(calibrationTimes(1):calibrationTimes(2)))
     SumSimulatedSquare=sum(SimulatedSquare(calibrationTimes(1):calibrationTimes(2)))
     sumSimMeas=sum(SimMeas(calibrationTimes(1):calibrationTimes(2)))
     correlation = sumSimMeas/(sqrt(SumMeasuredSquare)*sqrt(SumSimulatedSquare))

     subKGECal = 1 - sqrt((correlation-1)**2+((stdDevSim/stdDevMeas)-1)**2+((MeanSimulated/MeanMeasured)-1)**2)

     
     SubNSECal = 1- sum(SimuMeasSquare(calibrationTimes(1): calibrationTimes(2)))/Sum(MeasuredSquare(calibrationTimes(1):calibrationTimes(2)))
     temp1=sum(LogSimuMeasSquare(calibrationTimes(1): calibrationTimes(2)))
     temp2=Sum(LogMeasuredSquare(calibrationTimes(1):calibrationTimes(2)))
     SubLogNSECal = 1- sum(LogSimuMeasSquare(calibrationTimes(1): calibrationTimes(2)))/Sum(LogMeasuredSquare(calibrationTimes(1):calibrationTimes(2)))

   MeanMeasured=0
   MeanSimulated=0
   MeanLogMeasured=0
   MeanLogSimulated=0
   nvalues=0
   do z=ValidationTimes(1), ValidationTimes(2)
     if ((MeasuredDischarge(z)).GE.(0.0)) then
       MeanMeasured = MeanMeasured + MeasuredDischarge(z)
       MeanSimulated = MeanSimulated + SimulatedDischarge(z)
       MeanLogMeasured = MeanLogMeasured + LogMeasured(z)
       MeanLogSimulated = MeanLogSimulated + LogSimulated(z)
       nvalues=nvalues+1
     endif
   enddo
   MeanMeasured=MeanMeasured/nvalues
   MeanSimulated=MeanSimulated/nvalues
   MeanLogMeasured=MeanLogMeasured/nvalues
   MeanLogSimulated=MeanLogSimulated/nvalues
!   print*,mindischargevalues,spinupvalues,meanmeasured,meansimulated,MeasuredDisShort(1),SimulatedDisShort(1)
   SubBiasVal = 100* (MeanSimulated-MeanMeasured)/MeanMeasured
   do z=ValidationTimes(1), ValidationTimes(2)
     if ((MeasuredDischarge(z)).GE.(0.0)) then
       MeasuredSquare(z)=(MeasuredDischarge(z)-MeanMeasured)**2
       SimulatedSquare(z)=(SimulatedDischarge(z)-MeanSimulated)**2
       SimuMeasSquare(z)=(SimulatedDischarge(z)-MeasuredDischarge(z))**2
       SimMeas(z)=(MeasuredDischarge(z)-MeanMeasured)*(SimulatedDischarge(z)-MeanSimulated)
       LogMeasuredSquare(z)=(LogMeasured(z)-MeanLogMeasured)**2
       LogSimuMeasSquare(z)=(LogSimulated(z)-LogMeasured(z))**2
     else
        MeasuredSquare(z)=0
        SimulatedSquare(z)=0
        SimuMeasSquare(z)=0
        SimMeas(z)=0
        LogMeasuredSquare(z)=0
        LogSimuMeasSquare(z)=0
     endif
!  if (z.gt.10200) then 
!       print*,z,MeasuredDisShort(z),SimulatedDisShort(z),MeasuredSquare(z),SimuMeasSquare(z)
!   endif
   enddo

     stdDevMeas=sqrt(sum(MeasuredSquare(ValidationTimes(1):ValidationTimes(2)))/nvalues)
     stdDevSim=sqrt(sum(SimulatedSquare(ValidationTimes(1):ValidationTimes(2)))/nvalues)
     SumMeasuredSquare=sum(MeasuredSquare(ValidationTimes(1):ValidationTimes(2)))
     SumSimulatedSquare=sum(SimulatedSquare(ValidationTimes(1):ValidationTimes(2)))
     sumSimMeas=sum(SimMeas(ValidationTimes(1):ValidationTimes(2)))
     correlation = sumSimMeas/(sqrt(SumMeasuredSquare)*sqrt(SumSimulatedSquare))

     subKGEVal = 1 - sqrt((correlation-1)**2+((stdDevSim/stdDevMeas)-1)**2+((MeanSimulated/MeanMeasured)-1)**2)
      
    SubNSEVal = 1- sum(SimuMeasSquare(ValidationTimes(1):ValidationTimes(2)))/Sum(MeasuredSquare(ValidationTimes(1):ValidationTimes(2)))
     temp1=sum(LogSimuMeasSquare(ValidationTimes(1): ValidationTimes(2)))
     temp2=Sum(LogMeasuredSquare(ValidationTimes(1):ValidationTimes(2)))
    SubLogNSEVal = 1- sum(LogSimuMeasSquare(ValidationTimes(1):ValidationTimes(2)))/Sum(LogMeasuredSquare(ValidationTimes(1):ValidationTimes(2)))
 
!**********************  
!End of Discharge Caculations
!**********************  

    
   
return

9995 write (*,*) 'Error openinig file ',precipFilename
write(*,'(''paused, type [enter] to continue'')')
read (*,*)
stop 
9996 write (*,*) 'Error openinig file ',dischargefilename
write(*,'(''paused, type [enter] to continue'')')
read (*,*)
stop 

9988 write (*,*) 'Error openinig file ',SimulatedLevelFilepath
write(*,'(''paused, type [enter] to continue'')')
read (*,*)
stop 
9989 write (*,*) 'Error openinig file ',GroundLevelFilepath
write(*,'(''paused, type [enter] to continue'')')
read (*,*)
stop 


end subroutine WriteXmlFile

    
Subroutine orderdata(NumberOfValues,OrderArray,DeepSoilCond1,ShallowSoilCond1,ShallowSoilDepth1,AePeRatio1,Urban1,NSECal,BiasCal,KGECal,LogNSECal,RMSECal,PeaCal,SpeCal,NSEVal,BiasVal,KGEVal,LogNSEVal,RMSEVal,PeaVal,Speval)

implicit none

integer,intent(in)        :: NumberOfValues  
real,intent(inout)        :: OrderArray(NumberOfValues),DeepSoilCond1(NumberOfValues),ShallowSoilCond1(NumberOfValues),ShallowSoilDepth1(NumberOfValues),AePeRatio1(NumberOfValues),Urban1(NumberOfValues)
real,intent(inout)        :: NSECal(NumberOfValues),BiasCal(NumberOfValues),KGECal(NumberOfValues),LogNSECal(NumberOfValues),RMSECal(NumberOfValues),PeaCal(NumberOfValues),SpeCal(NumberOfValues)
real,intent(inout)        :: NSEVal(NumberOfValues),BiasVal(NumberOfValues),KGEVal(NumberOfValues),LogNSEVal(NumberOfValues),RMSEVal(NumberOfValues),PeaVal(NumberOfValues),Speval(NumberOfValues)

!local
integer                  :: i,j,itemp,OrigSim(NumberOfValues)
real                     :: temp
real                     :: temparray(NumberOfValues)
    
    !simple bubble loop[
do i=1,NumberOfValues
    OrigSim(i)=i
enddo
do i=1,NumberOfValues
    do j=1,NumberOfValues-1
     if (OrderArray(j).lt.OrderArray(j+1)) then
       temp=OrderArray(j+1)
       OrderArray(j+1)=OrderArray(j)
       OrderArray(j)=temp
       itemp=OrigSim(j+1)
       OrigSim(j+1)=OrigSim(j)
       OrigSim(j)=itemp
     endif
   enddo
enddo
temparray=DeepSoilCond1
do i=1,NumberOfValues
DeepSoilCond1(i)=temparray(OrigSim(i))
enddo
temparray=ShallowSoilCond1
do i=1,NumberOfValues
ShallowSoilCond1(i)=temparray(OrigSim(i))
enddo
temparray=ShallowSoilDepth1
do i=1,NumberOfValues
ShallowSoilDepth1(i)=temparray(OrigSim(i))
enddo
temparray=AePeRatio1
do i=1,NumberOfValues
AePeRatio1(i)=temparray(OrigSim(i))
enddo
!temparray=Strickler1
!do i=1,NumberOfValues
!Strickler1(i)=temparray(OrigSim(i))
!enddo
temparray=Urban1
do i=1,NumberOfValues
Urban1(i)=temparray(OrigSim(i))
enddo
temparray=NSECal
do i=1,NumberOfValues
NSECal(i)=temparray(OrigSim(i))
enddo
temparray=BiasCal
do i=1,NumberOfValues
BiasCal(i)=temparray(OrigSim(i))
enddo
temparray=KGECal
do i=1,NumberOfValues
KGECal(i)=temparray(OrigSim(i))
enddo
temparray=LogNSECal
do i=1,NumberOfValues
LogNSECal(i)=temparray(OrigSim(i))
enddo
temparray=RMSECal
do i=1,NumberOfValues
RMSECal(i)=temparray(OrigSim(i))
enddo
temparray=PeaCal
do i=1,NumberOfValues
PeaCal(i)=temparray(OrigSim(i))
enddo
temparray=SpeCal
do i=1,NumberOfValues
SpeCal(i)=temparray(OrigSim(i))
enddo
temparray=NSEVal
do i=1,NumberOfValues
NSEVal(i)=temparray(OrigSim(i))
enddo
temparray=BiasVal
do i=1,NumberOfValues
BiasVal(i)=temparray(OrigSim(i))
enddo
temparray=KGEVal
do i=1,NumberOfValues
KGEVal(i)=temparray(OrigSim(i))
enddo
temparray=LogNSEVal
do i=1,NumberOfValues
LogNSEVal(i)=temparray(OrigSim(i))
enddo
temparray=RMSEVal
do i=1,NumberOfValues
RMSEVal(i)=temparray(OrigSim(i))
enddo
temparray=PeaVal
do i=1,NumberOfValues
PeaVal(i)=temparray(OrigSim(i))
enddo
temparray=SpeVal
do i=1,NumberOfValues
SpeVal(i)=temparray(OrigSim(i))
enddo

return

end subroutine orderdata


    
subroutine reducedurbanrain(libraryfile,filelength,text,soildetend,NoDischargeValues,urbanRainFraction,nrows,ncols,Vegtypenew)

implicit none

character(300),intent(in) :: libraryfile
integer,intent(in)        :: filelength
character(300),intent(in) :: text(filelength)
integer,intent(in)        :: soildetend,NoDischargeValues
real,intent(in)           :: urbanRainFraction
integer,intent(in)        :: nrows,ncols,Vegtypenew(nrows,ncols)

!local 
real                      :: precip(NoDischargeValues)
character(300) :: text1,precipFile,precipFilename,precipfilenameUrban,header
integer :: namelength,vegorder,i,j
real :: vegurban(nrows*ncols)

!reduced precipitation in urban areas. Seperaqte sewers  
    text1=text(soildetend+2)
    namelength=len_trim(text1)
    precipFile=text1(30:namelength-30)
    precipFilename='../'//trim(libraryfile)//'/'//trim(precipFile)
    precipfilenameUrban='../'//trim(libraryfile)//'/'//'Urban_Precip.csv'
    open(16,FILE=trim(precipFilename),err=9995,status='old')
    Read(16,*) header
    open(18,FILE=trim(precipfilenameUrban))
    write(18,'(A50)') 'Urban rainfall. Reduced rainfall in urban areas'
    vegorder=0
     do i=1,nrows
        do j=1,ncols
            if (Vegtypenew(i,j).eq.1) then
               vegorder=vegorder+1
               vegurban(vegorder)=1.0
            elseif (Vegtypenew(i,j).eq.2) then
                !urban
               vegorder=vegorder+1
               vegurban(vegorder)=urbanRainFraction
            endif
        enddo
    enddo
    do i=1,NoDischargeValues
        Read(16,*) (precip(j),j=1,vegorder)
        write(18,'(*(f6.2, ", "))') (precip(j)*vegurban(j),j=1,vegorder)
    enddo
    close(16)
    close(18)

return    
    
9995 write (*,*) 'Error openinig file ',precipFilename
write(*,'(''paused, type [enter] to continue'')')
read (*,*)
stop 


    end subroutine reducedurbanrain

    
subroutine  objectivefunction(nvalues,NSE,Bias,KGE,logNSE,RMSE,Pea,Spe,NSEFrac,BiasFrac,KGEFrac,LogNSEFrac,RMSEFrac,PeaFrac,SpeFrac,Objfn)
implicit none

integer,intent(in)        :: nvalues
real,intent(in) :: NSE(nvalues),Bias(nvalues),KGE(nvalues),LogNSE(nvalues),RMSE(nvalues),Pea(nvalues),Spe(nvalues)
real,intent(in) :: NSEFrac,BiasFrac,KGEFrac,LogNSEFrac,RMSEFrac,PeaFrac,SpeFrac
real, intent(out) :: Objfn(nvalues)
real            :: NSETemp,BiasTemp,KGETemp,LogNSETemp,RMSETemp,PeaTemp,SpeTemp
integer         :: i

NSETemp = NSEFrac / (NSEFrac+BiasFrac+KGEFrac+LogNSEFrac+RMSEFrac+PeaFrac+SpeFrac)
BiasTemp = BiasFrac / (NSEFrac+BiasFrac+KGEFrac+LogNSEFrac+RMSEFrac+PeaFrac+SpeFrac)
KGETemp = KGEFrac / (NSEFrac+BiasFrac+KGEFrac+LogNSEFrac+RMSEFrac+PeaFrac+SpeFrac)
LogNSETemp = LogNSEFrac / (NSEFrac+BiasFrac+KGEFrac+LogNSEFrac+RMSEFrac+PeaFrac+SpeFrac)
RMSETemp = RMSEFrac / (NSEFrac+BiasFrac+KGEFrac+LogNSEFrac+RMSEFrac+PeaFrac+SpeFrac)
PeaTemp = peaFrac / (NSEFrac+BiasFrac+KGEFrac+LogNSEFrac+RMSEFrac+PeaFrac+SpeFrac)
SpeTemp = SpeFrac / (NSEFrac+BiasFrac+KGEFrac+LogNSEFrac+RMSEFrac+PeaFrac+SpeFrac)

Do i = 1,nvalues
    Objfn(i)=NSETemp*NSE(i) + BiasTemp*(1/abs(Bias(i)))+ KGETemp*KGE(i) + LogNSETemp*LogNSE(i) + RMSETemp*(1/RMSE(i)) + PeaTemp*Pea(i) + SpeTemp*Spe(i)
Enddo

return


    end subroutine objectivefunction


subroutine  objectivefunction1(NSE,Bias,KGE,logNSE,RMSE,Pea,Spe,NSEFrac,BiasFrac,KGEFrac,LogNSEFrac,RMSEFrac,PeaFrac,SpeFrac,Objfn)
implicit none

real,intent(in) :: NSE,Bias,KGE,LogNSE,RMSE,Pea,Spe
real,intent(in) :: NSEFrac,BiasFrac,KGEFrac,LogNSEFrac,RMSEFrac,PeaFrac,SpeFrac
real, intent(out) :: Objfn
real            :: NSETemp,BiasTemp,KGETemp,LogNSETemp,RMSETemp,PeaTemp,SpeTemp
integer         :: i

NSETemp = NSEFrac / (NSEFrac+BiasFrac+KGEFrac+LogNSEFrac+RMSEFrac+PeaFrac+SpeFrac)
BiasTemp = BiasFrac / (NSEFrac+BiasFrac+KGEFrac+LogNSEFrac+RMSEFrac+PeaFrac+SpeFrac)
KGETemp = KGEFrac / (NSEFrac+BiasFrac+KGEFrac+LogNSEFrac+RMSEFrac+PeaFrac+SpeFrac)
LogNSETemp = LogNSEFrac / (NSEFrac+BiasFrac+KGEFrac+LogNSEFrac+RMSEFrac+PeaFrac+SpeFrac)
RMSETemp = RMSEFrac / (NSEFrac+BiasFrac+KGEFrac+LogNSEFrac+RMSEFrac+PeaFrac+SpeFrac)
PeaTemp = peaFrac / (NSEFrac+BiasFrac+KGEFrac+LogNSEFrac+RMSEFrac+PeaFrac+SpeFrac)
SpeTemp = SpeFrac / (NSEFrac+BiasFrac+KGEFrac+LogNSEFrac+RMSEFrac+PeaFrac+SpeFrac)

Objfn=NSETemp*NSE + BiasTemp*(1/abs(Bias)) + KGETemp*KGE + LogNSETemp*LogNSE + RMSETemp*(1/RMSE) + PeaTemp*Pea + SpeTemp*Spe

return


    end subroutine objectivefunction1

    
FUNCTION julian_date (yyyy, mm, dd) RESULT (julian)
IMPLICIT NONE
! converts calendar date to Julian date
! cf Fliegel & Van Flandern, CACM 11(10):657, 1968
! example: julian_date(1970,1,1)=2440588
INTEGER,INTENT(IN) :: yyyy,mm,dd
INTEGER :: julian
julian = dd-32075+1461*(yyyy+4800+(mm-14)/12)/4 + &
367*(mm-2-((mm-14)/12)*12)/12- &
3*((yyyy + 4900 + (mm - 14)/12)/100)/4
END FUNCTION julian_date