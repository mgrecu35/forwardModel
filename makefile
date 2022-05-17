##  SFM  06/18/2013  Converted all libraries to those local to hrunting.
##                   Changed to Toolkit v3.31.4. Converted to task standard
##                   formats.  Added python scripts and missing libraries.
##  SFM  06/27/2013  Eliminated unused libraires
##  SFM  07/16/2013  Moved to new version Toolkit 3.40.4
##  SFM  07/17/2013  Added modules from M.Grecu
##  SFM  08/23/2013  Adapted to new gpmdev computer environment
##  SFM  10/24/2013  Moved to new version Toolkit 3.40.5
##  SFM  12/06/2013  Added header_getter code; Toolkit at 3.50.8
##  SFM  05/12/2014  Remove gfsModule.f90 & lib lnetcdff
##  SFM  06/16/2014  Substitute fModelFortran.NBF.cc for fModelFortran.cc
##  SFM  08/22/2014  Tweaked meta_assign options
##  SFM  08/22/2014  Added getbits.c module

#---Compilers
CC =gcc -fPIC -g
CPP=gcc -fPIC -g
F90=gfortran -fPIC
F77=gfortran -fPIC
#F90=ifort -fPIC -fstack-security-check -traceback -CB -ftrapuv -check uninit
#F77=ifort -shared-intel -fPIC -fstack-security-check -g -traceback -CB -ftrapuv -check uninit

#---Compiler Options
OPENMPF =
OPENMPG =
FOPTS   = -c -O2 $(OPENMPF)
F77OPTS = -c -O2 $(OPENMPF)
CFLAGS  = -O2 $(OPENMPF)
C++FLAGS= -O2 $(OPENMPG)

#---Compiler Directives
CPP_SWITCH = IFORT
CC_SWITCH  = IFORT

#---Specialty Library Locations
HDF4_INC=/usr/local/HDF4/include/
HDF5_INC=/usr/local/HDF5/include/
#HDF4_LIB=/usr/local/HDF4/lib
#HDF5_LIB=/usr/local/HDF5/lib
#TKIO=/PANFS/user/home/tkio/tkio-3.60.4
#---Required include files
HDF4_INC = /Users/mgrecu/miniconda3/include
HDF5_INC = /Users/mgrecu/miniconda3/include
INC=    -I$(HDF4_INC)\
        -I$(HDF5_INC)\
	-Iflann/include/\
        -I$(TKIO)/inc/ccode\
        -I/usr/local/include/\
        -I/home/grecu/cmbv6x/multiscatter-1.2.10/include

#---Required library files
LIB =	-L$(HDF5_LIB) \
	-L$(HDF4_LIB) \
	-L$(TKIO)/lib \
	-Llapack-3.4.1\
	-Lflann/lib -L/home/grecu/cmbv6x/multiscatter-1.2.10/lib \
        -L/PANFS/local/intel/Compiler/14.0.1/lib/intel64/ \

#---List of Normal Libraires to Search for References during Linkage

LIBES = -lmultiscatter -llapack -lblas


#---Compile Stages
OBJ =	f90DataTypes.o  random.o writegif.o nbin.o\
	weightModule.o lut_routines_sjm.o\
	emissModule.o missing_flags.o read_melt_percentages.o\
	globModule.o readTables_nonsph.o cloud.o\
	absorb.o gcloud.o  \
	geomask.o f90Types.o allocateMem.o ezlhconv.o readNh.o beamConvP.o beamConvSet.o\
	resampGMI.o interpol.o geophysEns.o gEnv.o outputminmax.o \
        dboux.o  fhb1_noNaN.o \
        retTablesInt.half.o multiscatter2_ascii.o band_dble.o\
	rterain.o \
        convPix2.py.o latlon.o setOpt.o convAllFreq.o asciiCon.o writeEnKF.o screenocean.o\
	filterUp1.o\
        kgainBLAS.o getScattProp.o armadillo_funcs.o

OBJEXE	=    $(OBJ)  mainfpy.o

#---Link and create Executable
combAlgLib:		$(OBJEXE)
	gcc -shared -o combAlgLib.so $(OBJEXE)

combAlg:	combAlgLib.so
		f2py -c -m combAlg --f90flags="-fopenmp" --arch="-fopenmp"\
		src_f90/fobjective.f90 src_f90/mainfpy.f90 \
                src_f90/hbprof_new.f90 \
		src_f90/bisection.f90  src_f90/simpleFModel.f90 \
                src_f90/forwardModel.f90 src_f90/emissivity-sp.f src_f90/radtran_tau_dble.f src_f90/rosen.f src_f90/multiscatter.f90\
                src_f90/absorption3D.f90 $(OBJ)  $(LIB) $(LIBES) -lstdc++

#---Auxiliary Functions
clean:
	rm -f *.o combAlg.exe *.mod nohup* intern*dat DUMMY*

#---Individual Object Modules Make-Formulae
mainfpy.o:	src_f90/mainfpy.f90
		$(F90) -c -fPIC -O2 src_f90/mainfpy.f90
read_melt_percentages.o:        src_f90/read_melt_percentages.f90
		$(F90) -c -O2 src_f90/read_melt_percentages.f90

interp_melt_percentages.o:      MeltFrac/interp_melt_percentages.f90
		$(F90) -c -O2 MeltFrac/interp_melt_percentages.f90

sst.rend.o:     src_f90/sst.rend.f90 src_f90/f90DataTypes.f90
		$(F90) -c  -O2 src_f90/sst.rend.f90

f90DataTypes.o:	src_f90/f90DataTypes.f90
		$(F90) -c  -O2 src_f90/f90DataTypes.f90

f90Types.o:	src_f90/f90Types.f90
		$(F90) -c  -O2 src_f90/f90Types.f90

pcaModule.o:	pcaModule.f90
		$(F90) -c  -O2 pcaModule.f90

outputminmax.o:	src_f90/outputminmax.f90
		$(F90) -c  -O2 src_f90/outputminmax.f90

l1l1l1.o:       l1l1l1.f90
		$(F90) -c  -O2 l1l1l1.f90

asciiCon.o:	src_f90/asciiCon.f90
		$(F90) -c  -O2 src_f90/asciiCon.f90

rteModule.o:	src_f90/rteModule.f90
		$(F90) -c  -O2 src_f90/rteModule.f90

armadillo_funcs.o:
		g++  -std=c++11 -c -fPIC -DARMA_DONT_USE_WRAPPER -I./src_c/armadillo-10.1.2/include src_c/armadillo_funcs.cpp

rainDisc.o:	rainDisc.f90
		$(F90) -c  -O2 rainDisc.f90

gEnv.o:		src_f90/gEnv.f90
		$(F90) -c  -O2 src_f90/gEnv.f90

pia35fromTbs.o:	pia35fromTbs.f
		$(F90) -c  -O2 pia35fromTbs.f

clearSky.Ocean.o:	clearSky.Ocean.f
		$(F90) -c  -O2 clearSky.Ocean.f

clearSkyTb.Land.o:	clearSkyTb.Land.f
		$(F90) -c  -O2 clearSkyTb.Land.f

sfcTpw.o:	sfcTpw.f
		$(F90) -c  -O2 sfcTpw.f

ffnet.o:	ffnet.f
		$(F90) -c  -O2 ffnet.f

ezlhconv.o:	src_f90/ezlhconv.f90
		$(F90) -c -O2 src_f90/ezlhconv.f90

readNh.o:	src_c/readNh.c
		$(CC) -c -O2 src_c/readNh.c

flannModule.o:	src_c/flannModule.c
		$(CC) $(INC) -c -O2 src_c/flannModule.c

dboux.o:	src_c/dboux.c
		$(CC) -c -O2 src_c/dboux.c

cmbAlg.o:	src_c/cmbAlg.c
		$(CC) -c -O2 src_c/cmbAlg.c

writeA.o:	writeA.c
		$(CC) -c -O2 writeA.c

getbits.o:	src_c/getbits.c
		$(CC) -c -O2 src_c/getbits.c

allocateMem.o:	src_f90/allocateMem.f90 \
		src_f90/f90Types.f90 src_f90/f90DataTypes.f90
		$(F90) -c  -O2 src_f90/allocateMem.f90

lbfgs.o:	toms778_f90/lbfgs.f90
		$(F90) -c -O2 toms778_f90/lbfgs.f90

emissModule.o:	src_f90/emissModule.f90
		$(F90) -c  -O2 $(OPENMP) src_f90/emissModule.f90

filterUp1.o:	src_f90/filterUp1.f90
		$(F90) -c  -O2  src_f90/filterUp1.f90

filterUp1.exp1.o:	filterUp1.exp1.f90
		$(F90) -c  -O2  filterUp1.exp1.f90

#kgainBLAS.Feb2015.o:	kgainBLAS.Feb2015.f90
#		$(F90) -c  -O2 kgainBLAS.Feb2015.f90

kgainBLAS.o:	src_f90/kgainBLAS.f90
		$(F90) -c  -O2 src_f90/kgainBLAS.f90

simCSTbs2.o:	simCSTbs2.f90
		$(F90) -c  -O2 simCSTbs2.f90

radarRet.o:	src_f90/radarRet.f90
		$(F90) -c  -O0 src_f90/radarRet.f90

radarRet_P1.o:	src_f90/radarRet_P1.f90
		$(F90) -c  -O0 src_f90/radarRet_P1.f90

radarRet_P2.o:	src_f90/radarRet_P2.f90
		$(F90) -c  -O0 -openmp src_f90/radarRet_P2.f90

radarRet_P4.o:	src_f90/radarRet_P4.f90 rteModule.o
		$(F90) -c  -O0 -openmp src_f90/radarRet_P4.f90

radarRet_P4_FS.o:	src_f90/radarRet_P4_FS.f90 rteModule.o
		$(F90) -c  -O0 -openmp src_f90/radarRet_P4_FS.f90

radarRet.Empty.o:	src_f90/radarRet.Empty.f90
		$(F90) -c  -O2  src_f90/radarRet.Empty.f90

convPix.o:	convPix.f90
		$(F90) -c  -O2 convPix.f90

globModule.o:	src_f90/globModule.f90
		$(F90) -c  -O2 src_f90/globModule.f90

writeCMB.rend.o:	writeCMB.rend.c
		$(CC) -c  -O2  $(INC) writeCMB.rend.c

random.o:	src_f90/random.f90
		$(F90) -c  -O2 src_f90/random.f90

emissEstim.o:	emissEstim.f90
		$(F90) -c  $(OPENMPF) -O2 emissEstim.f90

nbin.o:		src_f90/nbin.f90
		$(F90) -c  -O2 src_f90/nbin.f90

#begin  WSO 9/14/13 add missing flags module
missing_flags.o:	src_f90/missing_flags.f90
		$(F90) -c -O2 src_f90/missing_flags.f90
#end    WSO 9/14/13

main.o:		src_c/main.c
		$(CC) -c   -O2 src_c/main.c

inputParser.o:	src_c/inputParser.c
		$(CC) $(INC) -c -O2 src_c/inputParser.c

writeDiagnost.o:	writeDiagnost.c
		$(CC) $(INC) -c -O2 writeDiagnost.c

allocMemC.o:	RadarOnly/allocMemC.c
		$(CC) -D$(CC_SWITCH) -c -O2 RadarOnly/allocMemC.c

readAPriori.o:	RadarOnly/readAPriori.cc
		$(CPP) -c -D$(CPP_SWITCH)  $(C++FLAGS)	RadarOnly/readAPriori.cc

readHDF5.o:	readHDF5.c
		$(CC) $(INC) -c -O2 readHDF5.c

writeEnKF.o:	src_c/writeEnKF.c
		$(CC) $(INC) -c -O2 src_c/writeEnKF.c

weightModule.o:	src_f90/weightModule.f90
		$(F90) -c -O2 src_f90/weightModule.f90

convPix2.py.o:	src_f90/convPix2.py.f90
		$(F90) -c -O2  src_f90/convPix2.py.f90

getScattProp.o:	src_f90/getScattProp.f90
		$(F90) -c -O2  src_f90/getScattProp.f90

deconvTb2.o:	lbfgs.o	Fortran/deconvTb2.f90
		$(F90) -c -O2  Fortran/deconvTb2.f90

convAllFreq.o:	src_f90/convAllFreq.f90
		$(F90) -c -O2 src_f90/convAllFreq.f90

setOpt.o:	src_f90/setOpt.f90
		$(F90) -c -O2 src_f90/setOpt.f90

latlon.o:	src_f90/latlon.f90
		$(F90) -c -O2 src_f90/latlon.f90

mainf.o:	src_f90/mainf.f90 src_f90/f90DataTypes.f90
		$(F90) -c -O2 $(OPENMPF) src_f90/mainf.f90

get1Cprof.o:	get1Cprof.f90
		$(F90) $(OPENMPF) -c -O2 get1Cprof.f90

binaryOut.o:	binaryOut.c
		$(CC) -c -O2 binaryOut.c

multiscatter2_ascii.o:	src_c/multiscatter2_ascii.c
		$(CC) -I/home/grecu/cmbv6x/multiscatter-1.2.10/include -c -O2 src_c/multiscatter2_ascii.c

best2.o:	best2.c
		$(CC) -c -O2 best2.c

best2st.o:	best2st.c
		$(CC) -c -O2 best2st.c

pca.o:	pca.c
		$(CC) -c -O2 pca.c

interpol.o:	src_f90/interpol.f90
		$(F90) -c  -O2 $(FOPTS) src_f90/interpol.f90

beamConvP.o:	src_f90/beamConvP.f90
		$(F90) -c   -O2 $(FOPTS) src_f90/beamConvP.f90

beamConvSet.o:	src_f90/beamConvSet.f90
		$(F90) -c   -O2 $(FOPTS) src_f90/beamConvSet.f90

mie3.o:		RadarOnly/Mie/mie3.f90
		$(F90) $(FOPTS) RadarOnly/Mie/mie3.f90

ensFilter.o:	src_c/ensFilter.cc src_c/typedef.h
		$(CPP) -c -D$(CPP_SWITCH)  $(C++FLAGS)  \
		src_c/ensFilter.cc

clearSky.o:	RadarOnly/clearSky.cc RadarOnly/typedef.h
		$(CPP) -c -D$(CPP_SWITCH)  $(C++FLAGS) $(OPENMPC) \
		RadarOnly/clearSky.cc

#  SFM  start  06/16/2014 - added for M.Grecu
fModelFortran.NUBF.o:	src_c/fModelFortran.NUBF.cc
		$(CPP) -c -D$(CPP_SWITCH) $(C++FLAGS) \
		src_c/fModelFortran.NUBF.cc

fModelFortranCv.o:	RadarOnly/fModelFortranCv.cc
		$(CPP) -c -D$(CPP_SWITCH) $(C++FLAGS) \
		RadarOnly/fModelFortranCv.cc
#  SFM  end  06/16/2014

linsys.o	: RadarOnly/linsys.f90
		$(F90) $(FOPTS)  RadarOnly/linsys.f90

band_dble.o	: src_f90/band_dble.f90
		$(F90) $(FOPTS) $(OPENMP) src_f90/band_dble.f90

dsvdc.o	: 	RadarOnly/dsvdc.f90
		$(F90) $(FOPTS)  RadarOnly/dsvdc.f90

depAgg.o	: 	RadarOnly/depAgg.f90
		$(F90) $(FOPTS)  RadarOnly/depAgg.f90

writegif.o:	src_f90/writegif.f90
		$(F90) $(FOPTS)  src_f90/writegif.f90

cloud.o:	src_f90/cloud.f90
		$(F90) $(FOPTS)  src_f90/cloud.f90

absorb.o:	src_f90/absorb.f
		$(F90) $(F77OPTS)  src_f90/absorb.f

sort3a.o:	sort3a.f
		$(F90) $(F77OPTS)  sort3a.f

emissivity-sp.o:	RadarOnly/RTE/emissivity-sp.f
		$(F90) $(F77OPTS) $(OPENMP) RadarOnly/RTE/emissivity-sp.f

radtran_tau_dble.o:	src_f90/radtran_tau_dble.f
		$(F90) $(F77OPTS) $(OPENMP) src_f90/radtran_tau_dble.f

fhb1_noNaN.o:		src_f90/fhb1_noNaN.f90
		$(F90) $(F77OPTS)  $(OPENMP) src_f90/fhb1_noNaN.f90

bisection.o:	src_f90/bisection.f90
		$(F90) $(F77OPTS)  src_f90/bisection.f90

mie-melt_lam.o:	RadarOnly/MieStuff/mie-melt_lam.f
		$(F90) $(F77OPTS)  RadarOnly/MieStuff/mie-melt_lam.f

iceoptic.o:	RadarOnly/MieStuff/iceoptic.f
		$(F90) $(F77OPTS)  RadarOnly/MieStuff/iceoptic.f

watoptic.o:	RadarOnly/MieStuff/watoptic.f
		$(F90) $(F77OPTS)  RadarOnly/MieStuff/watoptic.f

mg-ell.o:	RadarOnly/MieStuff/mg-ell.f
		$(F90) $(F77OPTS)  RadarOnly/MieStuff/mg-ell.f

mie-sphere.o:	RadarOnly/MieStuff/mie-sphere.f
		$(F90) $(F77OPTS)  RadarOnly/MieStuff/mie-sphere.f

epswater.o:	RadarOnly/MieStuff/epswater.f
		$(F90) $(F77OPTS)  RadarOnly/MieStuff/epswater.f

melting.bauer.o:	RadarOnly/MieStuff/melting.bauer.f90
		$(F90) $(F77OPTS) RadarOnly/MieStuff/melting.bauer.f90

readTables_nonsph.o:	src_f90/readTables_nonsph.f90
		$(F90) $(F77OPTS) src_f90/readTables_nonsph.f90

retTablesInt.half.o:	src_f90/retTablesInt.half.f90 \
		src_f90/readTables_nonsph.f90
		$(F90) $(F77OPTS) $(OPENMP) src_f90/retTablesInt.half.f90

retTables2.o:	RadarOnly/retTables2.f90 \
		RadarOnly/readTables.Aug2015.f90
		$(F90) $(F77OPTS) $(OPENMP) RadarOnly/retTables2.f90

rterain.o:	gEnv.o src_f90/rterain.f90
		$(F90) $(FOPTS) $(OPENMP) src_f90/rterain.f90

rosen.o:	src_f90/rosen.f
		$(F90) $(FOPTS) $(OPENMP) src_f90/rosen.f

gaussNewton.o:	RadarOnly/gaussNewton.cc
		$(CC) -c -O2  RadarOnly/gaussNewton.cc

io_subs.o:	src_c/io_subs.c
		$(CC) -c  -D$(CC_SWITCH) -O2 $(INC) src_c/io_subs.c

read_GPM.o:	src_c/read_GPM.c
		echo $(INC)
		$(CC) -c  -D$(CPP_SWITCH) $(INC) src_c/read_GPM.c

read_GPM_FS.o:	src_c/read_GPM_FS.c
		echo $(INC)
		$(CC) -c  -D$(CPP_SWITCH) $(INC) src_c/read_GPM_FS.c

read_TRMM.o:	src_c/read_TRMM.c
		echo $(INC)
		$(CC) -c  $(INC) -D$(CPP_SWITCH) $(INC) src_c/read_TRMM.c

readENV.o:	src_c/readENV.c
		$(CC) -c -D$(CPP_SWITCH) $(INC) src_c/readENV.c

writeCMBOut_TK.o:	src_c/writeCMBOut_TK.c
		$(CC) -c -g -D$(CPP_SWITCH) -O2 $(INC) src_c/writeCMBOut_TK.c

writeCMBTOut_TK.o:	src_c/writeCMBTOut_TK.c
		$(CC) -c -g -D$(CPP_SWITCH) -O2 $(INC) src_c/writeCMBTOut_TK.c

writeCMBOut_TK_FS.o:	src_c/writeCMBOut_TK_FS.c
		$(CC) -c -g -D$(CPP_SWITCH) -O2 $(INC) src_c/writeCMBOut_TK_FS.c

screenocean.o:	src_c/screenocean.c
		$(CC) -c -O2 src_c/screenocean.c

dsphere2.o:	dsphere2.f
		$(F77) -c -O2 dsphere2.f

gcloud.o:	src_f90/gcloud.f
		$(F77) -c -O2 src_f90/gcloud.f

resampGMI.o:	src_f90/resampGMI.f90
		$(F90) -c $(FOPTS)  src_f90/resampGMI.f90

dPRgrid.o:	dPRgrid.f90
		$(F90) -c $(FOPTS) dPRgrid.f90

geophysEns.o:	src_f90/geophysEns.f90
		$(F90) -c $(FOPTS) src_f90/geophysEns.f90

read_ims.o:	read_ims.c
		$(CC) -c -O2 -Wall read_ims.c

autosnow_retrieval.o:	autosnow_retrieval.f90
	$(F77) -c -O2 autosnow_retrieval.f90

read_reynolds.o:read_reynolds.f90
		$(F77) -c -O2 read_reynolds.f90

meta_assign.o:	src_c/meta_assign.c
		$(CC) -c -O2 $(INC) src_c/meta_assign.c

geomask.o:	src_c/geomask.c
		$(CC) -c -O2 $(INC) src_c/geomask.c

#  SFM  start  12/06/2013 - added
header_getter_L2.o: src_f90/header_getter_L2.f90
		$(F77) -c -O2 $(INC) src_f90/header_getter_L2.f90
#  SFM  end    12/06/2013

file_check.o:	src_f90/file_check.f90
	$(F90) -c -O2 $(INC) src_f90/file_check.f90

lut_routines_sjm.o: src_f90/lut_routines_sjm.f90
		$(F90) -c $(FOPTS) src_f90/lut_routines_sjm.f90
