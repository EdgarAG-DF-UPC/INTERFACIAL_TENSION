objects = constants.o parametres.o distscompbc.o potencials.o integradors.o tensio.o altres.o
compilador = gfortran
OPT = -O3



#############################################################################################################
radius.x: radius.o $(objects)
	$(compilador) -o radius.x $(OPT) radius.o $(objects) 

spherical_tension.x: spherical_tension.o $(objects)
	$(compilador) -o spherical_tension.x $(OPT) spherical_tension.o $(objects) -fopenmp

planar_tension.x: planar_tension.o $(objects)
	$(compilador) -o planar_tension.x $(OPT) planar_tension.o $(objects) -fopenmp

#############################################################################################################

#############################################################################################################
radius.o: constants.mod
	$(compilador) -c $(OPT) radius.f90

spherical_tension.o: constants.mod parametres.mod spherical_tension.f90
	$(compilador) -c $(OPT) spherical_tension.f90 -fopenmp

planar_tension.o: constants.mod parametres.mod planar_tension.f90
	$(compilador) -c $(OPT) planar_tension.f90 -fopenmp

writepots.o: constants.mod parametres.mod writepots.f90
	$(compilador) -c $(OPT) writepots.f90

distscompbc.o: distscompbc.f90
	$(compilador) -c $(OPT) distscompbc.f90

potencials.o: potencials.f90 constants.mod parametres.mod
	$(compilador) -c $(OPT) potencials.f90

altres.o: altres.f90
	$(compilador) -c $(OPT) altres.f90

integradors.o: integradors.f90
	$(compilador) -c $(OPT) integradors.f90

tensio.o: constants.mod tensio.f90
	$(compilador) -c $(OPT) tensio.f90

constants.o: constants.f90
	$(compilador) -c $(OPT) constants.f90

parametres.o: parametres.f90
	$(compilador) -c $(OPT) parametres.f90
#############################################################################################################

#############################################################################################################
moduls: parametres.mod constants.mod

constants.mod: constants.f90
	$(compilador) -c $(OPT) constants.f90

parametres.mod: parametres.f90
	$(compilador) -c $(OPT) parametres.f90
#############################################################################################################

.PHONY: clean
clean:
	rm -f *.o
	rm -f *mod
	rm -f *.x

