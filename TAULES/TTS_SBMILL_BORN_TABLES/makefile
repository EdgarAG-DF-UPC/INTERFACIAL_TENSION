tts.born.x: Ulong.o Ushort.o
	gfortran TTS_Born.f90 -o tts.born.x -O3 Ulong.o Ushort.o
Ulong.o: Ulong.f90
	gfortran -c -O3 Ulong.f90 
Ushort.o: Ushort.f90
	gfortran -c -O3 Ushort.f90
sladek.born.x:
	gfortran Sladek_Born.f90 -o sladek.born.x


.PHONY: clean move
clean:
	rm -f *.o
	rm -f *mod
	rm -f *eam.alloy
	rm -f *table
	rm -f *x
move:
	mv *HeHe* He-He
	mv *LiHe* Li-He/
	mv *NaHe* Na-He/
	mv *K#He* K-He/
	mv *RbHe* Rb-He/
	mv *CsHe* Cs-He/
