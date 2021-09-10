gg = gfortran

he = Helium
all: $(he)

$(he): %: %.o 
	$(gg) -o $@.exe $<
%.o: %.f90
	$(gg) -c $<

clean:
	rm -f *.out *.exe *.o *.mod *.a
