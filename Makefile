# all:
# f90comp = gfortran

exe_name = test 
exe_obj = test.o 
exe_file = test.f90

path = .
src = $(path)/src
mod = $(path)/mod
obj = $(path)/obj
inc = $(path)/inc
case = $(path)/case
debug = $(path)/debug
objects = $(obj)/utils.o $(obj)/simulator.o \
		$(obj)/de_base.o $(obj)/sade.o \
		$(obj)/$(exe_obj)
compile = gfortran -J $(mod) -I $(inc)


# executable
$(debug)/$(exe_name): $(objects)
	$(compile) -o $(debug)/$(exe_name) $(objects)


# mod
$(mod)/random_numbers.mod: $(obj)/utils.o $(src)/utils.f90
	$(compile) -c $(src)/utils.f90 -o $(obj)/utils.o

$(mod)/utility.mod: $(obj)/utils.o $(src)/utils.f90
	$(compile) -c $(src)/utils.f90 -o $(obj)/utils.o

$(mod)/types.mod: $(obj)/simulator.o $(src)/simulator.f90
	$(compile) -c $(src)/simulator.f90 -o $(obj)/simulator.o

$(mod)/nosplit_simulator.mod: $(obj)/simulator.o $(src)/simulator.f90
	$(compile) -c $(src)/simulator.f90 -o $(obj)/simulator.o

$(mod)/de_base.mod: $(obj)/de_base.o $(src)/de_base.f90
	$(compile) -c $(src)/de_base.f90 -o $(obj)/de_base.o

$(mod)/sade.mod: $(obj)/sade.o $(src)/sade.f90
	$(compile) -c $(src)/sade.f90 -o $(obj)/sade.o


# obj
$(obj)/utils.o: $(src)/utils.f90
	$(compile) -c $(src)/utils.f90 -o $(obj)/utils.o

$(obj)/simulator.o: $(src)/simulator.f90
	$(compile) -c $(src)/simulator.f90 -o $(obj)/simulator.o

$(obj)/de_base.o: $(src)/de_base.f90
	$(compile) -c $(src)/de_base.f90 -o $(obj)/de_base.o

$(obj)/sade.o: $(src)/sade.f90
	$(compile) -c $(src)/sade.f90 -o $(obj)/sade.o

$(obj)/$(exe_obj): $(case)/$(exe_file)
	$(compile) -c $(case)/$(exe_file) -o $(obj)/$(exe_obj)
