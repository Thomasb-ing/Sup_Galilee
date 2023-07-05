#CC = g++
CC = mpic++
#CC = icc
CFLAGS = -Wall -fPIC -O2 -mtune=native -fopenmp
LIBS = -lm -lrt  -lmpi 

#CFLAGS = -Wall -fPIC -O2  -openmp
#CFLAGS = -fPIC -g
#LIBS = -lrt -lblas

#LIBS = -lrt -Wl,-rpath=$(BLASDIR) -L$(BLASDIR) -lopenblas

OBJS = Chronometer.o diffusion_sequentiel.o  

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

all:	diffusion.exe

clean:
	@rm -fr *.o *.exe *~

diffusion.exe: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBS)
