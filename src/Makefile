# Makefile for phylocom

CC = gcc -Wall -g -O1
HDRS	= phylocom.h nrutil.h ecovolve.h stats.h
OBJS	= main.o io.o new2fy.o aot.o traits.o comstruct.o nrutil.o fy2new.o bladj.o comnode.o combase.o prune.o comtrait.o stats.o
OBJSP	= phylomatic.o io.o new2fy.o nrutil.o fy2new.o prune.o combase.o
OBJSE	= ecomain.o prune.o nrutil.o fy2new.o combase.o
LIBS	= -lm


all: phylocom ecovolve phylomatic

phylocom: $(OBJS)
	$(CC) -o phylocom $(OBJS) $(LIBS)

main.o: main.c $(HDRS)

io.o: io.c $(HDRS)

fy2new.o: fy2new.c $(HDRS)

new2fy.o: new2fy.c $(HDRS)

aot.o: aot.c $(HDRS)

traits.o: traits.c $(HDRS)

comstruct.o: comstruct.c $(HDRS)

combase.o: combase.c $(HDRS)

bladj.o: bladj.c $(HDRS)

comnode.o: comnode.c $(HDRS)

nrutil.o: nrutil.c $(HDRS)

stats.o: stats.c $(HDRS)

comtrait.o: comtrait.c $(HDRS)

ecovolve: $(OBJSE)
	$(CC) -o ecovolve $(OBJSE) $(LIBS)

ecomain.o: ecomain.c  $(HDRS)

prune.o: prune.c $(HDRS)

phylomatic: $(OBJSP)
	$(CC) -o phylomatic $(OBJSP) $(LIBS)

phylomatic.o: phylomatic.c $(HDRS)

clean: 	
	rm *.o


# Compiling for Windows:
#   yum install mingw32-gcc
#   make w32
#   file phylocom.exe
#   wineconsole --backend=user ~/.wine/drive_c/windows/system32/cmd.exe
#   Z:\..\..> phylocom.exe

w32:
	mkdir w32
	cp *.[hc] Makefile w32
	# Fix lack of random and srandom in W32 
	sed -i -e 's/random\ *(/rand(/g' w32/*.c
	# Diffent compiler chain, replace with appropriate
	sed -i -e 's/gcc/i686-pc-mingw32-gcc/g' w32/Makefile
	cd w32 && $(MAKE)
	cp -f w32/phylocom phylocom.exe
	cp -f w32/phylomatic phylomatic.exe
	cp -f w32/ecovolve ecovolve.exe
	rm -rf w32
