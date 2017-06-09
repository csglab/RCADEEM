TARGET = ../../bin/RCADEEM
CPP01 = commandline
CPP02 = compare_PWMs
CPP03 = enrichment
CPP04 = hypergeometric
CPP05 = main
CPP06 = misc
CPP07 = motif_logo
CPP08 = open_read
CPP09 = optimize
CPP10 = postscript
CPP11 = PWM
CPP12 = read_fasta
CPP13 = read_rndforest
CPP14 = read_weights
CPP15 = scan
CPP16 = write

H01 = declarations
H02 = structures
H03 = commandline
H04 = graphics

OBJS = $(CPP01).o $(CPP02).o $(CPP03).o $(CPP04).o $(CPP05).o $(CPP06).o $(CPP07).o $(CPP08).o $(CPP09).o $(CPP10).o $(CPP11).o $(CPP12).o $(CPP13).o $(CPP14).o $(CPP15).o $(CPP16).o
CC = g++	# change this line if you are using a different compiler

$(TARGET): $(OBJS)
	$(CC) -O3 -o $(TARGET) $(OBJS)

$(CPP01).o: $(CPP01).cpp
	$(CC) -O3 -c $(CPP01).cpp

$(CPP02).o: $(CPP02).cpp
	$(CC) -O3 -c $(CPP02).cpp

$(CPP03).o: $(CPP03).cpp
	$(CC) -O3 -c $(CPP03).cpp

$(CPP04).o: $(CPP04).cpp
	$(CC) -O3 -c $(CPP04).cpp

$(CPP05).o: $(CPP05).cpp
	$(CC) -O3 -c $(CPP05).cpp

$(CPP06).o: $(CPP06).cpp
	$(CC) -O3 -c $(CPP06).cpp

$(CPP07).o: $(CPP07).cpp
	$(CC) -O3 -c $(CPP07).cpp

$(CPP08).o: $(CPP08).cpp
	$(CC) -O3 -c $(CPP08).cpp

$(CPP09).o: $(CPP09).cpp
	$(CC) -O3 -c $(CPP09).cpp

$(CPP10).o: $(CPP10).cpp
	$(CC) -O3 -c $(CPP10).cpp

$(CPP11).o: $(CPP11).cpp
	$(CC) -O3 -c $(CPP11).cpp

$(CPP12).o: $(CPP12).cpp
	$(CC) -O3 -c $(CPP12).cpp

$(CPP13).o: $(CPP13).cpp
	$(CC) -O3 -c $(CPP13).cpp

$(CPP14).o: $(CPP14).cpp
	$(CC) -O3 -c $(CPP14).cpp

$(CPP15).o: $(CPP15).cpp
	$(CC) -O3 -c $(CPP15).cpp

$(CPP16).o: $(CPP16).cpp
	$(CC) -O3 -c $(CPP16).cpp


$(OBJS): $(H01).h $(H02).h $(H03).h $(H04).h
