MAKE = make		#change this line if you are using a different GNU make software

dirFASTAtoRF = ./src/FASTAtoRF
dirRCADEEM = ./src/RCADEEM
dirRC = ./src/RC

all: MK_dir CC_FASTAtoRF CC_RCADEEM CC_RC RM_objectFiles

MK_dir:
	mkdir -p ./bin

CC_FASTAtoRF: $(dirFASTAtoRF)/Makefile
	$(MAKE) -C $(dirFASTAtoRF)
	
CC_RCADEEM: $(dirRCADEEM)/Makefile
	$(MAKE) -C $(dirRCADEEM)

CC_RC: $(dirRC)/Makefile
	$(MAKE) -C $(dirRC)

RM_objectFiles:
	rm -f $(dirFASTAtoRF)/*.o $(dirRCADEEM)/*.o $(dirRC)/*.o 

clean:
	rm -f $(dirFASTAtoRF)/*.o $(dirRCADEEM)/*.o $(dirRC)/*.o ./bin/*
