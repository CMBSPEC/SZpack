#============================================================================================
# definitions
#============================================================================================
DEV_DIR = ./Tools

MAKE_COMMAND = make -f Makefile.in

#============================================================================================
# rules
#============================================================================================
all: 
	$(MAKE_COMMAND) all

lib: 
	$(MAKE_COMMAND) lib

bin: 
	$(MAKE_COMMAND) bin
	
pub:
	make -f Makefile.pub all

pubtar:
	make -f Makefile.pub tarball

pySZpack:
	cd ./python; make;

#============================================================================================
# rules to clean up
#============================================================================================
.PHONY: clean cleanall cleanallDev cleanpy cleanallpy tidy wipeDS
clean:
	rm -f *.o

cleanall:
	rm -f ./src/*.o ./src/*~
	rm -f *.o *~ run_SZpack run_SZ_moment_method run_multiple_scatt libSZpack.a 

cleanallDEV: cleanall
	rm -f $(DEV_DIR)/*.o $(DEV_DIR)/*.~
	rm -f $(DEV_DIR)/Definitions/*.o $(DEV_DIR)/Definitions/*.~
	rm -f $(DEV_DIR)/Simple_routines/*.o $(DEV_DIR)/Simple_routines/*.~
	rm -f $(DEV_DIR)/Integration/*.o $(DEV_DIR)/Integration/*~
	rm -f $(DEV_DIR)/Output/*.o $(DEV_DIR)/Output/*~
	rm -f $(DEV_DIR)/Cosmology/*.o $(DEV_DIR)/Cosmology/*~

cleanpy:
	cd ./python; make clean;

cleanallpy:
	cd ./python; make cleanall;

tidy: cleanallDEV cleanallpy

wipeDS:
	find . -type f -name \.DS_Store -print | xargs rm

#============================================================================================
