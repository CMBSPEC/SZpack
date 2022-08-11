#============================================================================================
# definitions
#============================================================================================
DEV_DIR = ./Development

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

SZpack.py:
	cd ./python; make;
	cp ./python/SZpack.py .

#============================================================================================
# rules to clean up
#============================================================================================
clean:
	rm -f *.o

cleanall:
	rm -f ./src/*.o ./src/*~
	rm -f *.o *~ run_SZpack libSZpack.a 

cleanallDEV: cleanall
	rm -f $(DEV_DIR)/*.o $(DEV_DIR)/*.~
	rm -f $(DEV_DIR)/Definitions/*.o $(DEV_DIR)/Definitions/*.~
	rm -f $(DEV_DIR)/Simple_routines/*.o $(DEV_DIR)/Simple_routines/*.~
	rm -f $(DEV_DIR)/Integration/*.o $(DEV_DIR)/Integration/*~

cleanpy:
	cd python; make clean;
	rm -f SZpack.py*

tidy: cleanallDEV cleanpy

wipeDS:
	find . -type f -name \.DS_Store -print | xargs rm

#============================================================================================
