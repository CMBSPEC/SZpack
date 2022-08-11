#===================================================================================================
# rules 
#===================================================================================================

all: clean SZpack.py

SZpack.py: libSZpack.a
	     @echo "\n Creating SZpack Python routines\n" 
		 python setup.py build
		 python setup.py install --user
	     
libSZpack.a:
	     cd ../; make lib;
	     
clean:
	     rm -rf build
		 rm -f SZpack.py SZpack.pyc SZpack_wrap.cpp *.pyc
	     
cleanall: clean
	     cd ..; make tidy;
	      
#===================================================================================================
#===================================================================================================
