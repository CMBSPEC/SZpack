all:
	python setup.py build_ext --inplace
#	python setup.py install
#	pip install .

clean:
	rm -f *.so
	rm -rf ./build

cleanall:
	rm -f *.so
	rm -f  ./SZpack/*.so
	rm -rf ./build
	rm -rf ./dist
	rm -rf ./SZpack.egg-info
	rm -rf ./SZpack/__pycache__ ./SZpack/build