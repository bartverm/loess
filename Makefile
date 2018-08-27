all: loess.mexa64

ifndef MATLAB_PATH
MATLAB_PATH=/usr/local/MATLAB/R2018a
endif

loess.o: loess.cpp
	g++ -c  -I$(MATLAB_PATH)/extern/include -I$(MATLAB_PATH)/simulink/include -DMATLAB_MEX_FILE -ansi -D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -pthread -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -std=c++11 -frounding-math -O3 -DNDEBUG -pedantic -Wall  -O -DNDEBUG  "loess.cpp"	

loess.mexa64: loess.o
	g++ -O -pthread -shared -Wl,--version-script,$(MATLAB_PATH)/extern/lib/glnxa64/mexFunction.map -Wl,-Bsymbolic-functions -Wl,-z,relro -Wl,--no-undefined -frounding-math -o  "loess.mexa64" loess.o  -Wl,-rpath-link,$(MATLAB_PATH)/bin/glnxa64 -L$(MATLAB_PATH)/bin/glnxa64 -lmx -lmex -lmat -lm -rdynamic -lCGAL_Core -lCGAL -lgmpxx -lmpfr -lgmp -lboost_thread -lboost_system -std=c++0x
.PHONY : clean
clean:
	rm -rf *.mex* *.o
