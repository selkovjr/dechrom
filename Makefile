GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)

OBJ = dechrom.o survey.o plot.o util.o cfa_mask.o write_tiff.o progressbar.o

BIN = dechrom
LIBBIN=.

CXX=g++

HDRDIR=-I/opt/local/include/ -I/usr/local/include/
LIBDIR=-L/home/selkovjr/lib -L/opt/local/lib/ -L/usr/local/lib/


CFLAGS += -std=c++11 -g -O3 -D VERSION=\"$(GIT_VERSION)\" -fdiagnostics-color=auto \
  -fopenmp -funroll-loops -fomit-frame-pointer  -fno-tree-pre -falign-loops -ffast-math -ftree-vectorize \
  -Weffc++ -pedantic -Wall -Wextra  -Wno-write-strings -Wno-deprecated  $(HDRDIR) \
  `Magick++-config --cxxflags`

LDFLAGS += -g $(CFLAGS) $(LIBDIR) -lncurses -lgomp -lpthread \
  `pkg-config --static --libs opencv` \
  -lMagick++-7.Q16HDRI -lMagickWand-7.Q16HDRI -lMagickCore-7.Q16HDRI


LIBMX= survey.o find.o radial.o plot.o util.o cfa_mask.o write_tiff.o progressbar.o

default: $(OBJ) $(BIN)

$(OBJ) : %.o : %.cpp survey_args.h
	$(CXX) -c $(CFLAGS)  $< -o $@

$(BIN) : % : %.o  $(LIBMX)
	$(CXX) -o $(LIBBIN)/$@  $^ $(LDFLAGS)

radial.o: radial.cpp radial_args.h
	g++ -c -std=c++11 -g -O3 \
    -fopenmp -funroll-loops -fomit-frame-pointer -ffast-math \
    -DMAGICKCORE_HDRI_ENABLE=1 -DMAGICKCORE_QUANTUM_DEPTH=16 \
    -I/home/selkovjr/include/ImageMagick-7 \
    $<

find.o: find.cpp find_args.h
	g++ -c -std=c++11 -g -O3 -fopenmp -funroll-loops -fomit-frame-pointer -ffast-math -Wall -Wextra $<


.PHONY : clean
clean:
	$(RM) $(OBJ) radial.o find.o; cd $(LIBBIN); rm -f $(BIN)

