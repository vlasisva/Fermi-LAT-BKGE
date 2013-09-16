########                    Definitions

OPT = -O3 -fPIC

ROOTLIB    = $(shell root-config --libs)
ROOTFLAGS  = $(shell root-config --cflags)
ROOTINC    = $(shell root-config --incdir)

TARGET_OBJS = $(patsubst %.cxx,%.o,$(wildcard src/BackgroundEstimator/*.cxx)) $(patsubst %.cxx,%.o,$(wildcard src/BKGE_Tools/*.cxx)) 
TARGET_OBJS_r = $(patsubst %.cxx,%.o,$(wildcard src/rootIrfLoader/*.cxx))

EXTLIBS    = -L${FERMI_DIR}/lib -L..
INCLUDES   = -I. -I${FERMI_DIR}/include -Iinclude/rootIrfLoader -Iinclude -I$(ROOTINC)


LINK_TARGET   = g++ $(OPT) -shared -o ../libBKGE.so  $(TARGET_OBJS)  $(TARGET_OBJS_r) $(EXTLIBS) $(ROOTLIB) -lcfitsio_3.29 -lastro -lCLHEP-2.1.0.1 -lfacilities  -lCore -lRIO  


########                     Targets
all: rootIrfLoader BKGE

clean:
	@rm `find . -iname "*.o"`
	@rm ../libBKGE.so

BKGE:	$(TARGET_OBJS)
	rootcint -l -f src/BackgroundEstimator/Dict.cxx -c -p $(INCLUDES) include/*.h
	g++ $(OPT) $(ROOTFLAGS) -c src/BackgroundEstimator/Dict.cxx $(INCLUDES) -o src/BackgroundEstimator/Dict.o
	$(LINK_TARGET)

rootIrfLoader:	$(TARGET_OBJS_r)



########                 Default Rules
.SUFFIXES: .cxx .o

.cxx.o:;g++ $(OPT) $(ROOTFLAGS) -c $< $(INCLUDES) -o $@

src/BackgroundEstimator/Dict.o:
	
