 CC	=	g++ -Wall -std=c++11 -O3 -pg
UNAME = $(shell uname)
BINDIR = bin
OBJDIR = .obj
SRCDIR = src
DIRS = $(BINDIR) $(OBJDIR)
SOURCESXX := $(wildcard $(SRCDIR)/*.cxx)
SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
EXECUTABLES := $(SOURCESXX:$(SRCDIR)/%.cxx=$(BINDIR)/%)

ifeq ($(UNAME), Darwin)
   LIBS = -L/opt/local/lib -L/usr/lib -lgsl -lblas -I/opt/local/include
else
   LIBS	=  -L/usr/lib -lgsl -lblas -l armadillo -ljsoncpp -lopenblas
   INCS = -I /usr/local/include -I/usr/include -I/opt/local/include -I/usr/include/jsoncpp
endif

all: $(DIRS) $(EXECUTABLES)

directories : $(DIRS)

$(DIRS):
	mkdir -p $(DIRS)

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS) $(INCS)
$(EXECUTABLES): $(BINDIR)/% : $(SRCDIR)/%.cxx $(OBJECTS)
	$(CC) -o $@ $< $(OBJECTS) $(LIBS) $(INCS)

.PHONY: all clean directories

clean:
	rm $(EXECUTABLES)
	rm $(OBJECTS)
