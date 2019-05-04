# compiler:
CC = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall -O2 -fPIC -std=c++11
CFLAGS += -I /usr/local/include

# Options:
ifeq ($(ALLOWDIEOUT), 1)
	CFLAGS += -D_ALLOWDIEOUT_
endif

ifeq ($(DEBUG), 1)
	CFLAGS += -D_DEBUG_
endif

ifneq ($(NODISPLAY), 1)
	CFLAGS += -D_DISPLAY_
        CFLAGS_X11 += -I/opt/X11/include
	LFLAGS_X11 += -L/opt/X11/lib -lGL -lGLU -lGLUT
endif

# Classes:
SRC_DIR = src/
CLASSES = Phylogeny Universe Cell CellType Shapes SimulationParameterSet HSL2RGB CellSample
CLASS_FILES_BASE = $(addprefix $(SRC_DIR), $(CLASSES))
CLASS_FILES_COMP = $(addsuffix .o, $(CLASS_FILES_BASE))
CLASS_FILES_HEADER = $(addsuffix .hpp, $(CLASS_FILES_BASE))
CLASS_FILES_HEADER += $($(SRC_DIR), extern_global_variables.hpp)
CLASS_FILES_SOURCE = $(addsuffix .cpp, $(CLASS_FILES_BASE))

# Related to R package:
R_PKG_NAME = CHESS
R_PKG_VERS = 0.0.0.9000
R_PKG_DIR = r_package
R_PKG_SRC_DIR = $(R_PKG_DIR)/$(SRC_DIR)
R_PKG_SRC_FILES = $(CLASS_FILES_HEADER) $(CLASS_FILES_SOURCE) $(SRC_DIR)/extern_global_variables.hpp
#R_package_SSBP.cpp extern_global_variables.hpp
R_PKG_SRC_FILES_PATH = $(addprefix $(R_PKG_DIR)/, $(R_PKG_SRC_FILES))

# Define build targets:
TARGET = cancer_gillespie_simulation

OTHER_TARGETS += $(R_PKG_NAME)_$(R_PKG_VERS).tar.gz

# Make all targets:
all: $(TARGET) $(OTHER_TARGETS)

$(TARGET): src/main.cpp $(CLASS_FILES_COMP) $(CLASS_FILES_HEADER)
	$(CC) $(CFLAGS) $(CFLAGS_X11) $(LFLAGS) $(LFLAGS_X11) $(CLASS_FILES_COMP) src/main.cpp -o $@

%.o: %.cpp %.hpp
	$(CC) $(CFLAGS) $(CFLAGS_X11) -c $< -o $@

$(R_PKG_DIR)/% : %
	cp $< $@

$(R_PKG_NAME)_$(R_PKG_VERS).tar.gz: $(R_PKG_SRC_FILES_PATH)
	R $(RARGS) CMD build $(R_PKG_DIR) --no-resave-data --no-manual # source package

clean:
	$(RM) $(TARGET) $(CLASS_FILES_COMP) $(R_PKG_SRC_FILES_PATH) $(OTHER_TARGETS)
