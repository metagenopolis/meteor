
########################

OS := $(shell uname)
ifeq ($(OS), Darwin)
	CXX = clang++
	LINKER = $(CXX)
else
	CXX = g++
	LINKER = $(CXX)
endif

CXX_FLAGS = -O3 -Wall -std=c++98 

# Final binaries
BIN1 = meteor-counter
BIN2 = meteor-profiler
OUT ?= ./build
# Put all auto generated stuff to this build dir.
BUILD_DIR = $(OUT)
# List of all .cpp source files.
CPP := $(wildcard *.cpp)
CPP := $(filter-out counter_main.cpp, $(CPP))
CPP := $(filter-out profiler_main.cpp, $(CPP))
CPP1 = counter_main.cpp $(CPP)
CPP2 = profiler_main.cpp $(CPP)

# All .o files go to build dir.
OBJ1 = $(CPP1:%.cpp=$(BUILD_DIR)/%.o)
OBJ2 = $(CPP2:%.cpp=$(BUILD_DIR)/%.o)
# The compiler will create these .d files containing dependencies.
DEP1 = $(OBJ1:%.o=%.d)
DEP2 = $(OBJ2:%.o=%.d)

# Build meteor-counter and meteor-profiler
all : $(BUILD_DIR)/$(BIN1) $(BUILD_DIR)/$(BIN2)

# Build meteor-counter 
$(BIN1) : $(BUILD_DIR)/$(BIN1)
# Build meteor-profiler
$(BIN2) : $(BUILD_DIR)/$(BIN2)

# Actual target of the binaries 
$(BUILD_DIR)/$(BIN1) : $(OBJ1)
# 	@# Create build directories - same structure as sources.
# 	@mkdir -p $(@D)
	@# Just link all the object files.
	$(LINKER) $(CXX_FLAGS) $^ -o $@

$(BUILD_DIR)/$(BIN2) : $(OBJ2)
# 	@# Create build directories - same structure as sources.
# 	@mkdir -p $(@D)
	@# Just link all the object files.
	$(LINKER) $(CXX_FLAGS) $^ -o $@

# Include all .d files
-include $(DEP1)
-include $(DEP2)

# Build target for every single object file.
# The potential dependency on header files is covered
# by calling `-include $(DEP)`.
$(BUILD_DIR)/%.o : %.cpp
	mkdir -p $(@D)
	@# The -MMD flags additionaly creates a .d file with
	@# the same name as the .o file.
	$(CXX) $(CXX_FLAGS) -MMD -c $< -o $@

.PHONY : clean
clean :
	@# Remove all generated files.
	rm -f $(BUILD_DIR)/$(BIN1) $(OBJ1) $(DEP1) $(BUILD_DIR)/$(BIN2) $(OBJ2) $(DEP2)
	
#############

