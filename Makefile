#
# compiler
#
GRBROOT = /Library/gurobi1000
PLATFORM = macos_universal2
TWOUP    = $(GRBROOT)/$(PLATFORM)
INC      = $(TWOUP)/include/
CPP      = g++
CARGS    = -m64 -g
CPPLIB   = -L$(TWOUP)/lib -lgurobi_c++ -lgurobi100

#
GRBAPP   = dotnet
DOTNETFRAMEWORK = --framework=netcoreapp6.0
CC       	  = clang++
CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG

#
# scots 
#
SCOTSROOT		= ../../..
SCOTSINC		= -I$(SCOTSROOT)/bdd -I$(SCOTSROOT)/srcc -I$(SCOTSROOT)/utils 

#
# cudd 
#
CUDDPATH		=  $(SCOTSROOT)/cudd-3.0.0
CUDDINC 		= -I$(CUDDPATH)
CUDDLIBS		= -lcudd 
CUDDLPATH   	= -L$(CUDDPATH)/lib

# TARGET = vehicle

# all2: vehicle.cc
# 	$(CC) $(CARGS) $(CXXFLAGS) -o $@ $< -I$(INC) $(CUDDINC) $(SCOTSINC) $(CPPLIB) -lm $(CUDDLPATH) $(CUDDLIBS)

# $(TARGET): $(TARGET).o
# 	$(CC) $(CXXFLAGS) -o $(TARGET) $(TARGET).o $(CUDDLPATH) $(CUDDLIBS)

TARGET2 = init_abstr

ia: init_abstr.cc
	$(CC) $(CARGS) $(CXXFLAGS) -o $@ $< -I$(INC) $(CUDDINC) $(SCOTSINC) $(CPPLIB) -lm $(CUDDLPATH) $(CUDDLIBS)

$(TARGET2): $(TARGET2).o
	$(CC) $(CXXFLAGS) -o $(TARGET2) $(TARGET2).o $(CUDDLPATH) $(CUDDLIBS)


# TARGET1 = vehicleModel

# %.o:%.cc
# 	$(CC) -c $(CXXFLAGS) $(CUDDINC) $(SCOTSINC) $< -o $@

# $(TARGET1): $(TARGET1).o
# 	$(CC) $(CXXFLAGS) -o $(TARGET1) $(TARGET1).o $(CUDDLPATH) $(CUDDLIBS)


clean:
	rm -rf *.o *_c *_c++ *.class *.log *.rlp *.lp *.bas *.ilp *.mps *.prm *.dSYM; \
	if [ -d $(GRBAPP) ]; then \
		cd $(GRBAPP); \
		find . ! \( -name "gurobi*.csproj" -o -name . \) -exec rm -rf {} +; \
	fi


clean:
	rm  ./$(TARGET1)  ./$(TARGET1).o

