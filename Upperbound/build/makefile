

ROOT = ..
SRC = $(ROOT)/src
INC = $(ROOT)/include

# ITPP values
CFLAGS += `pkg-config --cflags itpp`
LDFLAGS += `pkg-config --libs itpp` -lpthread 

# GSL values
CFLAGS += `pkg-config --cflags gsl`
LDFLAGS +=  `pkg-config --libs gsl`

#BOOST values
#CFLAGS += -I<PATH/TO/BOOST/include>
#LDFLAGS += -L<PATH/TO/BOOST/lib>
LDFLAGS += -lboost_thread -lboost_program_options -lboost_system -lboost_date_time

#Code values
CFLAGS += -I$(INC)  

all : run

run : main.o UpperBound.o
	g++ main.o UpperBound.o -o run $(LDFLAGS)

UpperBound.o : $(SRC)/UpperBound.cpp $(INC)/UpperBound.h
	g++ $(CFLAGS) -c $(SRC)/UpperBound.cpp 

main.o : $(SRC)/main.cpp $(INC)/UpperBound.h
	g++ $(CFLAGS) -c $(SRC)/main.cpp  

clean :
	rm -rf *.o run debug 
