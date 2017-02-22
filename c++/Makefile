#
#
#
BOOST_INCLUDE = /home/fzeng/anaconda/include
BOOST_LIBRARY = /home/fzeng/anaconda/lib


#
#
#CC       = /home/fzeng/Tool/gnu/bin/g++
CC       = g++
CFLAGS   = -g -Wall -O3 -std=c++11 -fopenmp
INCLUDES = -I $(BOOST_INCLUDE)
DFLAGS   = 
LDFLAGS  = -L $(BOOST_LIBRARY)
LIBS     = -lz \
           -lboost_serialization \
           -lboost_iostreams \
           -lboost_program_options \
           -lboost_filesystem \
           -lboost_system


#
#
SRCS = InvertedIndex.cpp LSH.cpp SparseMatrix.cpp SparseVector.cpp \
       fesco_inverted_index_build.cpp \
       fesco_inverted_index_stats.cpp \
       fesco_inverted_index_sketch.cpp \
       fesco_inverted_index_cossim.cpp \
       Fesco.cpp
DEPS = kseq.h 
OBJS = $(SRCS:.cpp=.o)      


all: fesco

#
#
.PHONY: depend clean


fesco: $(OBJS) $(DEPS)
	$(CC) $(CFLAGS) $(INCLUDES) $(DFLAGS) $(LDFLAGS) $(LIBS) $(OBJS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) *.o *~ fesco

depend: $(SRCS) 
	makedepend $(INCLUDES) $^
