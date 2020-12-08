INCLUDE	:=	-I /usr/include/eigen3
CPP	:=	g++
LINK	:=	-lgtest -lpthread

#exercise:
#	$(CPP)  exercise.cpp $(LINK) $(INCLUDE)

algebra: algebra.hpp
	$(CPP)  algebra_test.cpp $(LINK) $(INCLUDE)
all:	algebra
