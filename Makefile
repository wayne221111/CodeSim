CPPFLAGS = -O3 -fopenmp

all: CodeSim randomc cmaes

CMAES.out: all CMAES_main.cpp
	$(CXX) $(CPPFLAGS) CMAES_main.cpp mersenne.o userintf.o cmaes.o Bit.o -o $@

CMAES2.out: all CMAES2.cpp
	$(CXX) $(CPPFLAGS) CMAES2.cpp mersenne.o userintf.o cmaes.o Bit.o -o $@

CMAES3.out: all CMAES3.cpp laguer.o mrqmin.o mrqcof.o gaussj.o covsrt.o histogram.out
	$(CXX) $(CPPFLAGS) CMAES3.cpp laguer.o mrqmin.o mrqcof.o gaussj.o covsrt.o mersenne.o userintf.o cmaes.o Bit.o -o $@

CMAES4.out: all CMAES4.cpp histogram.out
	$(CXX) $(CPPFLAGS) CMAES4.cpp mersenne.o userintf.o cmaes.o Bit.o -o $@

CMAES5.out: all CMAES5.cpp histogram.out
	$(CXX) $(CPPFLAGS) CMAES5.cpp mersenne.o userintf.o cmaes.o Bit.o -o $@

histogram.out: all LT_sim_histogram.cpp
	$(CXX) $(CPPFLAGS) LT_sim_histogram.cpp mersenne.o userintf.o Bit.o -o $@

exp_SVC.out : randomc CodeSim exp_SVC.cpp
	$(CXX) $(CPPFLAGS) exp_SVC.cpp mersenne.o userintf.o ConvoCode.o Bit.o -o $@

exp_LT_sim.out : randomc CodeSim exp_LT_sim.cpp
	$(CXX) $(CPPFLAGS) exp_LT_sim.cpp mersenne.o userintf.o ConvoCode.o Bit.o -o $@

decoding_pattern_test.out: CodeSim decoding_pattern_test.cpp
	$(CXX) $(CPPFLAGS) decoding_pattern_test.cpp mersenne.o userintf.o Bit.o -o $@

fitting_data_set_generator.out: fitting_data_set_generator.cpp
	$(CXX) $(CPPFLAGS) fitting_data_set_generator.cpp -o $@

LT_BER.out: LT_BER.cpp CodeSim randomc
	$(CXX) $(CPPFLAGS) LT_BER.cpp mersenne.o userintf.o Bit.o -o $@

LT_BER2.out: LT_BER2.cpp CodeSim randomc
	$(CXX) $(CPPFLAGS) LT_BER2.cpp mersenne.o userintf.o Bit.o -o $@

ConvoCodeSim.out: ConvoCodeSim.cpp CodeSim randomc
	$(CXX) $(CPPFLAGS) ConvoCodeSim.cpp mersenne.o userintf.o Bit.o ConvoCode.o -o $@

ConvoCodeSim_byte.out: ConvoCodeSim_byte.cpp CodeSim randomc
	$(CXX) $(CPPFLAGS) ConvoCodeSim_byte.cpp mersenne.o userintf.o Bit.o Byte.o ConvoCode.o -o $@

ConvoLTSim.out: ConvoCodeSim.cpp CodeSim randomc
	$(CXX) $(CPPFLAGS) ConvoLTSim.cpp mersenne.o userintf.o Byte.o Bit.o ConvoCode.o -o $@

SVC_UEP_exp_byte.out: ConvoCodeSim.cpp CodeSim randomc
	$(CXX) $(CPPFLAGS) SVC_UEP_exp_byte.cpp mersenne.o userintf.o Bit.o Byte.o ConvoCode.o -o $@

histogram_cumulator.out: histogram_cumulator.cpp
	$(CXX) -O3 histogram_cumulator.cpp -o $@

cmaes: cmaes.o cmaes.h cmaes_interface.h
	
randomc: randomc.h mersenne.o userintf.o

CodeSim: CodeSim.h Bit.o ConvoCode.o LT.h LTCode.h Byte.o

commitall:
	git commit -a

push:
	git push github master
pull:
	git pull github master

.PHONY: clean commitall push pull
clean: 
	@rm *.o *.out
