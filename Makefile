CXX = g++

EXEC = beda

EXECFINAL = ${EXEC}
SRC = src/

DIST = ${SRC}Distributions
KT = ${SRC}KinTran
LATT = ${SRC}Lattice
INIT = ${SRC}InitCond
KERN = ${SRC}Kernel
EL = ${SRC}Elastic
INELBASE = ${SRC}Inelastic
INEL = ${SRC}InelasticLPM
RUN = ${SRC}KTRun
KTDATA = ${SRC}KTData
MDATA = ${SRC}ManageData

OBJS = ${EXEC}.o ${RUN}.o ${KT}.o ${LATT}.o ${INIT}.o ${KERN}.o ${EL}.o  ${INEL}.o ${DIST}.o ${KTDATA}.o ${INELBASE}.o ${MDATA}.o

TKERN = testkernBH
KERNOBJS = ${TKERN}.o ${LATT}.o ${INIT}.o ${KERN}.o ${INEL}.o ${INELBASE}.o

KERNOBJSLPM = testLPM.o ${LATT}.o ${INIT}.o ${KERN}.o ${INEL}.o ${INELBASE}.o

${EXECFINAL}: ${OBJS}
	${CXX} -fopenmp -o $@ $^

${EXEC}.o: ${SRC}${EXEC}.cpp ${KT}.h ${INEL}.h ${RUN}.h
	${CXX} ${CFLAG} -fopenmp -c -o $@ $<

${MDATA}.o: ${MDATA}.cpp ${MDATA}.h
	${CXX} -c -o $@ $<

${LATT}.o: ${LATT}.cpp ${LATT}.h ${MDATA}.h
	${CXX} -c -o $@ $<

${INIT}.o: ${INIT}.cpp ${INIT}.h ${MDATA}.h
	${CXX} -c -o $@ $<

${KTDATA}.o: ${KTDATA}.cpp ${KTDATA}.h ${LATT}.h
	${CXX} -c -o $@ $<

${DIST}.o: ${DIST}.cpp ${DIST}.h ${KTDATA}.h ${LATT}.h ${INIT}.h
	${CXX} -c -o $@ $<

${KERN}.o: ${KERN}.cpp ${KERN}.h ${DIST}.h ${KTDATA}.h ${LATT}.h
	${CXX} -c -o $@ $<

${EL}.o: ${EL}.cpp ${EL}.h ${KERN}.h ${DIST}.h ${KTDATA}.h ${LATT}.h
	${CXX} -c -o $@ $<

${INELBASE}.o: ${INELBASE}.cpp ${INELBASE}.h ${KERN}.h ${DIST}.h ${KTDATA}.h ${LATT}.h
	${CXX} -fopenmp -c -o $@ $<

${INEL}.o: ${INEL}.cpp ${INEL}.h ${KERN}.h ${DIST}.h ${KTDATA}.h ${LATT}.h
	${CXX} -fopenmp -c -o $@ $<

${KT}.o: ${KT}.cpp ${KT}.h ${INEL}.h  ${EL}.h  ${KERN}.h ${KTDATA}.h
	${CXX} -fopenmp -c -o $@ $<

${RUN}.o: ${RUN}.cpp  ${KT}.h ${INEL}.h ${INELBASE}.h ${KERN}.h  ${EL}.h ${KERN}.h ${DIST}.h ${KTDATA}.h ${LATT}.h ${MDATA}.h
	${CXX} -fopenmp -c -o $@ $<


.PHONY: clean
clean:
	rm -f ${SRC}*.o *~ *# ${EXEC} test
