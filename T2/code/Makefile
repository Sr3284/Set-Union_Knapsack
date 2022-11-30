SCIP_LIB=/usr/local/lib
TRACE=DEBUG
CFLAGS=-D NO_CONFIG_HEADER -D SCIP_VERSION_MAJOR -D$(TRACE)
LDFLAGS=-L $(SCIP_LIB)

bin/mochila: bin/cmain.o bin/probdata_mochila.o bin/problem.o
	gcc $(CFLAGS) $(LDFLAGS) -o bin/mochila bin/cmain.o bin/probdata_mochila.o bin/problem.o -lscip -lm

bin/cmain.o: src/cmain.c
	gcc $(CFLAGS) -c -o bin/cmain.o src/cmain.c

bin/probdata_mochila.o: src/probdata_mochila.c src/probdata_mochila.h
	gcc $(CFLAGS) -c -o bin/probdata_mochila.o src/probdata_mochila.c

bin/problem.o: src/problem.c src/problem.h
	gcc $(CFLAGS) -c -o bin/problem.o src/problem.c

.PHONY: clean

clean:
	rm -f bin/*.o bin/mochila

