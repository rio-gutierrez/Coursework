SRCS = GF.c RK4.c funcs.c wave_periodic.c
OBJS = GF.o RK4.o funcs.o wave_periodic.o

PROG = wave_periodic
$(PROG): $(OBJS)
	gcc -Wall -g -fsanitize=address $(OBJS) -lm -o $@

%.o: %.c
	gcc -Wall -std=c99 -O2 -c $< 

clean:
	-rm $(OBJS) $(PROG)

cleanfull:
	-rm $(OBJS) $(PROG) *.asc

cleansuperfull:
	-rm $(OBJS) $(PROG) *.asc *.mp4 *.png
