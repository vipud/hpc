test: ep.c c_timers.c print_results.c randdp.c wtime.c
	pgcc -o ep -acc -ta=tesla:cc60 -Minfo=accel ep.c c_timers.c print_results.c wtime.c randdp.c -lm

orig: ep_orig.c c_timers.c print_results.c randdp.c wtime.c
	pgcc -o orig -acc -ta=multicore -Minfo=all ep.c c_timers.c print_results.c wtime.c randdp.c -lm
