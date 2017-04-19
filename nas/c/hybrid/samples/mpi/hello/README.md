testing pgprof
compile: mpicc/pgcc -Minfo=ccff
prof: pgprof -o a.prof ./(file).out
      pgprof -i a.prof
