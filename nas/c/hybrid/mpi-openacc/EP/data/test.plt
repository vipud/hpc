# data file
set xtics ("Sequential" 0.5, "1080" 2.5, "K40" 4.5, "K80" 6.5, "1080+K80" 8.5,)
set boxwidth 0.2
set style fill solid
plot "perf.dat" every 2 using 1:2 with boxes ls 1,\
"perf.dat" every 2::1 using 1:2 with boxes ls 2
# 1 - Sequential
# 2 - GTX 1080 (OpenACC)
# 3 - Tesla K40 (OpenACC)
# 4 - Tesla K80 (OpenACC + MPI)
# 5 - GTX 1080 + Tesla K80 (OpenACC + MPI)
