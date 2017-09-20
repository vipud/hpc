This program uses CUDA to calculate the Fibonacci sequence on a GPU.
To compile:

nvcc fib_cuda.cu -o fib_cuda

To run:

./fib_cuda x

where x is the desired numbers of the sequnce you wish to calculate.


Note: this program uses unsigned long integers, which overflow after the 93rd number in the sequence.



Analysis
Here is the output from nvprof of a kernel with only 1 block and 1 thread:

Time(%)      Time     Calls       Avg       Min       Max  Name
 84.98%  7.4240us         1  7.4240us  7.4240us  7.4240us  Fibonacci(double*, double*, int, double, double, double)
  9.52%     832ns         1     832ns     832ns     832ns  [CUDA memcpy HtoD]
  5.49%     480ns         1     480ns     480ns     480ns  [CUDA memcpy DtoH]

Here is the output from our modified kernel:

Time(%)      Time     Calls       Avg       Min       Max  Name
 82.16%  5.6000us         5  1.1200us     992ns  1.5040us  Fibonacci(unsigned long*, int)
 13.15%     896ns         1     896ns     896ns     896ns  [CUDA memcpy HtoD]
  4.69%     320ns         1     320ns     320ns     320ns  [CUDA memcpy DtoH]

Parallelizing the Fibonacci sequence shaved off roughly 2 microseconds, when calculating only the first 93 numbers.
