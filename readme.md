# Project Title

Project Descrption here

## Build
You can create a build folder and compile the project using CMake.
```bash
cmake -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cd build
make
```

## Benchmark
It is possible to benchmark the proejct using the main executable.
```bash
./main benchmark
```


## Dependencies for a Mac user:
> Flint with CMake is primarily for Windows users, I use the GNU autotools build. Install prerequisites:
> 
> `brew install autoconf automake libtool gmp mpfr`
>
> To build FLINT:
>
> ```bash
> cd flint-3.4.0/
>
> ./bootstrap.sh
>
> # Replace path correctly to your specification
> # March native as my CPU is an older generation
> ./configure \
>    --disable-assert \
>    --enable-avx2 \
>    --with-gmp-include=/usr/local/opt/gmp/include \
>    --with-gmp-lib=/usr/local/opt/gmp/lib \
>    --with-mpfr=/usr/local/opt/mpfr \
>    CC=clang \
>    CFLAGS="-Wall -O3 -march=native"
>
> make -j9          # I only have 4 cores that run 2 threads (+1)
> make -j 9 check   # Verify hardware <-> functions (8 threads +1)
> make install
> ```
> 

## Autors
- Ali Arda Barut
- Rafael LLamas
