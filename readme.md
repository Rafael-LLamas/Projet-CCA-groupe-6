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
> ```bash
> brew install autoconf automake libtool gmp mpf
> ```
>
> To build FLINT:
>
> ```bash
> cd flint-3.4.0/
>
> ./bootstrap.sh
> ```
> Replace path correctly to your specification
> This is an example on an intel based Mac using Homebrew (March native as my CPU is an older generation)
> ```bash
> ./configure \
>    --disable-assert \
>    --enable-avx2 \
>    --with-gmp-include=/usr/local/opt/gmp/include \
>    --with-gmp-lib=/usr/local/opt/gmp/lib \
>    --with-mpfr=/usr/local/opt/mpfr \
>    CC=clang \
>    CFLAGS="-Wall -O3 -march=native"
>   ```
> This is an example on an Apple Silicon based Mac using Homebrew
> ```bash
> ./configure \
>   --disable-assert \
>   --with-gmp-include=/opt/homebrew/opt/gmp/include \
>   --with-gmp-lib=/opt/homebrew/opt/gmp/lib \
>   --with-mpfr=/opt/homebrew/opt/mpfr \
>   CC=clang \
>   CFLAGS="-Wall -O3 -march=native"
>   ```
> Make the FLINT library and install it:
> ```bash
> make -j11          # I only have 10 cores
> make -j 11 check   # Verify hardware <-> functions (10 threads +1)
> make install       # sudo if necessary
> ```


## Autors
- Ali Arda Barut
- Rafael LLamas
