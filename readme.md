# Toeplitz-like matrices implementation

In this project, we show a way to implement operation for Toeplitz-like matrices in C with the lib FLINT.

## Build
You can create a build folder and compile the project using CMake.
```bash
cmake -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cd build
make
```

## Benchmark
It is possible to benchmark the project using the main executable.
```bash
./main <command> <options>
```
You can call the executable with no argument to get the usage.


## Dependencies for a non Windows user:
Flint with CMake is primarily for Windows users, I use the GNU autotools build. Install prerequisites:
### Mac:
 ```bash
 brew install autoconf automake libtool gmp mpf
 ```
### Fedora :
```bash
 sudo dnf install autoconf automake libtool gmp-devel mpfr-devel
 ```
 ### Arch :
 ```bash
 sudo pacman -S autoconf automake libtool gmp mpfr
 ```
 ### Debian-like :
  ```bash
 sudo apt install autoconf automake libtool libgmp-dev libmpfr-dev
 ```
 To build FLINT:
 ```bash
 cd flint-3.4.0/

 ./bootstrap.sh
 ```
 Replace path correctly to your specification
 This is an example on an intel based Mac using Homebrew (March native as my CPU is an older generation)
 ```bash
 ./configure \
    --disable-assert \
    --enable-avx2 \
    --with-gmp-include=/usr/local/opt/gmp/include \
    --with-gmp-lib=/usr/local/opt/gmp/lib \
    --with-mpfr=/usr/local/opt/mpfr \
    CC=clang \
    CFLAGS="-Wall -O3 -march=native"
   ```
 This is an example on an Apple Silicon based Mac using Homebrew
 ```bash
 ./configure \
   --disable-assert \
   --with-gmp-include=/opt/homebrew/opt/gmp/include \
   --with-gmp-lib=/opt/homebrew/opt/gmp/lib \
   --with-mpfr=/opt/homebrew/opt/mpfr \
   CC=clang \
   CFLAGS="-Wall -O3 -march=native"
   ```
 Make the FLINT library and install it:
 ### Linux:
 ```bash
 make -j$(nproc)
make -j$(nproc) check
sudo make install
 ```
 ### Mac:
 ```bash
 make -j$(sysctl -n hw.ncpu)
make -j$(sysctl -n hw.ncpu) check
sudo make install       
 ```

## Dependencies for a non Windows user:
Good Luck

## Autors
- Ali Arda Barut
- Rafael LLamas
