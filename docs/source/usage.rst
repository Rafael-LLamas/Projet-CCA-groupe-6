Usage
=====

This project provides an implementation of operations for Toeplitz-like matrices in C, utilizing the **FLINT** (Fast Library for Number Theory) mathematical library.

Building the Project
--------------------

You can compile the project using CMake by creating a build directory and running the following commands:

.. code-block:: console

   $ cmake -B build -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
   $ cd build
   $ make

Benchmarking
------------

You can benchmark the project using the main executable. Run the executable without any arguments to view the usage instructions:

.. code-block:: console

   $ ./main <command> <options>

Dependencies
------------

While FLINT with CMake is primarily geared toward Windows users, this project utilizes the GNU autotools build system for Unix-like environments. 

Prerequisites Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~

Install the required package management dependencies for your specific operating system:

**macOS (Homebrew):**

.. code-block:: console

   $ brew install autoconf automake libtool gmp mpfr

**Fedora:**

.. code-block:: console

   $ sudo dnf install autoconf automake libtool gmp-devel mpfr-devel

**Arch Linux:**

.. code-block:: console

   $ sudo pacman -S autoconf automake libtool gmp mpfr

**Debian/Ubuntu:**

.. code-block:: console

   $ sudo apt install autoconf automake libtool libgmp-dev libmpfr-dev

Building FLINT From Source
~~~~~~~~~~~~~~~~~~~~~~~~~~

Navigate to your FLINT directory and initialize the bootstrap script:

.. code-block:: console

   $ cd flint-3.4.0/
   $ ./bootstrap.sh

Configure the build according to your system architecture. Ensure you replace the paths to match your local specifications.

**Example for Intel-based macOS (using Homebrew and native architecture optimization):**

.. code-block:: console

   $ ./configure \
        --disable-assert \
        --enable-avx2 \
        --with-gmp-include=/usr/local/opt/gmp/include \
        --with-gmp-lib=/usr/local/opt/gmp/lib \
        --with-mpfr=/usr/local/opt/mpfr \
        CC=clang \
        CFLAGS="-Wall -O3 -march=native"

**Example for Apple Silicon macOS (M1/M2/M3 using Homebrew):**

.. code-block:: console

   $ ./configure \
        --disable-assert \
        --with-gmp-include=/opt/homebrew/opt/gmp/include \
        --with-gmp-lib=/opt/homebrew/opt/gmp/lib \
        --with-mpfr=/opt/homebrew/opt/mpfr \
        CC=clang \
        CFLAGS="-Wall -O3 -march=native"

Compile, test, and install the library based on your operating system:

**Linux:**

.. code-block:: console

   $ make -j$(nproc)
   $ make -j$(nproc) check
   $ sudo make install

**macOS:**

.. code-block:: console

   $ make -j$(sysctl -n hw.ncpu)
   $ make -j$(sysctl -n hw.ncpu) check
   $ sudo make install