Detailed installation for MAC
=============================

This document explains how to get the code to compile and run in MAC OS. First we have to install the dependencies, the steps for the same are as follows:

 * Install MKL ( or any other BLAS LAPACK library )
   * Download Intel Paralell Studio for MAC OS from [here](https://software.intel.com/en-us/intel-parallel-studio-xe).
   * Double click on the dmg file and follow the instructions to install it.
   * Add the following statement inside ~/.bash_profile file
```
   export DYLD_LIBRARY_PATH=/opt/intel/lib:/opt/intel/mkl/lib:$DYLD_LIBRARY_PATH
```
   * Close and reopen terminal.



* Install pkg-config
   * Download the latest version from [here](http://pkgconfig.freedesktop.org/releases/pkg-config-0.28.tar.gz).
   * Unzip the downloaded file.
   * Move to the unzipped directory.
   * Run 
```
    ./configure ( or './configure --with-internal-glib' if you encounter errors with glib)
    make
    sudo make install
```

 * Install cmake
   * Download the latest source from [here](http://www.cmake.org/download/).
   * Unzip the downloaded file.
   * Move to the unzipped directory.
   * Run 
```
    ./bootstrap
    make
    sudo make install
```

 * Install IT++
   * Download the latest source from [here](http://sourceforge.net/projects/itpp/files/).
   * Unzip the downloaded file.
   * Move to the unzipped directory.
   * Run
```
    mkdir build
    cd build
    cmake .. -DBLA_VENDOR=Intel11 ( change accordingly for different vendor )
    make
    sudo make install
```
 * Install GSL
   * Download the latest source from [here](http://www.gnu.org/software/gsl/).
   * Unzip the downloaded file.
   * Move to the unzipped directory.
   * Run
```
    ./configure
    make
    sudo make install
```
 * Install boost
   * Download latest source from [here](http://www.boost.org/doc/libs/1_57_0/more/getting_started/unix-variants.html).
   * Unzip the downloaded file.
   * Move to the unzipped directory.
   * Run
```
    ./bootstrap.sh --prefix=/usr/local
    sudo ./b2 install
```

Now we have completed installation of dependencies. Now inorder to compile and run the code to compute the Gaussian lower bound for 5 users, coherence interval 20, and SNR values from 0 to 10, run

```
    cd Gaussian/build
    make
    ./run --K 5 --T 20 --MCX 20000 --startSNR 0 --deltaSNR 1 --endSNR 10
```
