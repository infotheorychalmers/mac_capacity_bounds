Detailed installation for MAC
=============================

This document provides one way to get our code running on a MAC OS X. 

 * Install MKL ( or any other BLAS LAPACK library )
   * Download Intel Parallel Studio for MAC OS from [here](https://software.intel.com/en-us/intel-parallel-studio-xe).
   * Double click on the dmg file and follow the instructions to install it.
   * Add the following line in the .bash_profile file and restart the terminal.
```
   export DYLD_LIBRARY_PATH=/opt/intel/lib:/opt/intel/mkl/lib:$DYLD_LIBRARY_PATH
```



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

Our code should now work.
