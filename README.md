MAC capacity bounds
===================

Numerical routines for the computation of an upper bound and the lower bound on the MAC sum-rate capacity of a Rayleigh block-fading channel with no a priori CSI available at the transmitters and the receiver.

Dependencies
------------
The code is written in C/C++ using the following libraries 
 * GNU Scientific library (version 1.14)
 * IT++ (version 4.3.1)
 * boost (version 1.42)

Compiling
---------
In order to compile the code implementin each bound each bound
 * move to the specific directory.
 * move inside 'build' directory and run
```
	make all
```
This will compile the code. All the output files (.o files) and executable named 'run' will be created in the same directory

Running
-------
The options available for each bound can be consulted using the help option:
```
		run --help
```
For example, if you want to compute the Gaussian lower bound for 5 users, coherence interval 20, and SNR values from 0 to 10, run the following
```
		run --K 5 --T 20 --MCX 20000 --startSNR 0 --deltaSNR 1 --endSNR 10
```



	

