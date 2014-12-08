MAC capacity bounds
===================

Code for the upper bound and the lower bound on the MAC sum-rate capacity.

Dependencies
------------
The code is written in C/C++ using the following libraries 
 * GNU Scientific library (version 1.14)
 * IT++ (version 4.3.1)
 * boost (version 1.42)

Compiling
---------
Inorder to compile the code for each bound
 * cd to the specific directory.
 * cd inside to the 'build' directory and run
```
	make all
```
This will compile the code. All the output files (.o files) and executable named 'run' will be created in the same directory

Running
-------
The options for running each bound can be checked by using help option as
```
		run --help
```
For example if you have to get Gaussian lower bound for 5 user scenario, coherence interval 20 channel uses, and SNR values from 0 to 10 you run the following
```
		run --K 5 --T 20 --MCX 20000 --startSNR 0 --deltaSNR 1 --endSNR 10
```



	

