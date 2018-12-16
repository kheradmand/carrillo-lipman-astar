Carrillo-Lipman A*
===================


Build and Use
-------------

Prerequisites:
- C++ compiler with support for C++11

To build the project, run:
```
make
```

This will make create the `align` executable, which can be used as follows:
```
usage: ./align [OPTIONS]
	 OPTIONS:
		-h,--help		print this help message
		-a,--aligner <name>	set aligner type, one of
			astar: A* with Carrillo-Lipman (CL) heuristic
			cl-cs: original CL method, z from center star alignment
			dp: naive dynamic programming. default: astar
		-s,--score <file>	uses the provided scoring file.
			default: 1 for mismatch and gap, 0 otherwise

```


For example 
```
./align -a dp < test/in2
score: 1
abcd
a-cd

```

or 
```
./align --score score/blosum62.score <  evaluate/ACBP
starting to align
open set 178 closed set 90 sum 268
space size 7830 visited 90
score: 943
SQAE-FDKAAEEVKHL-KTKPADEEM-LFIYSHYKQATVGDINTERPGMLDFKGKAKWDAWNELKGTSKEDAMKAYIDKVEELKKKYGI
HMAQVFEECVSFINGLPRTINLPNELKLDLYKYYKQSTIGNCNIKEPSAHKYIDRKKYEAWKSVENLNREDAQKRYVDIVSEIFPYWQD
```