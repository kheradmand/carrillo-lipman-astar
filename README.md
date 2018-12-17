Carrillo-Lipman A*
===================

A central problem in multiple string comparison is the Multiple SequenceAlignment (MSA) problem.
An easy to work with and widely used scoring function for measuring the goodness of a multiple alignment is the SP-score.
SP-score of a multiple alignment is obtained from the sum of the scores of the pairwise global alignments induced by the multiple alignment

The MSA problem with SP-score (MSA-SP) is known to be NP-complete [1]. 
The naive dynamic programming approach to this problem has the running time of O(2^kk^2n^k) for k strings of length. 
The problem can also be viewed as a graph search problem in which each cell in the dynamic programming table is viewed as a node in a graph and the neighboring cells are connected through edges with weights that correspond to the SP-scores. 
The goal is to find the minimum weight path between the first node and the node corresponding to the alignment of all sequences.   

Carillo and Lipman[2] show an interesting observation that can be used to provide a proper heuristic function for an A* to solve the MSA-SP problem.
The observation is as follows. Let v_1, v_2, v_3 be the sequences of length n to be aligned.
Let D^+(i,j,k) be the minimum SP-cost of alignment of suffixes v_1[i..n], v_2[j..n], and $_3[k..n].
Let D^+_{a,b}(x,y) be the minimum cost of alignment of v_a[x..n] and v_b[y..n].
Now we define $(i,j,k) = D^+_{1,2}(i,j) + D^+_{2,3}(j,k) + D^+_{3,1}(k,i) as the heuristic used for A*.

Note that originally, Carillo and Lipman used this observation in a different optimization that uses a non-optimal path z to prune the dynamic programming.
z is obtained using some efficient approximate or heuristic approach. 

Here compare the performance of all three approaches: the naive dynamic programming approach, the original Carillo-Lipman speedup, and implementation of A* using Carillo-Lipman as the heuristic function. 
For the original Carillo-Lipman speedup implementation, z is obtained using the center star alignment method [3] which is a very efficient 2-approximation of the optimal answer.   

Note: The implementation actually solves the weighted edit distance problem (a minimization problem with non-negative "costs"),
rather than the sequence alignment problem (a maximization problem with possibly negative "scores"), although we refer to cost and score interchangeability.
 

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
        -h,--help        print this help message
        -a,--aligner <name>    set aligner type, one of
            astar: A* with Carrillo-Lipman (CL) heuristic (default)
            cl-cs: original CL method, z from center star alignment
            dp: naive dynamic programming.
        -s,--score <file>    uses the provided scoring file.
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

Use as a library:
-----------------
All implementation are derived from the the `aligner_t` class:
```
class aligner_t {
protected:
    scoring_function_t& scoring;
    
public:
    aligner_t(scoring_function_t& _scoring) : scoring(_scoring) {}
    virtual std::pair<sequences_t,score_t> get_alignment(const sequences_t& seqs) = 0;
    virtual ~aligner_t(){}

};

class astar_aligner_t : public aligner_t;
class cl_star_aligner_t : public aligner_t;
class dp_aligner_t : public aligner_t;
```

`setquece_t` is a vector of strings (as input and output of alignment).
`scoring_function_t` provides scoring. It can also load a scoring function from a file (e.g. `score/blosum62.score`).



