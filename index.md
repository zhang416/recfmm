recFMM is a program representation and implementation of a recursive scheme for parallelizing 
the adaptive fast multipole method (FMM) on shared-memory computers. It achieves remarkable 
high performance while maintaining mathematical clarity and flexibility. 
The parallelization scheme signifies the recursion feature that is intrinsic to the FMM but 
was not well exploited. The program modules of RECFMM constitute a map between numerical 
computation components and advanced architecture mechanisms. The mathematical structure is 
preserved and exploited, not obscured nor compromised, by parallel rendition of the recursion 
scheme. Modern software system—CILK in particular, which provides graph-theoretic optimal 
scheduling in adaptation to the dynamics in parallel execution—is employed. 
recFMM supports multiple algorithm variants that mark the major advances with low-frequency 
interaction kernels, and includes the asymmetrical version where the source particle ensemble 
is not necessarily the same as the target particle ensemble. 


### Reference
* B. Zhang, J. Huang, N. P. Pitsianis, X. Sun. **recFMM: Recursive Parallelization of the Adaptive
Fast Multipole Method for Coulomb and Screened Coulomb Interactions.** _Commun. Comput. Phys._ 
20 (2016) 534

### Acknoledgments
recFMM was developed under the support of NSF grants CCF-0905164, CCF-0905473, and ACI-1440396. 
