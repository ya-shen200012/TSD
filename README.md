# TSD
For repeating numerical experiments in the paper titled  "Triangle Steepest Descent: A Geometry-Driven Gradient Method Variation with R-Linear Convergence"

After running the "TSD_20250819.jl" file, you will get a csv file "combined_results_20250819.csv" and a log file "error_log_20250819.txt". "combined_results_20250819.csv" contains the raw data for table A.6-A.10 in the paper.

## Implementation Instruction 
However, for Problem with index 96 named "LRCOVTYPE", BBQ cannot stop in 4h21m, TSD10 will stop after running for 4h, TSD50 cannnot stop for 4h. So we delete this problem in our practical running, we recommend readers to do so. Specifically, change the line 500-502 in "TSD_20250819.jl" as the following code, it will run the first 95 problems. This process takes about 3.76 hours.
``` julia
# julia
for filtered_problems_i in eachindex(filtered_problems[1:95])
    println(filtered_problems_i)
    postinx=0
```
change the line 500-502 in "TSD_20250819.jl" as the following code, it will run the remaining problems. This process takes about 4.98 hours.
``` julia
# julia
for filtered_problems_i in eachindex(filtered_problems[97:280])
    println(filtered_problems_i)
    postinx=96
```
## Summary of the failed problems
Here, we give a detailed summary of the failed problems. We use the form "index problem-name“ to refer to a problem in the following. 

29 FLETCBV2, 33 S308NE, 208 MOREBV satisfies the termination condition from the beginning.

Both fail case 1, TSD has Inf, NaN in the gradient, BBQ reaches the maxium iteration number 200,000: 9  SBRYBND, 10  ARGLINC, 26  CERI651ELS, 42  STRATEC, 53  VESUVIALS, 61  PALMER1C, 79  PALMER2C, 86  10FOLDTRLS, 95  SCURLY20, 102  SCURLY30, 112  PALMER1D, 118  CERI651DLS, 124  PENALTY2, 127  VESUVIOLS, 130  CERI651ALS, 131  VIBRBEAM, 138  CERI651CLS, 171  BROWNBS, 205  ARGLINB, 221  SCURLY10, 244  PALMER3C, 262  PALMER4C, 273  SCOSINE, 279  CERI651BLS.

Both fail case 2, BBQ’s line search fails, and TSD reaches the maxium iteration number 10,000: 7  FBRAIN3LS, 60  LSC2LS, 72  VESUVIOULS, 151  INDEF, 172  SSI.

Both fail case 3, BBQ’s line search fails, and TSD has Inf, NaN in the gradient: 63  POWELLBSLS, 190  DANWOODLS, 193  RAT42LS, 207  PARKCH, 225  DEVGLA2, 251  ROSZMAN1LS, 265  MISRA1CLS.

TSD has NaN/Inf, but BBQ works case: 34  PENALTY1, 45  DQRTIC, 59  MARATOSB, 116  VARDIM, 123  POWER, 125  BA-L1SPLS, 134  DJTL, 142  OSBORNEA, 158  MEXHAT, 170  DENSCHND, 175  KSSLS, 187  BENNETT5LS, 196  DEVGLA1, 198  STREG, 200  AKIVA, 215  HIELOW, 235  QUARTC, 245  SSBRYBND.

