# False alarm coincidence probability

A general formula for probability that is used in the estimation of significance of the coincidences is implemented. The code is available at [github](https://github.com/mbejger/polgraw-allsky/tree/master). Run `git clone https://github.com/mbejger/polgraw-allsky.git` to get the repository.

#### Prerequisites 

The code is written in standard `C`. The only dependency is [GNU Scientific Library (GSL)](http://www.gnu.org/software/gsl/), used to manipulate the Fisher matrix (calculate the eigenvectors and eigenvalues), the $\Gamma$ function and for the combinations. 
 
### Theoretical description 

This description is a short Appendix version of [Implementation of an F-statistic all-sky search for continuous gravitational waves in Virgo VSR1 data](http://iopscience.iop.org/0264-9381/31/16/165014/). 

For a given frequency band we analyze $L$ non-overlapping time segments: the search in the $l$th segment produces $N_l$ candidates. The size of the parameter space for each time segment is the same and it can be divided into the number $N_{\rm cell}$ of independent cells. The code tests the null hypothesis that the coincidences among candidates from $L$ segments are accidental.
The probability for a candidate event to fall into any given coincidence cell is equal to $1/N_{\rm cell}$. The probability $\epsilon_l$ that a given coincidence cell is populated with one or more candidate events is given by
$$
\epsilon_l = 1 - \Big(1 - \frac{1}{N_{\rm cell}}\Big)^{N_l},  
\quad\text{and for independent candidates (one candidate in one cell) it is}\quad  
\epsilon_l = \frac{N_l}{N_{\rm cell}}. 
$$
For two or more candidates within a given cell we choose the one with the highest signal-to-noise ratio. The probability $p_F(N_{\rm cell})$ that any given coincidence cell out of the total of $N_{\rm cell}$ cells contains candidate events from $C_{max}$ or more distinct data segments is given by a generalized binomial distribution: 
\begin{eqnarray}
p_F(N_{\rm cell}) &=& \sum_{n=C_{max}}^{L} \frac{1}{n!(L-n)!} \times \sum_{\sigma\in\Pi(L)} \epsilon_{\sigma(1)}\ldots\epsilon_{\sigma(n)}(1-\epsilon_{\sigma(n+1)})\ldots(1-\epsilon_{\sigma(L)}),
\end{eqnarray}
where $\sum_{\sigma \in \Pi(L)}$ is the sum over all permutations of $L$ data sequences.
Finally the probability $P_F$ that there is ${\mathcal C}_{max}$ or more coincidences
in one or more of the $N_{\rm cell}$ cells is
\begin{equation}
P_F = 1 - \left(1 - p_F(N_{\rm cell})\right)^{N_{\rm cell}}.
\end{equation} 

In order to find coincidences the entire cell coincidence grid is shifted by half a cell width in all possible $2^4 = 16$ combinations of the four parameter-space dimensions of $(f, \dot{f}, \delta, \alpha)$, and coincidences are searched in all the 16 coincidence grids. It does not account for cases when candidate events are located on opposite sides of cell borders, edges, and corners. This leads to a higher number of accidental coincidences, and consequently it underestimates the false alarm probability.

In the four dimension parameter space of $(f, \dot{f}, \delta, \alpha)$ the formula for the probability $P^{\rm shifts}_F$ that there are ${\mathcal C}_{\mathrm{max}}$ or more independent coincidences in one or more of the $N_{\rm cell}$ cells in all 16 grid shifts is 
\begin{eqnarray}
\label{eq:FAPs}
P^{\rm shifts}_F = 1 - \bigg[ 1 - \Big( 2^4 p_F(N_c) 
- \Big( {4 \choose 1} p_F(2 N_c) + {4 \choose 2} p_F(2^2 N_c) 
+ {4 \choose 3} p_F(2^3 N_c) + {4 \choose 4} p_F(2^4 N_c)   \Big) \\ \nonumber  
- \Big({4 \choose 2} p_F(2^2 N_c) + {4 \choose 3} p_F(2^3 N_c) 
+ {4 \choose 4} p_F(2^4 N_c) \Big) 
- \Big( {4 \choose 3} p_F(2^3 N_c) + {4 \choose 4} p_F(2^4 N_c) \Big) 
- {4 \choose 4} p_F(2^4 N_c)\Big)  \bigg]^{N_c}.
\end{eqnarray}

By choosing a certain false alarm probability $P_F$, we can calculate the threshold number ${\mathcal C}_{\mathrm{max}}$ of coincidences. If we obtain more than ${\mathcal C}_{\mathrm{max}}$ coincidences, the null hypothesis that coincidences are accidental is rejected at the significance level of $P_F$.


### Compilation

Run `make fap`; resulting binary is called `fap` (modify the `Makefile` to fit your system).

### How to run the program?

Sample call of `fap` is
```
> ./fap -nod 6 -band 0600 -data test-fap/summary_0600_1.txt -grid test-fap/ -vetofrac 0.222837 -cellsize 4 -threshold 0.2 
```
Number of days in the time segment `nod` equals 6, fraction of the band vetoed `vetofrac` is 0.222837 and the cell size `cellsize` is 4. Directory containing the grid matrix file `grid.bin` should be pointed out by the `grid` switch. The input data file `test-fap/summary_0600_1.txt` is a line from the `coincidences` output: 
```
0600_1 1110 155.312500     8     5  1.46199508e+00 1.30224869e-09 1.09355692e+00 2.08138626e+00 1.198182e+01 3 2020252 1969390 6 2137490 2038513 14 1752761 1666583 8 5526511 3334399 5 1827438 1790517
```

#### Full list of switches 

Type 
```
> ./fap --help 
```
to obtain the following description: 

| Switch          | Description       |
|-----------------|:------------------|
|-band            | Band number
|-cellsize        | Cell size (default value: 4)
|-data            | Coincidence summary file
|-grid            | Grid matrix directory (default value: .)
|-dt              | Data sampling time dt (default value: 2)
|-threshold       | FAP threshold (default value: 0.1)
|-nod             | Number of days
|-vetofrac        | Vetoed fraction of the band (default value: 0)

Also:

|                 |             | 
|-----------------|:------------|
| --help          |This help    |

### Output

Output to the screen (`stdout`) in the case of 
```
> ./fap -nod 6 -band 0600 -data test-fap/summary_0600_1.txt -grid test-fap/ -vetofrac 0.222837 -cellsize 4 -threshold 0.2 
```
is
```
Number of days in time segments: 6
Input data: test-fap/summary_0600_1.txt
Grid matrix data directory: test-fap/
Band number: 0600 (veto fraction: 0.222837)
The reference frequency fpo: 155.312500
The data sampling time dt: 2.000000
FAP threshold: 0.200000
Cell size: 4
```
whereas the result (if `threshold` is reached) is printed to `stderr`: 
```
#band f_min        f_max        PFshifts      noc Nkall    f s d a hemisphere 
0600  1.553125e+02 1.555625e+02 1.810747e-01  5   10799402 1.461995e+00 1.302249e-09 1.093557e+00 2.081386e+00 1.198182e+01 1
```
