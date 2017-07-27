# False alarm probability of coincidences

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

### Full list of switches 

To obtain the full list of options, type 
```
% ./fap --help 
```

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

### Example 

Using the software injection added to 2-day Gaussian noise data segments (see [minimal example of the pipeline](../polgraw-allsky/pipeline_script)): 
```bash 
% ./fap -nod 2 -band 1234 -data <(sort -gk5 -gk10 summary | tail -1) -grid ../../testdata/2d_0.25/004 -vetofrac 0.0 -cellsize 4 -threshold 1.0 
```
or, with the auxilary `fap.sh` script, 
```bash 
% band=1234; bash fap.sh <(sort -gk5 -gk10 summary | tail -1) <(echo $band 0.0) ../../testdata/2d_0.25/004 
```
Number of days in the time segment `nod` equals 2, fraction of the band vetoed `vetofrac` is 0 (no lines, Gaussian data) and the cell size scalling factor `cellsize` is 4. Directory containing the grid matrix file `grid.bin` of the reference frame (in this case frame `004`) should be given by the `grid` switch. The input data is the last line of a sorted `summary` file to select the shift giving the best coincidence with the highest signal-to-noise ratio: 
```
1234_2 1111 308.859375     8     5  9.95663703e-01 -1.10830358e-09 -1.12585347e-01 1.97463002e+00 1.246469e+01 5 2040 1987 1 2483 2419 4 2384 2193 3 2247 2137 8 2408 2363 2 2249 2172 6 2305 2220 7 2226 2191 6 2 8 3 5
```
(see the [coincidences](../polgraw-allsky/coincidences) section for details). 

### Output

```bash 
% ./fap -nod 2 -band 1234 -data <(sort -gk5 -gk10 summary | tail -1) -grid ../../testdata/2d_0.25/004 -vetofrac 0.0 -cellsize 4 -threshold 1.0 
```
is
```bash 
Number of days in time segments: 2
Input data: /dev/fd/63
Grid matrix data directory: ../../testdata/2d_0.25/004
Band number: 1234 (veto fraction: 0.000000)
The reference frequency fpo: 308.859375
The data sampling time dt: 2.000000
FAP threshold: 1.000000
Cell size: 4
1234 3.088594e+02 3.091094e+02 7.665713e-08 5 17682 9.956637e-01 -1.108304e-09 -1.125853e-01 1.974630e+00 1.246469e+01 2
```
The last line (in case the probability `PFshifts` is lower than the `threshold`) is printed to `stderr`. The meaning of this output is the following:  
```
#band f_min       f_max        PFshifts     noc Nkall  f s d a hemisphere 
1234 3.088594e+02 3.091094e+02 7.665713e-08 5   17682  9.956637e-01 -1.108304e-09 -1.125853e-01 1.974630e+00 1.246469e+01 2
```
