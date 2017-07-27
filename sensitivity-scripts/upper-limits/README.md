The summary of simulations for a given band (`0165`, say) processed by the [summary.sh](https://github.com/mbejger/polgraw-allsky/blob/master/sensitivity-scripts/summary.sh) result in the following list of `h0` amplitudes followed by the corresponding fractions of significant coincidences (`N_coin/N`):

```bash
band h   ul 
0165 0.150 0.61
0165 0.200 0.78
0165 0.250 0.95
0165 0.300 0.99
```

We are interested in a 95% upper limit i.e. the `h0` corresponding to the fraction 0.95 of significant coincidences in the simulation (`N_coin/N=0.95`). This is obtained by fitting a sigmoid function

```python
def sigmoid(x, x0, k):
     y = 1.0 / (1.0 + np.exp(k*(x0-x)))
     return y
```
to the above data. Fitting is done by [ul.py](https://github.com/mbejger/polgraw-allsky/blob/master/sensitivity-scripts/upper-limits/ul.py):

```bash
% python ul.py 0165_results 0165 0.01 test.pdf
0165 2.7026e-01
```
The output is the band number and `h0` corresponding to the 95% upper limit. Last command-line option [test.pdf](https://github.com/mbejger/polgraw-allsky/blob/master/sensitivity-scripts/upper-limits/test.pdf) is optional. It produces the auxiliary plot, with the 95% upper limit is denoted by red circle.  
