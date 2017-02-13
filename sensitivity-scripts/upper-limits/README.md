The summary of simulations for a given band (`0666`, say) processed by [mcsummary.sh](https://github.com/mbejger/polgraw-allsky/blob/master/sensitivity-scripts/mcsummary.sh) result in the following list of `h0` amplitudes followed by the corresponding fractions of significant coincidences (`N_coin/N`):
```bash
0666 8e-02 9e-02 1e-01 1.25e-01 1.5e-01 2e-01 3e-01 0.59 0.74 0.79 0.94 1.0 1.0 1.0
```
or, equivalently, in a column format, processed by [rewrite_and_sort.py](https://github.com/mbejger/polgraw-allsky/blob/master/sensitivity-scripts/upper-limits/rewrite_and_sort.py):
```bash
$ python rewrite_and_sort.py <(echo "0666 8e-02 9e-02 1e-01 1.25e-01 1.5e-01 2e-01 3e-01 0.59 0.74 0.79 0.94 1.0 1.0 1.0") > 0666_results
```
```bash
$ cat 0666_results
band h   ul
0666 0.080 0.59
0666 0.090 0.74
0666 0.100 0.79
0666 0.125 0.94
0666 0.150 1.00
0666 0.200 1.00
0666 0.300 1.00
```
We are interested in a 95% upper limit i.e. the `h0` corresponding to the fraction 0.95 of significant coincidences in the simulation (`N_coin/N=0.95`). This is obtained by fitting a sigmoid function

```python
def sigmoid(x, x0, k):
     y = 1.0 / (1.0 + np.exp(k*(x0-x)))
     return y
```
to the above data. 'Good' initial guesses of `x0` (`x` value where the slope is steep) and `k` are in this case `x0=0.01` and `k=50`. For bands lower than `0100`, `x0` is around 1. Fitting is done by [ul.py](https://github.com/mbejger/polgraw-allsky/blob/master/sensitivity-scripts/upper-limits/ul.py):

```bash
$ python ul.py 0666_results 0666 0.01 test.pdf
0666 1.2816e-01
```
The output is the band number and `h0` corresponding to the 95% upper limit. Last command-line option [test.pdf](https://github.com/mbejger/polgraw-allsky/blob/master/sensitivity-scripts/upper-limits/test.pdf) is optional. It produces the auxiliary plot, with the 95% upper limit is denoted by red circle.

