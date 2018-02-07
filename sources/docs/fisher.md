# Fisher matrix calculation

The Fisher matrix associated with the signal model and its inversion is calculated using [this code](https://github.com/mbejger/polgraw-allsky/tree/master/search/network/src-cpu/fisher.c). 

### Prerequisites

The code is written in standard `C` and it's mostly based on functions used in [search](../polgraw-allsky/search_for_candidates/). Arbitrary-precision interval arithmetic [Arb](http://arblib.org) library is used to invert the (usually) not-very-well posed Fisher matrix, so it has to be installed beforehand. [Arb](http://arblib.org) requires [FLINT](http://www.flintlib.org), [MPFR](http://www.mpfr.org), and either [MPIR](http://www.mpir.org) or [GMP](http://www.gmplib.org). 

### Compilation 

Run  `make fisher` in `search/network/src-cpu` - resulting  binary is called `fisher`. Modify the `Makefile` (especially the variable `ARB_DIR`) to fit your system. 

####
#### Full list of switches

For the full list of options, type

```bash 
% ./fisher --help 
```

| Switch      | Description       |
|-------------|:------------------| 
|-data        | Data directory (default is `.`)
|-ident       | Frame number
|-band        | Band number
|-fpo         | Reference band frequency `fpo` value
|-dt          | Data sampling time dt (default value: `0.5`)
|-usedet      | Use only detectors from string (default is `use all available`)
|-addsig      | Add signal with parameters from `<file>`

Also: 

|                 |             | 
|-----------------|:------------|
| --help          | This help | 


### Example

Minimal call to `fisher` is as follows: 

```
% ./fisher -data 2d_0.25 -ident 001 -band 1234 -usedet H1 -dt 2 -nod 2 -addsig sigfile
```

where
 
* `data` is the base directory of input data files (e.g., [this Gaussian data](https://polgraw.camk.edu.pl/H1L1_2d_0.25.tar.gz)),
*  Sampling time `dt` is $2 s$, 
* `ident` is the number of time frame to be analyzed ($001$),
* `nod` number of days is $2$, 
* `band` is the number of the frequency band (see the [input data structure](../polgraw-allsky/input_data) for details). 
* `usedet` switch to chose a detector (here $H1$) 
* `addsig` switch to chose a file with signal data

The `sigfile` file consists of 8 numbers: 

- frequency [radians, between 0 and $\pi$] above `fpo`  
- spindown (frequency time derivative) [$\mathrm{Hz/s}$]  
- declination [radians, between $\pi$ and $-\pi$]
- right ascension [radians, between 0 and $2\pi$]
- 4 amplitudes $a_1, a_2, a_3, a_4$ 

e.g., 
```bash 
1.431318175386891
-7.9539e-9
0.6363615896875658
4.396884357060633
7.764354801848407e-3
-1.422468474545797e-2
-1.559826840666228e-2
-8.623005535014139e-3
```

The amplitudes $a_1, a_2, a_3, a_4$ correspond to the signal 
amplitude model 

$$ 
h = a_1 h_1 + a_2 h_2 + a_3 h_3 + a_4 h_4, 
$$ 


where 

$$
h(t) = \left(a_1 a(t) + a_2 b(t)\right)\cos(\psi) 
+ \left(a_3 a(t) + a_4 b(t)\right)sin(\psi), 
$$

with $\psi(f, \dot{f}, \delta, \alpha, t)$ being the phase of the signal, and $a$ and $b$ the amplitude modulation functions (calculated in the [modvir](https://github.com/mbejger/polgraw-allsky/tree/master/search/network/src-cpu/settings.c) function). 

### Example output 

```bash 
Number of days is 2
Input data directory is 2d_0.25
Frame and band numbers are 1 and 1234
The reference frequency fpo is 308.859375
The data sampling time dt is 2.000000
Adding signal from 'sigfile'
Settings - number of detectors: 1
Using H1 IFO as detector #0... 2d_0.25/001/H1/xdatc_001_1234.bin as input time series data
Using 2d_0.25/001/H1/DetSSB.bin as detector H1 ephemerids...
The Fisher matrix:
1.4602194451385117e+10 9.6224528459395950e+14 1.0472943290223141e+10 1.9109038902196317e+11 -5.6465717962684985e+06 -4.1172082782989331e+06 -2.9517981884375727e+06 7.0470722915177587e+06 
9.6224528459395950e+14 6.7134859540063060e+19 6.5156368686556762e+14 1.1289208033400466e+16 -3.3340781322676825e+11 -2.4401458951787103e+11 -1.7381344007331134e+11 4.1673668536337537e+11 
1.0472943290223141e+10 6.5156368686556762e+14 8.0545511145922174e+09 1.5540065972734366e+11 -4.6112260504154256e+06 -3.4059196396034071e+06 -2.3414307185722976e+06 5.7018637665277841e+06 
1.9109038902196317e+11 1.1289208033400466e+16 1.5540065972734366e+11 3.1204237328935298e+12 -9.2852217833914995e+07 -6.9173973293523684e+07 -4.6213130625602841e+07 1.1410003283718885e+08 
-5.6465717962684985e+06 -3.3340781322676825e+11 -4.6112260504154256e+06 -9.2852217833914995e+07 7.6814199126393523e+03 -1.6833193661163169e-03 0.0000000000000000e+00 0.0000000000000000e+00 
-4.1172082782989331e+06 -2.4401458951787103e+11 -3.4059196396034071e+06 -6.9173973293523684e+07 -1.6833193661163169e-03 1.0351065428550139e+04 0.0000000000000000e+00 0.0000000000000000e+00 
-2.9517981884375727e+06 -1.7381344007331134e+11 -2.3414307185722976e+06 -4.6213130625602841e+07 0.0000000000000000e+00 0.0000000000000000e+00 7.6814199126393523e+03 -1.6833193661163169e-03 
7.0470722915177587e+06 4.1673668536337537e+11 5.7018637665277841e+06 1.1410003283718885e+08 0.0000000000000000e+00 0.0000000000000000e+00 -1.6833193661163169e-03 1.0351065428550139e+04 
Inverting the Fisher matrix...
Diagonal elements of the covariance matrix:
2.561275e-04 3.867944e-18 2.363103e-03 9.137957e-04 1.343874e+05 4.107046e+04 3.329715e+04 1.117611e+05 
```


