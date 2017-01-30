# Coincidences between candidates 

In order to establish the probability of detection of a real gravitational wave, after finding the candidate signals with the [F-statistic candidate signal search](../polgraw-allsky/candidate_search/), the pipeline searches for coincidences between candidates found in different time frames. 

The coincidences code is available at [github](https://github.com/mbejger/polgraw-allsky/tree/master). To get it, run `git clone https://github.com/mbejger/polgraw-allsky.git`.

#### Prerequisites 

The code is written in standard `C`. The only dependency is [GNU Scientific Library (GSL)](http://www.gnu.org/software/gsl/), used to manipulate the Fisher matrix (calculate the eigenvectors and eigenvalues). [GNU struct dirent](http://www.gnu.org/software/libc/manual/html_node/Accessing-Directories.html#Accessing-Directories) objects are used to read the directories. 
 
#### The idea behind coincidences 

After finding the candidate signals in different time frames (`search`), we want to confirm the existence of signals with the same parameters along the span of time segments. to further perform a validation search for high-coincidence, or otherwise interesting candidates (3rd stage, to be done). To do this, the candidates from a list of trigger files (time frames) are read, and for each trigger file

* a candidate is transformed to a defined time moment (frequency shifted to a reference frame), 
* translated into linear (integer) coordinates, 
* duplicates within each frame are removed, 
* list of unique candidates from all the frames is created and sorted, 
* duplicates are counted - these are the coincidences. 
```
[-x----][-x----][-w----][-x----][-----t][-x----][-x----][---f--]
5 coincidences of x in 8 data segments, not bad...
```
*TODO: describe cell shifts (16 different shifts in total: 0101, 0110, 0010 etc. in f, s, d, a directions) and scaling of cells (used to define the linear parameters for a given cell to subsequently compare the candidate values)* 

### 7. Compilation

Run  `make coincidences`. Resulting binary is called `coincidences`. Modify the `Makefile` to fit your system. By default the `YEPPP!` library is selected, but `coincidences` does not need it (select `GNUSINCOS`, for example). 

### 8. How to run the program?

Sample call of `coincidences` is as follows (code compiled with the `GNUSINCOS` option):
```
> ./coincidences -data ./triggers_J0154+2706 -trigname J0154+2706 -refloc ./coinc-testdata -fpo 687.656250 -shift 0010 -scale 4844 -refr 100 -dt 2 -mincoin 15 -output results -narrowdown 0.2
```
#### 8.1. Full list of switches 
Type 
```
> ./coincidences --help 
```
to obtain the following description: 

| Switch          | Description       |
|-----------------|:------------------|
|-data            | Data directory (default is `./candidates`)
|-output          | Output directory (default is `./coinc-results`)
|-shift           | Cell shifts in `fsda` directions (4 digit number, e.g. `0101`, default `0000`)
|-scale           | Cell scaling in `fsda` directions (4 digit number, e.g. `4824`, default `1111`)
|-refr            | Reference frame number
|-fpo             | Reference band frequency `fpo` value
|-dt              | Data sampling time dt (default value: `0.5`)
|-trigname        | Part of triggers' name (for identifying files)
|-refloc          | Location of the reference grid.bin and starting_date files
|-mincoin         | Minimal number of coincidences recorded
|-narrowdown      | Narrow-down the frequency band (range [0, 0.5] +- around center)
|-snrcutoff       | Signal-to-noise threshold cutoff (default value: 6)

Also:

|                 |             |
|-----------------|:------------|
| --help          |This help    |

### 9. Output
Output to the screen (`stdout`) in the case of 
```
> ./coincidences -data ./triggers_J0154+2706 -trigname J0154+2706 -refloc ./coinc-testdata -fpo 687.656250 -shift 0010 -scale 4844 -refr 100 -dt 2 -mincoin 15 -output results -narrowdown 0.2
```
is
```
Reading the reference grid.bin at ./coinc-testdata
fftpad from the grid file: 1
Settings dt: 2.000000, oms: 8641.343293
Reference frame number: 100
Cell shifts  in f, s, d, a directions: 0 0 1 0 
Cell scaling in f, s, d, a directions: 4 8 4 4
Reading triggers_090_000_J0154+2706_1.bin... Frame 90: 3765/7533
Reading triggers_178_000_J0154+2706_1.bin... Frame 178: 10934/15676
Reading triggers_055_000_J0154+2706_1.bin... Frame 55: 4045/4568
Reading triggers_186_000_J0154+2706_1.bin... Frame 186: 27045/39588
Reading triggers_066_000_J0154+2706_1.bin... Frame 66: 15250/40642
Reading triggers_091_000_J0154+2706_1.bin... Frame 91: 10655/20453
Reading triggers_074_000_J0154+2706_1.bin... Frame 74: 13139/22727
Reading triggers_129_000_J0154+2706_1.bin... Frame 129: 11026/27644
Reading triggers_164_000_J0154+2706_1.bin... Frame 164: 13731/18744
Reading triggers_058_000_J0154+2706_1.bin... Frame 58: 28268/113006
Reading triggers_147_000_J0154+2706_1.bin... Frame 147: 8715/10321
Reading triggers_163_000_J0154+2706_1.bin... Frame 163: 16692/27502
Reading triggers_063_000_J0154+2706_1.bin... Frame 63: 19239/65766
Reading triggers_131_000_J0154+2706_1.bin... Frame 131: 14185/27568
Reading triggers_187_000_J0154+2706_1.bin... Frame 187: 22569/38268
Reading triggers_001_000_J0154+2706_1.bin... Frame 1: 16258/23323
Reading triggers_075_000_J0154+2706_1.bin... Frame 75: 13374/22863
Reading triggers_002_000_J0154+2706_1.bin... Frame 2: 19019/25287
Reading triggers_204_000_J0154+2706_1.bin... Frame 204: 16244/20967
Reading triggers_005_000_J0154+2706_1.bin... Frame 5: 20565/28894
Reading triggers_198_000_J0154+2706_1.bin... Frame 198: 12896/23195
Reading triggers_140_000_J0154+2706_1.bin... Frame 140: 12201/28702
Reading triggers_128_000_J0154+2706_1.bin... Frame 128: 28354/43979
Reading triggers_188_000_J0154+2706_1.bin... Frame 188: 14069/19071
Reading triggers_079_000_J0154+2706_1.bin... Frame 79: 9647/21067
Reading triggers_044_000_J0154+2706_1.bin... Frame 44: 19454/25459
Reading triggers_118_000_J0154+2706_1.bin... Frame 118: 7177/11786
Reading triggers_195_000_J0154+2706_1.bin... Frame 195: 11244/16600
Reading triggers_041_000_J0154+2706_1.bin... Frame 41: 13505/16568
Reading triggers_123_000_J0154+2706_1.bin... Frame 123: 19189/35112
Reading triggers_057_000_J0154+2706_1.bin... Frame 57: 26405/136897
Reading triggers_106_000_J0154+2706_1.bin... Frame 106: 28089/90068
Reading triggers_207_000_J0154+2706_1.bin... Frame 207: 14298/17108
Reading triggers_190_000_J0154+2706_1.bin... Frame 190: 24976/39798
Reading triggers_199_000_J0154+2706_1.bin... Frame 199: 10825/14307
Reading triggers_021_000_J0154+2706_1.bin... Frame 21: 21062/31870
Reading triggers_051_000_J0154+2706_1.bin... Frame 51: 32053/148466
Reading triggers_153_000_J0154+2706_1.bin... Frame 153: 17958/26136
Reading triggers_004_000_J0154+2706_1.bin... Frame 4: 32965/58003
Reading triggers_159_000_J0154+2706_1.bin... Frame 159: 18799/30658
Reading triggers_125_000_J0154+2706_1.bin... Frame 125: 13853/22953
Reading triggers_132_000_J0154+2706_1.bin... Frame 132: 9116/16824
Reading triggers_167_000_J0154+2706_1.bin... Frame 167: 10366/14077
Reading triggers_023_000_J0154+2706_1.bin... Frame 23: 30762/48433
Reading triggers_183_000_J0154+2706_1.bin... Frame 183: 22550/43246
Reading triggers_120_000_J0154+2706_1.bin... Frame 120: 7671/12633
Reading triggers_197_000_J0154+2706_1.bin... Frame 197: 15657/23986
Reading triggers_087_000_J0154+2706_1.bin... Frame 87: 4580/7618
Reading triggers_148_000_J0154+2706_1.bin... Frame 148: 12294/19698
Reading triggers_116_000_J0154+2706_1.bin... Frame 116: 5102/7457
Reading triggers_061_000_J0154+2706_1.bin... Frame 61: 33001/116930
Reading triggers_089_000_J0154+2706_1.bin... Frame 89: 9472/29583
Reading triggers_160_000_J0154+2706_1.bin... Frame 160: 19104/34216
Reading triggers_012_000_J0154+2706_1.bin... Frame 12: 12069/20261
Reading triggers_157_000_J0154+2706_1.bin... Frame 157: 15953/28418
Reading triggers_136_000_J0154+2706_1.bin... Frame 136: 10848/17265
Reading triggers_068_000_J0154+2706_1.bin... Frame 68: 21157/43332
Reading triggers_115_000_J0154+2706_1.bin... Frame 115: 8311/16111
Reading triggers_210_000_J0154+2706_1.bin... Frame 210: 26596/37238
Reading triggers_014_000_J0154+2706_1.bin... Frame 14: 11475/17309
Reading triggers_003_000_J0154+2706_1.bin... Frame 3: 20306/30310
Reading triggers_077_000_J0154+2706_1.bin... Frame 77: 12905/59221
Reading triggers_176_000_J0154+2706_1.bin... Frame 176: 26268/65780
Reading triggers_110_000_J0154+2706_1.bin... Frame 110: 5381/13017
Reading triggers_015_000_J0154+2706_1.bin... Frame 15: 9253/13996
Reading triggers_169_000_J0154+2706_1.bin... Frame 169: 44430/83664
Reading triggers_047_000_J0154+2706_1.bin... Frame 47: 9216/21530
Reading triggers_149_000_J0154+2706_1.bin... Frame 149: 14636/32072
Reading triggers_145_000_J0154+2706_1.bin... Frame 145: 7953/10547
Reading triggers_082_000_J0154+2706_1.bin... Frame 82: 9759/19030
Reading triggers_024_000_J0154+2706_1.bin... Frame 24: 26444/45609
Reading triggers_093_000_J0154+2706_1.bin... Frame 93: 4030/9100
Reading triggers_150_000_J0154+2706_1.bin... Frame 150: 12250/16337
Reading triggers_114_000_J0154+2706_1.bin... Frame 114: 7470/13938
Reading triggers_080_000_J0154+2706_1.bin... Frame 80: 12001/23822
Reading triggers_177_000_J0154+2706_1.bin... Frame 177: 15090/26795
Reading triggers_117_000_J0154+2706_1.bin... Frame 117: 12307/20300
Reading triggers_175_000_J0154+2706_1.bin... Frame 175: 4524/5422
Reading triggers_062_000_J0154+2706_1.bin... Frame 62: 19430/43939
Reading triggers_205_000_J0154+2706_1.bin... Frame 205: 14109/18186
Reading triggers_138_000_J0154+2706_1.bin... Frame 138: 9164/15275
Reading triggers_121_000_J0154+2706_1.bin... Frame 121: 7055/16649
Reading triggers_013_000_J0154+2706_1.bin... Frame 13: 13437/18421
Reading triggers_053_000_J0154+2706_1.bin... Frame 53: 14543/20233
Reading triggers_154_000_J0154+2706_1.bin... Frame 154: 12496/16592
Reading triggers_018_000_J0154+2706_1.bin... Frame 18: 12559/19815
Reading triggers_042_000_J0154+2706_1.bin... Frame 42: 14353/18348
Reading triggers_011_000_J0154+2706_1.bin... Frame 11: 13332/19033
Reading triggers_185_000_J0154+2706_1.bin... Frame 185: 15324/22440
Reading triggers_105_000_J0154+2706_1.bin... Frame 105: 3828/7795
Reading triggers_144_000_J0154+2706_1.bin... Frame 144: 13529/26921
Reading triggers_191_000_J0154+2706_1.bin... Frame 191: 23632/40532
Reading triggers_065_000_J0154+2706_1.bin... Frame 65: 34496/98915
Reading triggers_209_000_J0154+2706_1.bin... Frame 209: 16456/20729
Reading triggers_045_000_J0154+2706_1.bin... Frame 45: 12352/25900
Reading triggers_130_000_J0154+2706_1.bin... Frame 130: 7874/15641
Reading triggers_161_000_J0154+2706_1.bin... Frame 161: 13463/19649
Reading triggers_181_000_J0154+2706_1.bin... Frame 181: 13530/19242
Reading triggers_022_000_J0154+2706_1.bin... Frame 22: 12618/18092
Reading triggers_141_000_J0154+2706_1.bin... Frame 141: 26407/56758
Reading triggers_173_000_J0154+2706_1.bin... Frame 173: 13824/20601
Reading triggers_137_000_J0154+2706_1.bin... Frame 137: 10673/20096
Reading triggers_206_000_J0154+2706_1.bin... Frame 206: 6985/8183
Reading triggers_056_000_J0154+2706_1.bin... Frame 56: 9745/12755
Reading triggers_184_000_J0154+2706_1.bin... Frame 184: 14660/26148
Reading triggers_072_000_J0154+2706_1.bin... Frame 72: 18750/52603
Reading triggers_078_000_J0154+2706_1.bin... Frame 78: 12204/19397
Reading triggers_007_000_J0154+2706_1.bin... Frame 7: 19871/27320
Reading triggers_059_000_J0154+2706_1.bin... Frame 59: 23289/107891
Reading triggers_162_000_J0154+2706_1.bin... Frame 162: 27353/65940
Reading triggers_170_000_J0154+2706_1.bin... Frame 170: 11087/14251
Reading triggers_180_000_J0154+2706_1.bin... Frame 180: 9352/11327
Reading triggers_124_000_J0154+2706_1.bin... Frame 124: 7468/19457
Reading triggers_189_000_J0154+2706_1.bin... Frame 189: 19688/36281
Reading triggers_168_000_J0154+2706_1.bin... Frame 168: 20320/27800
Reading triggers_119_000_J0154+2706_1.bin... Frame 119: 11197/21976
Reading triggers_203_000_J0154+2706_1.bin... Frame 203: 13601/17558
Reading triggers_108_000_J0154+2706_1.bin... Frame 108: 4874/12027
Reading triggers_142_000_J0154+2706_1.bin... Frame 142: 9806/21425
Reading triggers_171_000_J0154+2706_1.bin... Frame 171: 21826/39658
Reading triggers_165_000_J0154+2706_1.bin... Frame 165: 14057/20140
Reading triggers_094_000_J0154+2706_1.bin... Frame 94: 4326/7880
Reading triggers_158_000_J0154+2706_1.bin... Frame 158: 19553/28134
Reading triggers_174_000_J0154+2706_1.bin... Frame 174: 20472/28609
Reading triggers_016_000_J0154+2706_1.bin... Frame 16: 12721/17353
Reading triggers_071_000_J0154+2706_1.bin... Frame 71: 15500/26422
Reading triggers_133_000_J0154+2706_1.bin... Frame 133: 12714/23266
Reading triggers_088_000_J0154+2706_1.bin... Frame 88: 4691/8545
Reading triggers_196_000_J0154+2706_1.bin... Frame 196: 19311/26083
Reading triggers_006_000_J0154+2706_1.bin... Frame 6: 15138/22097
Reading triggers_208_000_J0154+2706_1.bin... Frame 208: 7803/9154
Reading triggers_084_000_J0154+2706_1.bin... Frame 84: 6187/10823
Reading triggers_083_000_J0154+2706_1.bin... Frame 83: 8113/17695
Reading triggers_155_000_J0154+2706_1.bin... Frame 155: 15076/21385
Reading triggers_143_000_J0154+2706_1.bin... Frame 143: 11628/19918
Reading triggers_135_000_J0154+2706_1.bin... Frame 135: 9137/15006
Reading triggers_107_000_J0154+2706_1.bin... Frame 107: 26982/77757
Reading triggers_172_000_J0154+2706_1.bin... Frame 172: 23742/44957
Reading triggers_092_000_J0154+2706_1.bin... Frame 92: 3740/11404
Reading triggers_054_000_J0154+2706_1.bin... Frame 54: 28481/125688
Reading triggers_017_000_J0154+2706_1.bin... Frame 17: 22677/38941
Reading triggers_073_000_J0154+2706_1.bin... Frame 73: 15433/35465
Reading triggers_067_000_J0154+2706_1.bin... Frame 67: 15925/32033
Reading triggers_010_000_J0154+2706_1.bin... Frame 10: 18828/29511
Reading triggers_113_000_J0154+2706_1.bin... Frame 113: 9137/21032
Reading triggers_008_000_J0154+2706_1.bin... Frame 8: 24052/34056
Reading triggers_152_000_J0154+2706_1.bin... Frame 152: 13927/25378
Reading triggers_146_000_J0154+2706_1.bin... Frame 146: 14085/20278
Reading triggers_086_000_J0154+2706_1.bin... Frame 86: 9693/23688
Reading triggers_182_000_J0154+2706_1.bin... Frame 182: 11481/16302
Reading triggers_139_000_J0154+2706_1.bin... Frame 139: 10928/21291
Reading triggers_134_000_J0154+2706_1.bin... Frame 134: 14389/23714
Reading triggers_050_000_J0154+2706_1.bin... Frame 50: 29735/127712
Reading triggers_111_000_J0154+2706_1.bin... Frame 111: 9733/19731
Reading triggers_201_000_J0154+2706_1.bin... Frame 201: 18172/26948
Reading triggers_112_000_J0154+2706_1.bin... Frame 112: 9651/21172
Reading triggers_095_000_J0154+2706_1.bin... Frame 95: 7536/16646
Reading triggers_156_000_J0154+2706_1.bin... Frame 156: 18557/29548
Reading triggers_109_000_J0154+2706_1.bin... Frame 109: 7876/17800
Reading triggers_019_000_J0154+2706_1.bin... Frame 19: 25890/51813
Reading triggers_060_000_J0154+2706_1.bin... Frame 60: 10356/17781
Reading triggers_046_000_J0154+2706_1.bin... Frame 46: 19287/25444
Reading triggers_076_000_J0154+2706_1.bin... Frame 76: 11102/18072
Reading triggers_166_000_J0154+2706_1.bin... Frame 166: 14167/29895
Reading triggers_202_000_J0154+2706_1.bin... Frame 202: 17350/28904
Reading triggers_052_000_J0154+2706_1.bin... Frame 52: 9601/14333
Reading triggers_127_000_J0154+2706_1.bin... Frame 127: 9316/14658
Reading triggers_151_000_J0154+2706_1.bin... Frame 151: 48093/96198
Reading triggers_193_000_J0154+2706_1.bin... Frame 193: 22411/31960
Reading triggers_070_000_J0154+2706_1.bin... Frame 70: 12978/20700
Reading triggers_179_000_J0154+2706_1.bin... Frame 179: 23480/55171
Reading triggers_081_000_J0154+2706_1.bin... Frame 81: 9145/23443
Reading triggers_126_000_J0154+2706_1.bin... Frame 126: 9155/15582
Reading triggers_122_000_J0154+2706_1.bin... Frame 122: 25245/43311
Reading triggers_192_000_J0154+2706_1.bin... Frame 192: 16635/31789
Reading triggers_069_000_J0154+2706_1.bin... Frame 69: 20628/35691
Reading triggers_020_000_J0154+2706_1.bin... Frame 20: 14195/22248
Reading triggers_049_000_J0154+2706_1.bin... Frame 49: 17332/35069
Reading triggers_200_000_J0154+2706_1.bin... Frame 200: 13238/23902
Reading triggers_048_000_J0154+2706_1.bin... Frame 48: 23554/94544
Reading triggers_009_000_J0154+2706_1.bin... Frame 9: 14132/21967
Reading triggers_085_000_J0154+2706_1.bin... Frame 85: 7641/13944
Reading triggers_194_000_J0154+2706_1.bin... Frame 194: 16616/22392
Reading triggers_043_000_J0154+2706_1.bin... Frame 43: 21511/74531
Total number of candidates from all frames: 2838849
J0154+2706 0010 687.656250   184   174    1.440594e+00 -1.220905e-09 4.334456e-01 5.040727e-01 2.190372e+02
```
Last line is printed to `stderr` (can be streamed to a `summary` file). The line contains the `trigname` identifier, the `shift` value, the `fpo` band frequency, the number of files read (184), and the maximum coincidence found (174). Next 5 numbers are arithmetic mean values of the frequency $\,f$ (in radians), frequency derivative $\,\dot{f}$ (spindown), sky positions $\delta$ and $\alpha$, and the mean signal-to-noise ratio, $\widetilde{\mathrm{snr}}=\sqrt{\sum_i \mathrm{snr}_i^2}$. 

Coincidences above `mincoin` are recorded in a binary file `.coi`. Each coincidence is a set of following numbers: 
$$
N,\quad\bar{f},\quad\bar{s},\quad\bar{d},\quad\bar{a},\quad\widetilde{\mathrm{snr}},\quad\mathrm{fr}_{1},\,\dots\,\mathrm{fr}_{N},\quad\mathrm{p}_{1},\,\dots\,\mathrm{p}_{N}
$$
where 

* $N$ is the size of coincidence (written as one `unsigned short int`), 
* $\bar{f}$, $\bar{a}$, $\bar{d}$, $\bar{a}$ and $\widetilde{\mathrm{snr}}$ are the mean parameters of the signal ($5\times$`float`),
* $\mathrm{fr}_{1},\,\dots\,\mathrm{fr}_{N}$ are frame numbers ($N\times$`unsigned short int`), 
* $\mathrm{p}_{1},\,\dots\,\mathrm{p}_{N}$ are positions of signals that take part in the coincidences, in their corresponding trigger files ($N\times$`int`).

The `master` branch contains an auxiliary `read_coi.c` code that reads the binary `.coi` files and prints the output to `stdout`.  
