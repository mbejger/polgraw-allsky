polgraw-allsky
==============

All-sky almost-monochromatic gravitational-wave pipeline (Polgraw group)

See the [documentation page](http://mbejger.github.io/polgraw-allsky/) for more details (documentation in progress).

Contributors
------------ 
In alphabetic order: 
* Michał Bejger
* Jan Bolek
* Paweł Ciecieląg
* Orest Dorosh
* Aleksander Garus
* Andrzej Królak
* Maciej Piętka
* Andrzej Pisarski 
* Gevorg Poghosyan
* Rafał Skrzypiec


This repository contains several elements of the pipeline: 

* Generation of the initial data (gen2day)
* Parameter space grid generation (gridgen)
* Search for candidate signals (search) 

This code searches for continuous gravitational wave candidate 
signals in time-domain data using the F-statistic on data 
from a network of detectors (CPU implementation). Currently the 
available detectors are H1 (Hanford), L1 (Livingston) of LIGO 
and V1 (Cascina) of Virgo.  

Other versions are the network-of-detectors [spotlight](https://github.com/mbejger/polgraw-allsky/tree/master/search/spotlight) search, 
[one detector](https://github.com/mbejger/polgraw-allsky/tree/master/search/one-detector) version 
and the [GPU](https://github.com/mbejger/polgraw-allsky/tree/master/search/network/src-gpu) version. 

* Coincidences between candidates in different time frames (coincidences)
* Followup of promissing signals (followup) 
* test data 


