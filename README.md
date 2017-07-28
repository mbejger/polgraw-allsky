polgraw-allsky
==============

All-sky almost-monochromatic gravitational-wave pipeline (Polgraw group)

See the [documentation page](http://mbejger.github.io/polgraw-allsky/) for more details.

Contributors
------------ 
In alphabetic order: 

* Pia Astone 
* Michał Bejger
* Jan Bolek
* Paweł Ciecieląg
* Orest Dorosh
* Aleksander Garus
* Andrzej Królak
* Máté Ferenc Nagy-Egri
* Maciej Piętka
* Andrzej Pisarski 
* Gevorg Poghosyan
* Magdalena Sieniawska 
* Rafał Skrzypiec


This repository gathers several elements of the pipeline: 

* Generation of the initial data: `gen2day`
* Parameter space grid generation:`gridgen`
* Coherent search for candidate signals in narrowband time segments: serial CPU version (`gwsearch-cpu`) and openMP version 
* Search for coincidences among different time segments: `coincidences` 
* Estimation of false alarm probability of coincidences: `fap`
* Followup of interesting outliers: `followup`
* Sensitivity upper limits: in `sensitivity-scripts` 
* [Test data](https://polgraw.camk.edu.pl/H1L1_2d_0.25.tar.gz).  

This pipeline searches for continuous gravitational wave signals in time-domain data using the F-statistic on data from a network of detectors (currently the available detectors are H1 (Hanford), L1 (Livingston) of LIGO and V1 (Cascina) of Virgo)  

