_A project of the [Environmental Forecasting Initiative 2023 Unconference](https://ecoforecast.org/efi-rcn-unconference-2023/) (Original issue description [here](https://github.com/eco4cast/unconf-2023/issues/9))_


## Overview

The goal of this project was to outline and demonstrate methods of supporting interpretation of uncertainty in model forecasts.

### What and Why

 - _Uncertainty Attribution_ decomposes uncertainty in model forecasts to different components contributing to it, which can include
     - Parameters (individal or collective)
     - Drivers / Covariates (individal or collective)
     - Initial Conditions
     - Process Error (Structural or Stochastic)
     - Random effects
     - Observation Errors

Uncertainty attribution can be performed for a single forecast point, a group of points, or performed _across the range of a variable_ (for instance, how far ahead the forecast is made) to examine how contribution to uncertainty varies. It is akin to sensitivity analysis to where the sensitivities are to inputs of the model (parameters, drivers, etc.) and their range is drawn from their own uncertainty.

   
-  _Post-hoc uncertainty modeling_ analyzes how uncertainy in forecasts varies across the predictor space. 

In both cases, these analyses apply to _trained_ models, but they can be used as part of feedback to model updating. For instance, uncertainty attribution can identify parameters that contribute heavily to uncertainty and thus would be good targets to better identify. Post-hoc uncertainty modeling can identify conditions of poor performance that more data could improve. Both can be used together - areas of high uncertainty can then be analyzed to attribute the biggest contributors.

Both tools may also be used with _scenario_ or _simulation_ analyses, where forecast models can be re-fit under simulated data to determine how different training data changes the amount, patternm abd contributors to forecast uncertainty.


### How

We performed _uncertainty attribution_ analyses on two forecasts from the EFI Theory Working Groups's [repository of machine-learning forecast models](https://github.com/eco4cast/Forecast_submissions): An ARIMA and a random forest model of [phenology](https://projects.ecoforecast.org/neon4cast-dashboard/phenology), both at the Harvard Forest NEON site.

We use two methods: _one-at-a-time_ and _Sobol_.  One-at-a-time ....

_Sobol_ is a sensitivity analysis method that determines the model sensitivity to each input by running the model with all inputs except one perturbed, holding each input fixed in turn.  


### Repo structure

 - `randfor.R` performs an uncertainty attribution analysis of a random forest phenology forecast model.

 - `ARIMA.R` performs an uncertainty attribution analysis of an ARIMA phenology forecast model.


### Some whiteboard snapshots


Mapping the problem:

![](https://hackmd.io/_uploads/r1HqVHXd3.jpg)

Next steps:

![](https://hackmd.io/_uploads/HJzFdSm_3.jpg)


### References

 - Dietze, M. 2017. Prediction in ecology: a first-principles framework. Ecological Applications 27: 2048–2060. DOI: 10.1002/eap.1589
 - Lewis ASL, CR Rollinson, AJ Allyn, J Ashander, S Brodie, CB Brookson, E Collins, M Dietze, AS Gallinat, N Juvigny-Khenafou, G Koren, D McGlinn, JA Peters, NR Record, CJ Robbins, J Tonkin, GM Wardle. 2023. “The power of forecasts to advance ecological theory” Methods in Ecology and Evolution 14(3):746-756  http://doi.org/10.1111/2041-210X.13955