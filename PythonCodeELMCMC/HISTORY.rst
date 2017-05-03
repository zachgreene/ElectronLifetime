from Zach Greene

##################
2017/4/3
##################
Removed getter defficiency and change in LXe outgassing rate, and put in new impurity change at 1485451100.


##################
2017/4/2
##################
Changed number of parametsr from 14 to 15 by including a new getter defficiency following the January 18, 2017 earthquake.  Kept version at 14.  Later changed from 15 to 16 to include change in LXe outgassing rate.


##################
2017/3/31
##################
Changed number of parameters from 15 to 14 by elimnating OutgassingRateLXeDecreasingLinearAdditionalConsts that occurred around the end of SR0.


##################
2017/3/29
##################
Comparing pax_v6.2.0 with pax_v6.4.2 shows lifetimes from the earlier version are incorrectly low.  Rather than reprocessing all data with newer pax, a conversion factor is found, and added in Tools.py.  Rn lifetimes before 1478000000 (November 1, 2016) and between 1484900000 and  1486100000 (January 18 to February 2, 2017)  are corrected by this.


from Qing Lin

@ 2017-02-13
We need another parameter to address for the outgassing rising during the Rn calibration at Christmas. I would just add another getter low efficiency period.


@ 2017-02-07
We need to address for the larger slope of the "slow" increasing.


@ 2016-12-12
We need to address for the new impurity dropping and the getter defficiency found afterwards


@ 2016-11-15
The code is getting slower and slower with more and more time.
1) Creating global SC values instead of requiring them at each calculation step
2) Model the outgassing level chaning during gas-only as exponential instead of linear


@ 2016-10-21
According to the trend, V4 fitting cannot explain that we see the continuous increasing of electron lifetime after 10-15. One explanation is that actually the outgassing levels of both the gas and liquid volume are decreasing over time. 
In this version, we assume that the outgassing level is decreasing following an exponential trend. So there are two extra parameters for the decreasing of gas/liquid volume outgassing.


@ 2016-10-25
Use linear instead of exponential decreasing for V6 fitting.
Last parameters are then the days taken to decrease outgassing level to 0


@ 2016-10-03
Change V4 to fit with outgassing level change


@ 2016-09-27
V4 is inherited from V3, and is a branch fitting. 
The fitting used the result from V3, but varying the outgassing level after gas-only flow.


@ 2016-09-02
Differences of V3 compare to V2:
1) Abandon the drop at Kr83m flushing. It cannot be explained by O2 trend.


@ 2016-08-08
The V2 code makes the two modifications:
1) Instead of using the very slow and unstable Historian database, the new module switch to using the SC Chris made: 

https://xecluster.lngs.infn.it/dokuwiki/lib/exe/fetch.php?media=xenon:xenon1t:analysis:meetings:my_sc_adventure.html

Also the flow of FCV103 need to be added.

2) Also now the gas manipulation is started, so instead of using the simplified model, for the electron lifetime it is also good to use the full model with gas and liquid concentration separated:

https://xecluster.lngs.infn.it/dokuwiki/doku.php?id=mcfate:refinedmodellingevolutionv1


@ 2016-07-11
It is basically based on the same model used in PythonCodeV3 for electron lifetime trend fitting.

The only difference is that we saw the error matrix obtained with previous bounded fitting (even with a bounded fitting algorithm). Basically in our bounded condition, the gaussian assumption for error does not apply any more. That's the reason we want to switch to MCMC fitting in this code.

The algorithm using in this code is the emcee
