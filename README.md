# Data-Driven-Estimation
Neural model state and parameter estimation from data

Copyright (C) 2018 Philippa J. Karoly, Dean R. Freestone - All Rights Reserved

You may use, distribute and modify this code under the terms of the GNU General Public license v2.0. You can find a copy of this license in file LICENSE.md. Note than when distributing derived works, the source code of the work must be made available under the same license.

Figures in this repository are licensed under Creative Commons with conditions:
- Attribution: please cite Seizure Pathways: A model-based investigation, P.J. Karoly, L. Kuhlmann, D. Soudry, D.B. Grayden, M.J. Cook and D.R. Freestone (2018) *PLoS Comp. Bio.* \[under review\]
- Share-Alike: you must distribute derived works under the same license
- Non-Commercial

To cite this code in academic research please reference:
Seizure Pathways: A model-based investigation, P.J. Karoly, L. Kuhlmann, D. Soudry, D.B. Grayden, M.J. Cook and D.R. Freestone (2018) *PLoS Comp. Bio.* \[under review\]

## Code
Source code is located in src folder
- **set_params.m**: A function that sets the parameters of the neural mass model (Jansen & Rit 1995)
- **g.m**: A function that computes the erf sigmoid
- **prop_mean_covariance.m**: A function that computes the posterior mean and covariance of the state/parameter distribution propagated through the neural mass model

Example code is located in Example folder
- **Example.m**: Runs the estimation for a single seizure
- **generateData.m**: Runs the model at different values of input and generates simulated data

## Data
One example seizure is provided in data folder. Additional seizure data can be found [here (NOTE DATA TO BE UPLOADED AFTER PUBLICATION)](https://figshare.com/s/733ed9c6b3f71875410a)
Data is a single .mat file with a variable, *Seizure* that has dimension T x N, where T is the number of samples, and N is the number of electrode channels (16). The seizure onset is 5 minutes from the start of the data and seizure offset is 1 minute from the end of the data. Data is sampled at 400Hz.

## Supplementary Figures
These figures relate specifically to the results presented in Karoly et al (2018). We provide additional figures for connectivity parameter estimation and signal energy showing all 16 channels.

## Notes on Filter Implementation
If you are not familiar with Kalman filtering, a review of the recommended resources (or similar) is strongly advised before implementing this code. Density filters can be plagued by numerical instability and in practise fine-tuning of filter parameters is often required to run the filter on real-world data. We provide a non-exhaustive list of some gotchas and heuristics that we have observed over the years.

- Data range is inconsistent with the model (scale your data so it lies not to far from the bounds of what your model can simulate)
- State/parameter values differ by many orders of magnitude (scale your model)
- Mismatched DC between model and data (add one constant offset parameter to the estimate)
- Covariance, P becomes asymmetric (enforce symmetry with P = (P + P') / 2) 
- Covariance becomes too small, breaking the numerics of the filter (add a very small amount to the diagonal of Q, initialise P larger)
- Estimates don't converge (try using annealing to gradually increase R)
- Too many parameters to estimate (maybe you need to lump them together)

## Further References
More information can be found in the following refrences
1. [Freestone, D. R., Karoly, P. J., Neši?, D., Aram, P., Cook, M. J., & Grayden, D. B. (2014). 
Estimation of effective connectivity via data-driven neural modeling. Frontiers in neuroscience, 8, 383](https://www.frontiersin.org/articles/10.3389/fnins.2014.00383/full)

2. [Ahmadizadeh, S., Karoly, P. J., Nešic, D., Grayden, D. B., Cook, M. J., Soudry, D., & Freestone, D. R. (2018). 
Bifurcation analysis of two coupled Jansen-Rit neural mass models. PloS one, 13(3)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0192842)

3. [Bifurcation-Estimation (Repository)](https://github.com/pkaroly/Bifurcation-Estimation)
 
4. Kuhlmann, L., Freestone, D. R., Manton, J. H., Heyse, B., Vereecke, H. E., Lipping, T., ... & Liley, D. T. (2016). 
Neural mass model-based tracking of anesthetic brain states. NeuroImage, 133, 438-456.

## Recommended Resources

1. Optimal State Estimation: Kalman, H Infinity, and Nonlinear Approaches, Dan Simon and [the authors example code!](http://academic.csuohio.edu/simond/)

3. Neural Control Engineering: The Emerging Intersection Between Control Theory, Steve Schiff (esp Chapter 2 and Chapter 5) and [all the example code](https://www.dropbox.com/sh/b23je0226el37wx/AABQJlWFxiI36u_cJba33NeXa?dl=0&preview=Code+Archives+Neural+Control+Engineering+062512.zip)
