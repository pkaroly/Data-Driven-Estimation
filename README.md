# Data-Driven-Estimation
Neural model state and parameter estimation from data by Karoly et al (2018)

This code is licensed under the MIT license (see LICENSE.md)

Figures in this repository are copyright (2018) Philippa Karoly

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
One example seizure is provided in data folder. Additional seizure data can be found [LINK TO BE ADDED AFTER PUBLICATION](www.google.com).
Data is a single .mat file with a variable, *Seizure* that has dimension T x N, where T is the number of samples, and N is the number of electrode channels (16). The seizure onset is 5 minutes from the start of the data and seizure offset is 1 minute from the end of the data. Data is sampled at 400Hz.

## Supplementary Figures
Additional figures for connectivity parameter estimation and signal energy showing all channels.

## Further References
More information can be found in the following refrences
1. [Freestone, D. R., Karoly, P. J., Neši?, D., Aram, P., Cook, M. J., & Grayden, D. B. (2014). 
Estimation of effective connectivity via data-driven neural modeling. Frontiers in neuroscience, 8, 383](https://www.frontiersin.org/articles/10.3389/fnins.2014.00383/full)

2. [Ahmadizadeh, S., Karoly, P. J., Nešic, D., Grayden, D. B., Cook, M. J., Soudry, D., & Freestone, D. R. (2018). 
Bifurcation analysis of two coupled Jansen-Rit neural mass models. PloS one, 13(3)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0192842)

3. [This repository](https://github.com/pkaroly/Bifurcation-Estimation)
 
4. Kuhlmann, L., Freestone, D. R., Manton, J. H., Heyse, B., Vereecke, H. E., Lipping, T., ... & Liley, D. T. (2016). 
Neural mass model-based tracking of anesthetic brain states. NeuroImage, 133, 438-456.
