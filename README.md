# Data-Driven-Estimation
Neural model state and parameter estimation from data by Karoly et al (2018)

This code is licensed under the MIT license (see LICENSE.md)

Figures in this repository are copyright (2018) Philippa Karoly

To cite this code in academic research please reference:
*Seizure Pathways: A model-based investigation* (2018) P.J. Karoly, L. Kuhlmann, D. Soudry, D.B. Grayden, M.J. Cook and D.R. Freestone \[under review\]

## Code
Source code is located in src folder
- **set_params.m**: A function that sets the parameters of the neural mass model (Jansen & Rit 1995)
- **g.m**: A function that computes the erf sigmoid
- **prop_mean_covariance.m**: A function that computes the posteiror mean and covarnacie of the state/parameter distribution propagated through the neural mass model

Example code is located in Example folder
- **Example.m**
- **generateData.m**

## Data
One example seizure is provided in data folder. Additional seizure data can be found [here](www.google.com).
Data is a single .mat file with a variable, *Seizure* that has dimension T x N, where T is the number of samples, and N is the number of electrode channels (16). The seizure onset is 5 minutes from the start of the data and seizure offset is 1 minute from the end of the data. Data is sampled at 400Hz.

## Supplementary Figures
