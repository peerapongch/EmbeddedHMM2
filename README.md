# EmbeddedHMM2

This repository contains the code relevant to my master's project on inference in high-dimensional non-linear non-Gaussian state space models. This project is co-supervised by Alexander Shestopaloff and Arnaud Doucet.

Brief description: State space models are an important class of models with wide use in a variety of disciplines such as biology, econometrics and signal processing. This project will consider Bayesian inference in state space models where the state can be high-dimensional using Markov Chain Monte Carlo (MCMC) methods that have recently been developed in Shestopaloff and Neal (2017). The student working on the project will first reproduce the results of the methods presented in this paper, implementing the samplers for latent states, and then perform an extension of the results to do Bayesian inference for both parameters and latent states by adding a sampling step for the parameters, which are currently taken to be known.

The code focuses on two experimental models described in the dissertation. Model 1 is a unimodal model with a Gaussian latent process with a Poisson observation density at each time. Model 2 is a bimodal model for which the Poisson observations depend on the absolute value of the latent states at the corresponding time. The scripts relevant to sampling and analysis of these models are located in **cases/model1/** and **cases/model2/** respectively. Note that each model has two sample instances used in the project, a smaller 3-dimensional variant used for testing and the full 10-dimensional variant described in the dissertation. The smaller variants for Model 1 and 2 are **Poisson1** and **Poisson3** respectively. The full Models 1 and 2 are **Poisson2** and **Poisson4** respectively. Each model directory also contains data sufficient to reproduce the results, as well as the scripts used in plotting figures shown in the dissertation. 

Notably, for both models, scripts with the keyword **long_run** are used to execute long sampler runs to generate a large sample of latent state sequences from the posterior smoothing distribution. These scripts can be run straightaway but it would take a long time to complete. It would be more efficient to break up parts of the scripts and run them in parallel, for example on a virtual machine with multiple processors and large memory as it was performed in this project. Using these samples, the autocorrelation time of each sampler and the smoothing mean of the latent sequence can be computed using functions in the **evaluation/** directory. Lastly, and most importantly, when using RStudio to run these scripts, please set the working directory as the source file directory. This can be done by Session --> Set Working Directory --> To Source File Location. 

The samplers are coded up form scratch to compare the effectiveness of three sampling methods for the latent sequence, namely Metropolis-Hastings (here referred simply as Metropolis), Particle Gibbs with Backward Simulation (PGBS), and Embedded Hidden Markov Model (Embedded HMM). For the Embedded HMM sampler, we investigate separately the forward and backward pool state selection scheme, resulting into different samplers, Lambda and Gamma samplers respectively. Each sampler is implemented as a function and separated into different files. In particular, each file is prefixed with the model name followed by an abbreviation of the sampler name. For example, the Metropolis sampler for model 1 is implemented in the file **samplers/model1_met.R**. In addition, investigation is made into combination of these samplers resulting in variants of abovementioned samplers which are described in the dissertation. The abbreviations are summarised below. 

| Sampler Name | Abbreviation |
| ------------- | ------------- |
| Metropolis  | met  |
| PGBS  | pgbs  |
| PGBS + Metropolis| pgmet |
| Lambda | lambda |
| Gamma | gamma |
| Combined | combined |

Futhermore, different proposals for parameter sampling via Metropolis are also implemented alongside the Lambda sampler. All these files are also located in **samplers/**. Again, the samplers are implemented for Model 1 and 2 separately and they can be found in the following files: 
* model1_lambda_param_RW.R : random walk Metropolis for all the model parameters 
* model1_lambda_param_RW_SHSC.R : random walk Metropolis for all the model parameters, followed by shift and scale updates for the parameters of the observation distribution
* model1_lambda_param_RW_SHSC_2.R : random walk Metropolis for the parameters of the transition distribution, followed by shift an scale updates for the parameters of observation distribution
* model2_lambda_param_RW.R : random walk Metropolis for all the model parameters 
* model2_lambda_param_RW_SC.R : random walk Metropolis for all the model parameters, followed by a scale update for the parameters of the observation distribution
* model2_lambda_param_RW_SC_2.R : random walk Metropolis for the parameters of the transition distribution, followed by a scale update for the parameters of observation distribution

## Main reference 
A.Y. Shestopaloff and R. M. Neal (2017). Sampling latent states for high-dimensional nonlinear state space models with the embedded HMM method. Bayesian Analysis.
