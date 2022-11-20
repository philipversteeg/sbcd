# SBCD: Local Constraint-Based Causal Discovery under Selection Bias

Code accompanying the simulation experiments in "Local Constraint Based Causal Discovery under Selection Bias". Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

### Installation

We use R for simulations. The required packages can be installed from CRAN and BioConductor with the following R commands.
```{r}
install.packages(c("pcalg", "InvariantCausalPrediction", "expm", "PRROC"))
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RBGL")
```

### Experiments

The experiments are run by calling the **Simulate.R** script with the appropriate command line arguments.
```shell script
Rscript Simulate.R [exp] [nseed] [nsamples] [ancestors/patterns]
```
Here **nseed** and **nsamples** are the number of random models and the number of samples gathered in each model. The **ancestors/patterns** argument indicates the true condition that is used to compare predictions to, and the **exp** parameter represents with model type are sampled: 
- 1: fixed graph
- 2: random small graphs
- 3: random large graphs

In **run_simulate.sh** we gathered the calls to the **simulate.R** script that are used to produce the figures in the paper. Results and log files are saved in the **results** and **logs** folders respectively. 
For details we refer to the main paper and appendix.

### Citation

The full paper can be found in the Proceedings of Machine Learning Research [here](https://proceedings.mlr.press/v177/versteeg22a.html), with the following citation.
```latex
@inproceedings{versteeg2022local,
  title = {Local Constraint-Based Causal Discovery under Selection Bias},
  author = {Versteeg, Philip and Mooij, Joris and Zhang, Cheng},
  booktitle = {Proceedings of the First Conference on Causal Learning and Reasoning},
  pages = {840--860},
  year = {2022},
  volume = {177},
  series = {Proceedings of Machine Learning Research},
  publisher = {PMLR}
}
```
