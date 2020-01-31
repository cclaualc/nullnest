![](/logo.png) 


# Overview

The `nullnest` repository provides a set of programs developed for constructing a **null model** of bipartite networks which constraints the degree sequences -on average-, along with the measurement of the degree of **nestedness** of bipartite networks using a large variety of metrics allowing for the comparison of this value against the null expectation. The programs included in this package were used in the works [Payrato2019](#references) and [preprint-Payrato2019](#references), where the interested reader will find more information about the characteristics, functioning and scope of both the null model and each nestedness metrics. 

The package is divided into two main blocks. On the one hand, we provide a [program](#simulated-annealing) to compute the null model, for any bipartite network introduced by the user, which keeps the original degree sequence constant on average while maximizing the entropy of the null ensemble. We also give the ready-to-use results of the [null model](#probability-matrices) for an important number of real networks available online. On the other hand, the package contains the programs to measure the degree of nestedness of any network, measuring as well the first two moments of its null distribution, either by using [analytical](#analytical-nestedness-measures) expressions or by numerically [sampling](#sampled-nestedness-measures) the null ensemble. 

When using the `nullnest` repository please acknowledge the authors by citing [Payrato2019](#references) and/or [preprint-Payrato2019](#references), together with the authors of any external package or library employed in a particular program (for instance, when calling `R` packages or the random number generator `dranxor.f90` by Toral and Chakrabarti ([Toral1993](#references))). 

<br/><br/>

# Index


This documentation is organized as follows:


- [Examples](#examples) 
- [Null Model](#null-model)
	- [Simulated Annealing](#simulated-annealing)
	- [Probability Matrices](#probability-matrices)
- [Analytical Nestedness Measures](#analytical-nestedness-measures)
	- [Analytical measures with NODF](#analytical-measures-with-NODF)
	- [Analytical measures with the Spectral Radius](#analytical-measures-with-the-spectral-radius)
- [Sampled Nestedness Measures](#sampled-nestedness-measures)
	- [Temperature](#temperature)
	- [Manhattan Distance (NMD)](#manhattan-distance)
	- [NODF](#nodf)
	- [Discrepancy](#discrepancy)
	- [NIR](#NIR)
	- [Spectral Radius](#spectral-radius)

Each of the sections corresponds to one of the folders of the repository and provides information about its content, format and functioning.

<br/><br/>

# Examples 

Several programs of the package share the same type of input/output files. We provide below the details of the content and format of these files. In all three cases, the index *i* stands for the index of the system studied, to be specified by the user. In this folder, we illustrate the format with a few examples files indexed by *i=0*.

| File name | Description |
| ---------- | ---------- |
|  **matrix*i*ord.txt** | It contains the real bi-adjacency  *n x m* interaction matrix, for which we intend to construct the particular null ensemble. Importantly, both rows and columns are ordered by decreasing degree, and interactions are binary (1's represent presence of interactions and 0's absences). Species with no interactions (zero degree) must be removed from the matrix before being used. Moreover, the matrix elements should be separated by spaces, without commas. |
| **matrix*i*rand.txt** | It contains an *n x m* bi-adjacency matrix, representing the expected probabilities of interaction between nodes in the null ensemble. That is, each matrix element corresponds to the probability *p<sub>i,j</sub>* that node *i* from the row's guild interacts with node *j* from the column's guild. Matrix elements are separated by spaces, without commas. In section [Probability Matrices](#probability-matrices) we provide a large repository containing the probability matrices that allow to generate the random ensembles corresponding to 231 real networks. |
| **general*i*.txt** | It contains the dimension of the bi-adjacency interaction matrices above: the first line contains *n*, the number of rows, while the second line contains *m*, the number of columns. |

<br/><br/>

# Null Model

## Simulated Annealing

This folder contains a `FORTRAN90` program named `simulated_annealing.f90`. This program implements a null model based on the maximum entropy and maximum likelihood ensemble proposed by Squartini and Garlaschelli [(Squartini2011](#references)), and later extended to bipartite economic networks by Saracco et al. ([Saracco2015](#references)). Specifically, we require the degree sequences to be kept only in average, while the probability of finding a given network in the statistical ensemble is constructed using the Exponential random graph model. We follow the numerical implementation made in [Payrato2019](#references), where we applied this null model to the study of nestedness in ecological networks. For more details on the aim and construction of this null model we refer the interested reader to the summarized review in [Payrato2019](#references) or the extensive discussion in [Squartini2011](#references). 

The program `simulated_annealing.f90` constructs the mentioned null model for a given bipartite network provided by the user, typically representing a real system. In particular, the program performs a pseudo-random search using simulated annealing in order to obtain the global maximum of the probability of finding the real degree sequences in the ensemble. This probability is determined by the Lagrange multipliers (in principle one per node), which eventually define the corresponding null ensemble for the specific network under study. We adapt to our purpose and particular problem the code for simulated annealing developed by Goffe et al. ([Goffe1994](#references)), which is in turn an implementation of the algorithm proposed by Corana et al. ([Corana1987](#references)). The random number generator employed is the one by Toral and Chakrabarti ([Toral1993](#references)), included in the repository and named `dranxor.f90`.

Importantly, and compared to other implementations of the null model that can be found in the literature (see for example the [bicm](#https://github.com/tsakim/bicm) module), the `simulated_annealing.f90` program is specially concerned with finding the *global* maximum and hence avoiding getting trapped in local maxima. This means that it can be more computational demanding than alternative packages, but the optimization method is highly robust.

In very general terms, the code we provide performs the following **tasks**:

- First, it tackles the nodes with identical degree. This allows converting the actual optimization problem to an equivalent problem with smaller dimension, that this, having fewer variables to optimize. This reduction of dimensionality boosts both the accuracy and the efficiency of the global search.
- Second, it performs a global optimization of the log-likelihood of finding the real degree sequences in the ensemble, by means of a pseudo-random search. The algorithm is characterized by accepting uphill moves as well as downhill moves following the Metropolis criteria. This non-deterministic approach performs an exhaustive exploration of the search space and avoids getting trapped in local maximums. As a parameter named *temperature* decreases, the probability of accepting downhill moves declines and eventually the algorithm converges to the global maximum.
- Finally, using the optimal Lagrange multipliers, it constructs the matrix containing the probability of interactions between nodes of different guilds. This matrix and other outputs are produced as we specify below.


The following subsections aim at providing a practical guide for easily using the program. Though not indispensable for running the program, we highly recommend to the user to browse through the commented code to obtain a more detailed insight on how the program works. When using this program, please acknowledge the authors by citing [Payrato2019](#references), [Goffe1994](#references) and [Toral1993](#references).



### Input

The program takes two input files characterizing the real network we wish to randomize: the file **matrix*i*ord.txt**, containing its corresponding bi-adjacency matrix, and the file **general*i*.txt**, containing its dimension. The format of both files is detailed in section [Examples](#examples) of this documentation.


Furthermore, the program also accepts two arguments by command line: the integer index *i* identifying the real network for which we construct the null model, and the integer index *iloop* which identifies the particular run and permits performing different, independent runs. In the subsections below we specify how to introduce these arguments. 

### Output

As output, the program produces three files: 

| File name | Description |
| ---------- | ---------- |
|**matrix*i*rand.txt** | It contains the *n x m* bi-adjacency matrix containing the expected probabilities of interaction between nodes in the null ensemble, calculated using the optimal Lagrange multipliers. Its format is specified in section [Examples](#examples). |
| **max\_loglikelihood*i*.txt** | It contains the index *iloop* to identify the run, followed by the maximum of the log-likelihood found by the optimization algorithm. |
| **lagrange\_multipliers*i*.txt** | It contains the optimal Lagrange multipliers found by the simulated annealing. |


### Parameters and options

The code is prepared to run properly without need of adjusting any parameters. Indeed, the parameters concerning the simulated annealing search are set following the recommendations of Corana et al. ([Corana1987](#references)) and Goffe et al. ([Goffe1994](#references)), and, at the same time, aiming at a reliable but efficient calculation. However, a user might tune some of the parameters freely by checking the indications given by Goffe et al. ([Goffe1994](#references)). In order to facilitate this, we kept whenever possible the same notation as Goffe et al.

We next specify a number of options that might be set by commenting or de-commenting some sections of the code. These options are not active by default, so if needed they have to be switched on by the user. Although the details of each option are to be found in the code, we briefly enumerate them below: 

| Option | Description |
| ---------- | ---------- |
| Print the evolution of the simulated annealing | It is possible to print the temperature and global maximum found at each iteration in order to check by terminal the state of the optimization procedure while the program is running. |
| Print the accuracy of the solution | After converging to a solution, the program calculates the average degree sequences in the ensemble. It is possible to print the differences between this expected degree sequences and the real ones by terminal in order to check the results. |
| Change the output format to avoid over-writing | As mentioned above, it is possible to perform independent runs of the whole program by changing the argument *iloop*. By default, the program overwrites previous output files at each run. However, this can be changed so that new runs are appended to an existing file in order to eventually compare the different outputs. |




### Compiling and running the program

The `simulated_annealing.f90` program needs to be compiled before being executed. Using the open-source compiler `gfortran`, the program can be compiled as follows:

```bash
gfortran -Og -o "simulated_annealing.exe" simulated_annealing.f90
```

When executing the program, we must provide two integer arguments in the following form:

```bash
.\simulated_annealing.exe i iloop
```

These indexes *i* and *iloop* correspond, as mentioned before, to the index of the system under study and to the index identifying the run.  We provide a complete example on how to compile and execute the program, using the example files, in the file **bash\_example.sh**.

Given the characteristics of the exhaustive search performed by simulated annealing, it can be quite computationally demanding. To provide a benchmark, the programs takes around 26 minutes to solve a network of 84 *x* 101 species, running in a processor CORE i5.


## Probability Matrices

In this folder we provide, for a large set of real networks, the main results of constructing the null model proposed by Squartini and Garlaschelli ([Squartini2011](#references)) and applied in [Payrato2019](#references) to bipartite ecological networks. In particular, the folder includes the matrix containing the probability of interaction among species in the null ensemble, corresponding to each real network in the dataset. 

These probability matrices are given for 222 ecological networks belonging to the [Web of life](http://www.web-of-life.es/) public database, one ecological network by Burkle et al. ([Burkle2013](#references)) available in [dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.rp321) and 8 economic networks studied in Hernández et al. ([Hernandez2018](#references)) available in [figshare](https://figshare.com/articles/data_used_for_article_zip/6080396). The purpose of this folder is to constitute a repository where the user can rapidly find the ready-to-use probability matrices corresponding to a large group of empirical networks, without need of running the randomizing algorithm.

The folder contains three types of files: 

| File name | Description |
| ---------- | ---------- |
| **references.csv** | It contains a legend to identify the real networks for which the null model is constructed. In detail, it provides the index *i* used to identify each network in our folder (column *Index*), followed by the name of the corresponding database (*Database*), the original name given to the network in that dataset (*ID*), the type of interaction composing the real network, and finally the corresponding reference. |
| **matrix*i*rand.txt** | It contains a *n x m* matrix representing the average probability of interaction between the species of the two guilds, calculated over the maximum-likelihood and maximum-entropy null ensemble. Each matrix corresponds to one of the real networks, as coded in the file *references.csv*, which may be identified by its index *i*. The format of this file is specified in section [Examples](#examples) of this documentation. |
| **general*i*.txt** | It contains the dimension of each bi-adjacency matrix written in matrix*i*rand.txt. Its format is specified in section [Examples](#examples) .|


Please take into account that in the construction of the null model, we take as input an ordered matrix where species with no interactions have been removed. This has to be considered when comparing the present results with the original real matrices in the corresponding database, since some of them occasionally contain zero-degree species.

<br/><br/>

# Analytical Nestedness Measures

This folder contains two programs that perform analytical measures of nestedness using a pair of nestedness metrics: the NODF ([Almeida2008](#references)) and the Spectral Radius ([Staniczenko2013](#references)). In particular, we implement the analytical expressions derived in [Payrato2019](#references) that allow to directly calculate the first two moments of the nestedness ditribution across the ensemble without numerically sampling it.


## Analytic measures with NODF

This folder contains a program named `NODF_analytic.f90` which calculates three quantities using the metrics NODF ([Almeida2008](#references)).


1. The real value of nestedness, measured by NODF, of the real network provided by the user.
2. The average expectation of NODF in the null ensemble corresponding to the real network. This quantity is computed using the analytical expression for the average derived in [Payrato2019](#references), which uses solely the probability matrix of interactions among species.
3. The standard deviation of the distribution of NODF in the null ensemble, derived analytically in [Payrato2019](#references).


### Input

The program takes the following input: the file  **matrix*i*ord.txt**, which contains the real matrix to study; the file  **matrix*i*rand.txt**, which contains the corresponding probability matrix of interaction determined by the null model; and the file  **general*i*.txt**, which provides information about the matrices' size. All the three files have the format specified in section [Examples](#examples) of this documentation. 

### Output

As output, the program produces the following quantities measured by NODF: the nestedness of the real network, the expected average in the ensemble and the standard deviation of the nestedness distribution in the null ensemble. It writes the three quantities in a file called **analytic\_NODF.txt**.

### Compiling and running the program

To compile the program using the compiler `gfortran`, one may use the following command:

    gfortran -Og -o "NODF_analytic.exe" NODF_analytic.f90

To run the compiled program on a network indexed by *i*, simply use the syntax:


    .\NODF_analytic.exe i


An example bash script that performs these actions on the example files (indexed with *i=0*) is included in the folder and named **bash\_NODF\_analytic.sh**.


## Analytic measures with the Spectral Radius

We provide an `R` program named `spectral_radius_analytic.R` that computes the following quantities using the spectral radius ([Staniczenko2013](#references)):

1. The real value of the spectral radius calculated for a matrix of interactions provided by the user.
2. The first two moments (average value and standard deviation) of the distribution of spectral radius over the null ensemble corresponding to the given network. This is calculated using solely the probability matrix of interactions with the theoretical expression derived in [Payrato2019](#references).


The program requires to have installed `R` ([RCore](#references)) as well as the following `R` packages: the `MASS` package ([massRpackage](#references)) and the `RSpectra` package ([RSpectrapackage](#references)).

### Input

The program takes the following input, defined as specified in section [Examples](#examples) of this documentation: the file **matrix*i*ord.txt**, which contains the real matrix to study and the file **matrix*i*rand.txt**, which contains the corresponding probability matrix of interaction determined by the null model.

### Output

As output, the program produces the following quantities measured by the spectral radius: the nestedness of the real network, the expected average in the ensemble and the standard deviation of the nestedness distribution in the null ensemble. It writes the three quantities in a file called **analytic\_spectral\_radius.txt**.


### Running the program

To run the program on a network indexed by *i*, simply use the syntax:

```bash
Rscript spectral_radius_analytic.R i
```

We provide an example script in the file **bash\_spectral\_radius\_analytic.sh**.

<br/><br/>

# Sampled Nestedness Measures

In this folder we provide the programs to perform nestedness measures using the six metrics studied in [preprintPayrato2019](#references). In particular, we give the codes to calculate the following nestedness indexes: temperature, nestedness metrics based on Manhattan Distance (NMD), NODF, Discrepancy, Nestedness Index based on network Robustness (NIR) and Spectral Radius. 

Each program can be used to calculate the following quantities, using the corresponding metrics:

1. The real value of the nestedness index, calculated for a matrix of interactions provided by the user.
2. The first two moments (average value and standard deviation) of the distribution of nestedness over the null ensemble corresponding to the given network. 


In order to calculate the quantities in *ii*, the program performs a numerical sampling of the null ensemble. In detail, each program generates a number *Nsample* of null networks (decided by the user) by sampling the probability of interaction between nodes in the ensemble, and then computes the nestedness index on each null network as well as the resulting average and standard deviation. 

Given that some of the programs may be computationally demanding, we recommend to perform a trial run using a small sample size (for example, *Nsample*=10²) to monitor the computational time, then run the program definitely using the wished sample size. We recommend a minimum sampling size of 10³ networks, while an interesting balance between accuracy and computational speed is *Nsample*=10⁴.

Moreover, note that for very small networks the sampling procedure may be problematic. In fact, the sampling of the probability matrix of interactions may generate a null network where no node has interactions. We have observed this phenomena for networks having *<20* nodes, which therefore we recommend should be analysed with special care. 


## Temperature

The program `temperature.R` measures the nestedness quantities detailed above using the temperature metrics ([Atmar1993](#references)), following the implementation proposed by Rodríguez and Santamaria ([Rodriguez2006](#references)). 
We invert the metrics so that a larger value of the index means more nested. The program requires to have installed `R` ([RCore](#references)) as well as the following `R` packages: the `MASS` package ([massRpackage](#references)) and the `bipartite` package ([BipartitePackage](#references)). 

### Input

The program takes the following input, defined as specified in section [Examples](#examples) of this documentation: the file **matrix*i*ord.txt**, which contains the real matrix to study, and the file **matrix*i*rand.txt**, which contains the probability matrix of interactions given by the null model.

### Output

As output, the program produces the following quantities measured by the temperature metrics: the nestedness of the real network, the expected average in the ensemble and the standard deviation of the nestedness distribution in the null ensemble. It writes these three quantities in a file called **temperature.txt**.


### Running the program

To run the program on a network indexed by *i* and generate a sample size of a total of *Nsample* null networks, use the following command:

    Rscript temperature.R i Nsample

We provide an example script in the file **bash\_temperature.sh**. 

## Manhattan Distance

The program `NMD.R` measures the nestedness quantities specified before using the nestedness idnex based on Manhattan distance (NMD from now on) ([Corso2008](#references)). We invert the metrics so that a larger value of the index means more nested, as well as rescale it so that it varies between 0 and 100. The program requires to have installed `R`([RCore](#references)] as well as the following `R` packages: the `MASS` package ([massRpackage](#references)) and the `bipartite` package ([BipartitePackage](#references)). Importantly, the `nestedness.corso` function from `bipartite` used to calculate the NMD is no longer included in the recent versions of the package, which is why to run the program it is required to have installed the bipartite version number 0.90 from its [archive of old sources](https://cran.r-project.org/src/contrib/Archive/bipartite/). 

Furthermore, the nestedness metrics NMD has the peculiarity of involving the comparison with a null model which keeps constant the size and the fill of the network. Hence, besides the size of the sampling (*Nsample*) the user also needs to decide how many random networks *Nreps* should be produced -either for the real calculation or the sample- in order to compute the NMD. As we explicite below, this information is introduced by terminal when running the program.

### Input

The program takes the following input, defined as specified in section [Examples](#examples)  of this documentation: the file **matrix*i*ord.txt**, which contains the real matrix to study, and the file **matrix*i*rand.txt**, which contains the probability matrix of interactions given by the null model.

### Output

As output, the program produces the following quantities measured by the NMD metrics: the nestedness of the real network, the expected average in the ensemble and the standard deviation of the nestedness distribution in the null ensemble. It writes these three quantities in a file called **NMD.txt**.


### Running the program

To run the program on a network indexed by *i*, generate a sample size of a total of *Nsample* null networks and evaluate the metrics with *Nreps* randomized networks, use the following command:

    Rscript NMD.R i Nsample Nreps


We provide an example script in the file **bash\_NMD.sh**. 


## NODF

The `FORTRAN90` program `NODF.f90` measures the nestedness quantities specified above, using the NODF metrics ([Almeida2008](#references)).
The program calls an external random number generator named `dranxor.f90`by Toral and Chakrabarti ([Toral1993](#references)), which is included in the folder.


### Input

The program takes the following input, defined as specified in section [Examples](#examples)  of this documentation: the file **matrix*i*ord.txt**, which contains the real matrix to study, the file **matrix*i*rand.txt**, which contains the probability matrix of interactions given by the null model, and the file **general*i*.txt**, which contains the bi-adjacency matrix's dimension.

### Output

As output, the program produces the following quantities measured by the NODF metrics: the nestedness of the real network, the expected average in the ensemble and the standard deviation of the nestedness distribution in the null ensemble. It writes these three quantities in a file called **NODF.txt**.


### Running the program


To compile the program using the compiler `gfortran`, we recommend using the following command:

    gfortran -Og -o "NODF.exe" NODF.f90

To run the program on a network indexed by *i* and generate a sample size of a total of *Nsample* null networks, use the following command:

    .\NODF.exe i Nsample

We provide an example script in the file **bash\_NODF.sh**.


## Discrepancy

The program `discrepancy.R` measures the nestedness quantities detailed above using the discrepancy metrics ([Brualdi1999](#references)), in a standarized version which implies normalizing by the total number of links and rescaling it between 0 and 100. We also invert the metrics so that the larger the index the more nested the network.

 The program requires to have installed `R` ([RCore](#references)) as well as the following `R` packages: the `MASS` package ([massRpackage](#references)) and the `bipartite` package ([BipartitePackage](#references)). 

### Input

The program takes the following input, defined as specified in section [Examples](#examples)  of this documentation: the file **matrix*i*ord.txt**, which contains the real matrix to study, and the file **matrix*i*rand.txt**, which contains the probability matrix of interactions given by the null model.

### Output

As output, the program produces the following quantities measured by the discrepancy metrics: the nestedness of the real network, the expected average in the ensemble and the standard deviation of the nestedness distribution in the null ensemble. It writes these three quantities in a file called **discrepancy.txt**.


### Running the program

To run the program on a network indexed by *i* and generate a sample containing *Nsample* null networks, use the following command:

    Rscript discrepancy.R i Nsample


We provide an example script in the file **bash\_discrepancy.sh**. 


## NIR

The `FORTRAN90` program `NIR.f90` measures the nestedness quantities specified above, using the *nestedness index based on robustness* (NIR) ([Burgos2009](#references)).
The program calls an external random number generator named `dranxor.f90`by Toral and Chakrabarti ([Toral1993](#references)), which is included in the folder.


### Input

The program takes the following input, defined as specified in section [Examples](#examples)  of this documentation: the file **matrix*i*ord.txt**, which contains the real matrix to study, the file **matrix*i*rand.txt**, which contains the probability matrix of interactions given by the null model, and the file **general*i*.txt**, which contains the bi-adjacency matrix's dimension.

### Output

As output, the program produces the following quantities measured by the NIR metrics: the nestedness of the real network, the expected average and the standard deviation of the nestedness distribution in the null ensemble. It writes these three quantities in a file called **NIR.txt**.


### Running the program


To compile the program using the compiler `gfortran`, we recommend using the following command:

    gfortran -Og -o "NIR.exe" NIR.f90

To run the program on a network indexed by *i* and generate a sample size of a total of *Nsample* null networks, use the following command:

    .\NIR.exe i Nsample

We provide an example script in the file **bash\_NIR.sh**.


## Spectral Radius

The program `spectral_radius.R` measures the nestedness quantities detailed before using the spectral radius metrics ([Staniczenko2013](#references)). The program requires to have installed `R` ([RCore](#references)) as well as the following `R` packages: the `MASS` package ([massRpackage](#references)), the `RSpectra` package ([RSpectra](#references)) and the `bipartite` package ([BipartitePackage](#references)). 


### Input

The program takes the following input, defined as specified in section [Examples](#examples) of this documentation: the file **matrix*i*ord.txt**, which contains the real matrix to study, and the file **matrix*i*rand.txt**, which contains the matrix of interaction probabilities given by the null model.

### Output

As output, the program produces the following quantities measured by the spectral radius: the nestedness of the real network, the expected average and the standard deviation of the nestedness distribution in the null ensemble. It writes these three quantities in a file called **spectral\_radius.txt**.


### Running the program

To run the program on a network indexed by *i* and generate a sampling of *Nsample* null networks, use the following command:

    Rscript spectral_radius.R i Nsample

We provide an example script in the file **bash\_spectral\_radius.sh**. 

<br/><br/>

# Authors

Codes by Clàudia Payrató Borràs (you may find me [here](https://www.researchgate.net/profile/Claudia_Payrato_Borras)). For any doubt, suggestion or correction on the code, you can contact me by mail at the following address: claudiapb13@gmail.com

<br/><br/>

# License


To be selected yet.

<br/><br/>

# References

[[Almeida2008](https://doi.org/10.1111/j.0030-1299.2008.16644.x)] M. Almeida-Neto, P. Guimarães, P. R. Guimarães, R. D. Loyola, and W. Ulrich, “A consistent metric for nestedness analysis in ecological systems: reconciling concept and measurement,” Oikos, vol. 117, no. 8, pp. 1227–1239, 2008.

[[Atmar1993](https://doi.org/10.1007/BF00317508)] W. Atmar and B. D. Patterson, “The measure of order and disorder in the distribution of species in fragmented habitat,” Oecologia, vol. 96, no. 3, pp. 373–382, 1993.

[[Brualdi1999](https://doi.org/10.1007/s004420050784)] R. A. Brualdi and J. G. Sanderson, “Nested species subsets, gaps, and discrepancy,” Oecologia, vol. 119, no. 2, pp. 256–264, 1999.

[[Burgos2009](https://doi.org/10.1016/j.cpc.2008.11.007)] E. Burgos, H. Ceva, L. Hernández, and R. Perazzo, "Understanding and characterizing nestedness in mutualistic bipartite networks", Computer Physics Communications , vol. 180, no. 4, pp. 532-535, 2009.

[[Burkle2013](https://science.sciencemag.org/content/339/6127/1611)] Burkle, L. A., Marlin, J. C., & Knight, T. M. (2013). Plant-pollinator interactions over 120 years: loss of species, co-occurrence, and function. Science, 339(6127), 1611-1615.

[[Corana1987](https://doi.org/10.1145/29380.29864)] A. Corana, M. Marchesi, C. Martini, and S. Ridella, “Minimizing multimodal functions of continuous variables with the “simulated annealing” algorithm corrigenda for this article is available here,” ACM Transactions on Mathematical Software (TOMS), vol. 13, no. 3, pp. 262–280, 1987.

[[Corso2008](arXiv:0803.0007)] G. Corso, A. I. Araujo, and A. M. Almeida, “A new nestedness estimator in community networks,” arXiv preprint arXiv:0803.0007, 2008.

[[Dormann2008](http://www.biom.uni-freiburg.de/Dateien/PDF/dormann2008rnews.pdf)] C. F. Dormann, B. Gruber, and J. Fruend, “Introducing the bipartite package: Analysing ecological networks.,” R News, vol. 8, no. 2, pp. 8–11, 2008.

[[Goffe1994](https://doi.org/10.1016/0304-4076(94)90038-8)] W. L. Goffe, G. D. Ferrier, and J. Rogers, “Global optimization of statistical functions with simulated annealing,” Journal of econometrics, vol. 60, no. 1-2, pp. 65–99, 1994.

[[Hernandez2018](https://doi.org/10.1371/journal.pone.0196206)] L. Hernández, A. Vignes, and S. Saba, "Trust or robustness? an ecological approach to the study of auction and bilateral markets", PloS one , vol. 13, no. 5, p. e0196206, 2018

[[Payrato2019](https://doi.org/10.1103/PhysRevX.9.031024)] C. Payrató-Borràs, L. Hernández, and Y. Moreno, “Breaking the spell of nestedness: The entropic origin of nestedness in mutualistic systems,” Physical Review X, vol. 9, no. 3, p. 031024, 2019.

[preprint-Payrato2019] "Lost in nestedness? ..." To be included.

[[rARPACK](https://CRAN.R-project.org/package=rARPACK)] rARPACK: Solvers for Large Scale Eigenvalue and SVD Problems. 

[[RCore](https://cran.r-project.org/package=CORE)] R Core Team, R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria, 2013.

[[Rodriguez2006](https://doi.org/10.1111/j.1365-2699.2006.01444.x)] M. A. Rodrı́guez-Gironés and L. Santamarı́a, “A new algorithm to calculate the nestedness temperature of presence–absence matrices,” Journal of Biogeography, vol. 33, no. 5, pp. 924–935, 2006.

[[RSpectra](https://CRAN.R-project.org/package=RSpectra)] G. G. J. N. Yixuan Qiu, Jiali Mei, RSpectra: Solvers for Large-Scale Eigenvalue and SVD Problems.

[[Saracco2015](https://doi.org/10.1038/srep10595)] F. Saracco, R. Di Clemente, A. Gabrielli, and T. Squartini, “Randomizing bipartite networks: the case of the world trade web,” Scientific Reports, vol. 5, no. 10595, 2015.

[[Squartini2011](https://doi.org/10.1088/1367-2630/13/8/083001)] T. Squartini and D. Garlaschelli, “Analytical maximum-likelihood method to detect patterns in real networks,” New Journal of Physics, vol. 13, no. 8, p. 083001, 2011.

[[Staniczenko2013](https://doi.org/10.1038/ncomms2422)] P. P. Staniczenko, J. C. Kopp, and S. Allesina, “The ghost of nestedness in ecological networks,” Nature communications, vol. 4, p. 1391, 2013.

[[Toral1993](https://doi.org/10.1016/0010-4655(93)90016-6)] R. Toral and A. Chakrabarti, “Generation of gaussian distributed random numbers by using a numerical inversion method,” Computer physics communications, vol. 74, no. 3, pp. 327–334, 1993.

