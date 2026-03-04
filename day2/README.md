# Day 2 - phylogeny reconstruction

## Introduction

In the theoretical lectures, you have seen how you can reconstruct phylogenies under various approaches: data-based matrices, parsimony, maximum likelihood, and Bayesian. During this practical session, we shall focus on the last two!

## Phylogeny reconstruction under maximum likelihood (`IQ-TREE`)

After inspecting the different alignments we generated during the first practical session, we can see how important it is to maintain correspondence between DNA alignments and AA alignments -- we would always prefer working with a codon-aware alignment!

We will be using `IQ-TREE` for inferring the best-scoring maximum-likelihood tree with two different alignments: the codon-aware alignment `aln_nuc_against_protsuper5.fasta` (**DNA alignment**) and the structure-based alignment `aln_prot_usalign_clean.fasta` (**AA alignment**). If time allows, we will see how to partition the DNA alignment according to codon positions -- if not, all the information is here and you can revise this at home!

Please note that, at this stage, we only have sequence alignments, but we still do not know which model of evolution will better fit these datasets. Do not worry: `IQ-TREE` has a model selection algorithm that can help you find **the best-fitting substitution model** for both your DNA and AA alignments! In addition, you will be able to simultaneously run a bootstrap analysis to calculate bootstrap support values for each of the clades in the inferred phylogeny.

> [!IMPORTANT]
> Remember that bootstrap support values (or bootstrap proportions) do not reflect uncertainty: **they are not confidence intervals**. Bootstrap support values can help you interpret the support for the clades inferred in your phylogeny: clades with higher values (e.g., >60% or >70%) will have stronger support than those with lower values.

Let's run `IQ-TREE`!

```sh
# Run from my_session
# Change directories if you are not
# there yet

# Copy data and scripts for day2
cp -R ~/biol0033-tutorial/day2 .

# Move to directory `day2` and create
# directories for phylogeny inference
cd day2
mkdir -p iq-tree/{dna,aa}

# Let's start with the DNA aln
cd iq-tree/dna
cp ../../../day1/aln/aln_nuc_against_protsuper5.fasta aln_dna.fasta
# Options:
## -s   Path to aln
## -B   Enable bootstraping with 1000 pseudosamples
## -T   Determine the best number of CPU cores to speed up analyses
##
## This command will do the following:
##  - Select the best-fitting nucleotide substitituion model for your data
##  - Reconstruct the phylogeny (best-scoring ML tree)
##  - Run a bootstrap analyses and add bootstrap support values to
##    each clade
iqtree2 -s aln_dna.fasta -B 1000 -T AUTO

# Now, let's run another analysis for the AA alignment
cd ../aa
cp ../../../day1/aln/aln_prot_usalign_clean.fasta aln_aa.fasta
iqtree2 -s aln_aa.fasta -B 1000 -T AUTO
```

The log output files are self explanatory and consist of the following:

* `<aln_file_name>.iqtree`: all your results will be summarised and explained in this file, include a textual representation of the best-scoring ML tree and consensus tree.
* `<aln_file_name>.treefile`: the best-scoring ML tree is saved in Newick format (unrooted). You will see both branch lengths and bootstrap support values printed on the tree.
* `<aln_file_name>.contree`: the consensus tree with both branch lengths and bootstrap support values are saved in Newick format (unrooted).
* `<aln_file_name>.bionj`: starting tree based on an improved neighbour-joining (NJ) tree that `IQ-TREE` uses for ML optimisation. There are no support values, only branch lengths that are not ML-optimised.
* `<aln_file_name>.log`: same content as the screen output.
* `<aln_file_name>.mldist`: ML pairwise distance matrix estimated under the best-fitting substitution model. You may want to use this output file to check whether taxa are extremely divergent, check whether there are any LBA artifacts, or obtain ML-based distances required for running other programs.
* `<aln_file_name>.model.gz`: compressed file with information regarding all tested substitution models and their statistics. E.g.: log-likelihood, number of parameters, AIC/AICc/BIC, estimated model parameters, rate-heterogeneity settings. You should not delete it if you want to re-run or extend the analysis.
* `<aln_file_name>.splits`: bipartition support file used for generating the consensus tree and calculating bootstrap support values for each branch. You should not delete it if you need to restart your analysis.
* `<aln_file_name>.ckp.gz`: if your phylogenetic analysis were interrupted, this is the checkpoint file that you could use to resume it.

> [!IMPORTANT]
> Discuss the results you have obtained for both datasets and how they may complement each other. You may visit the [`IQ-TREE` website](https://iqtree.github.io/doc/Substitution-Models) to learn more about additional nucleotide and amino acid models that have not been discussed during the lectures. You can also navigate the different sections in [their documentation](https://iqtree.github.io/doc/) to get familiar with the settings and other FAQ.

Once you have spent some time discussing the results, you may reveal the following section:

<details>
<summary><b>[ Click here to learn more about the model selection results ]</b></summary>
<br>

### Interpreting model selection results with `IQ-TREE`

#### DNA alignment

The preferred model for the DNA alignment is `TIM3+F+G4`. Let's evaluate each part:

* `TIM3` is a **restricted GTR-like model** that allows for different transition/transversion rates. While you may think this is similar to the HKY85 model, the TIM3 model differs because transitions (i.e., A↔G and C↔T) are different from each other. In other words, the number of A↔G transitions is different from the number of C↔T transitions. Nevertheless, the number of transversions is equal. Consequently, instead of having one parameter $\kappa$, there are three rate parameters: one for A↔G, another for C↔T, and another for transversions. If we were to think about nested models, the order would be as follows (from simpler to more complex): HKY  ⊂  TIM3  ⊂  GTR. HKY is a special case of TIM3 and TIM3 is a restricted version of GTR.
* `F` means that the observed base frequencies of the sequence alignment fit better than assuming equilibrium, reflecting lineage-specific nucleotide composition.
* `G4` relates to the Gamma distribution with 4 categories, thus assumes site-rate heterogeneity. Each category can accommodate the rate variability across a dynamic GTPase (i.e., highly conserved sites will evolve very slowly, while loop regions and variable surface residues will evolve much faster).

Most protein-coding genes in vertebrates tend to present different rates for A↔G and C↔T due to CpG deamination, codon structure, and GC-biased gene conversion (strong in birds, which are actually present in our dataset!). Transversions, however, tend to be responsible for a change in a given AA, thus not occurring that regularly (i.e., nonsynonymous). Consequently, transversion rates my be low and similar to each other. These processes are better fitted by TIM3 than HKY or GTR. In particular, **TIM-type models are selected when genes are under purifying selection, taxa moderately diverged, and codon structure dominates substitution patterns**.

#### AA alignment

The preferred model for the AA alignment is `Q.plant+I+R3`. Let's evaluate each part:

* `Q.plant`: this model is not linked to plant datasets. Instead, the AA substitution matrix designed for this model ([Minh et al. 2021](https://doi.org/10.1093/sysbio/syab010)) was inferred from large sets of plant nuclear proteins, many of which are structurally constrained (e.g., suppression of radical changes; preference for conservative substitutions such as I↔L, D↔E, K↔R; etc.). In other words, **`Q.plant` captures fold-driven evolution, not taxonomy**, which is exactly what fits a folded, dynamic GTPase domain with extremely conserved core residues such as the one we are analysing! A structure-based alignment tends to focus on better aligning secondary structures such as helices and sheets. Nevertheless, spurious homology in loops will be removed. To this end, there may be a strong signal of structural constraint and fewer apparent radical AA changes.  This is the main reason why this matrix outperforms other general models.
* `I`: it seems that there are many positions that remain invariant across the alignment. `IQ-TREE` models these positions as a truly invariant class and not just slow rates.
* `R3`: instead of using a Gamma distribution to model rate evolution, `RX` is related to the free-rate model ([Yang 1995](http://www.genetics.org/content/139/2/993.abstract)). In our case, it looks like a total of 3 categories for classifying rate variation seem to be a better bit than Gamma-distributed rates (e.g., "core", "flexible", "very flexible").

#### Summary of model selection for both datasets

| Aspect              | Nucleotide alignment       | Structure-based protein alignment |
| ------------------- | -------------------------- | --------------------------------- |
| Dominant constraint | Codons & GC bias           | Protein fold & mechanics          |
| Best-fitting model  | TIM3+F+G4                  | Q.plant+I+R3                      |
| Rate heterogeneity  | Smooth, continuous         | Discrete classes                  |
| Invariant sites     | Implicit (slow gamma)      | Explicit (+I)                     |
| Evolutionary signal | Mixed functional + neutral | Mostly functional                 |

#### How should we interpret these results?

* The MX dynamin like GTPase 1 (MX1) protein-coding gene seems to be a conserved catalytic core.
* There seems to be flexible, lineage-specific surface regions.
* Structural constraints seem to be dominating AA evolution, and thus selection may be acting differently from the DNA level.

#### Notes on the `IQ-GPT`

Rob Lanfear has recently developed [a GPT that tries to help users with `IQ-TREE` analyses](https://chatgpt.com/g/g-aZvnPPUW1-iq-gpt). Please feel free to try this AI tool but, as always with AI tools, please proceed with care and be critical of the results obtained. Should you have any queries and/or experience any issue, please check the documentation. If you cannot find the answer to your question there, you may want to check whether your issue has already been posted and discussed in the [`IQ-TREE GitHub`](https://github.com/iqtree/iqtree2/issues). If not, you may want to open a new issue.

</details>

### Comparing tree topologies

Now, we want to test whether the best-scoring ML tree inferred for the AA alignment fits the DNA alignment better than the best-scoring ML tree inferred for said alignment, and vice versa. In order to do this, we need to create a file with the two best-scoring ML trees in Newick format and run AU / SH / KH tests for each sequence alignment. What do theses tests do?

* **Kishino-Hasegawa (KH) test**: this test compares only two trees and tests whether their likelihood difference is significant. This test assumes that the two trees were specified _a priori_. This test will not be valid for us because we are comparing best-scoring ML trees, so we may ignore the result printed on the screen.
* **Shimodaira-Hasegawa (SH) test**: this test compares multiple trees at once and corrects for selection bias. This test assumes that your candidate trees come from the same dataset and trees may have been select _a posteriori_. Our trees are generated with two different datasets, and thus we cannot use this test.
* **Approximately Unbiased (AU) test**: this test compares multiple trees and uses multiscale bootstrap. This test assumes that trees may have been selected _a posteriori_ and that site likelihood are approximately independent. AU is less conservative than SH, but it is the best compromise to balance type I and type II errors.

```sh
# Run from my_session/day2/iq-tree/
# Change directories if you are not
# there yet

# Create file with competing topologies
# Tree1: DNA | Tree2: AA
cat dna/*treefile aa/*treefile > candidates.tree
mkdir -p {aa,dna}/topology_test

# Run first topology test with DNA alignment
cd dna/topology_test
## New options:
## -n 0       Tree searches will not be carried out. Instead,
##            trees passed via option `-z` will be read for 
##            comparison tests and branch length optimisation
## -z         File with candidate tree topologies to compare
## -zb 10000  Number of RELL (Resampling Estimated Log-Likelihoods)
##            bootstrap replicates used to approximate
##            the distribution of likelihood differences
##            RELL is a computational shortcurt to increase
##            bootstrap efficienct -- each RELL replicate is one
##            resampled pseudosample and is used to estimate the 
##            the variance of the likelihood differences
##            The AU test will use the multiscale bootstrap
##            resampling and estimate how often each topology is
##            supported by resampled datasets -- the more replicates
##            the more stable AU p-values
## -pre       Name to give to output files, otherwise it will not 
##            continue as files have been generated for this alignment
iqtree2 -s ../aln_dna.fasta -m TIM3+F+G4 -z ../../candidates.tree -n 0 -zb 10000 -au -pre dna_AUtest

# Do the same for AA data
cd ../../aa/topology_test
iqtree2 -s ../aln_aa.fasta -m Q.plant+I+R3 -z ../../candidates.tree -n 0 -zb 10000 -au -pre aa_AUtest
```

> [!IMPORTANT]
> Discuss the results you have obtained for both datasets and how they may complement each other -- you may want to check the `*.iqtree` files! You may visit the [`IQ-TREE` website](https://iqtree.github.io/doc/Substitution-Models) to learn more about additional nucleotide and amino acid models that have not been discussed during the lectures. You can also navigate the different sections in [their documentation](https://iqtree.github.io/doc/) to get familiar with the settings and other FAQ.

Once you have spent some time discussing the results, you may reveal the following section:

<details>
<summary><b>[ Click here to learn more about the results when comparing tree topologies ]</b></summary>
<br>

Inside the relevant `*.iqtree` files, you should see these tables:

> DNA alignment

```text
Tree      logL    deltaL  bp-RELL    p-KH     p-SH       c-ELW       p-AU
-------------------------------------------------------------------------
  1  -15216.0045       0    0.83 +  0.828 +      1 +     0.818 +    0.843 + 
  2  -15221.3781  5.3736    0.17 +  0.172 +  0.172 +     0.182 +    0.157 + 
```

First, you have the maximum log-likelihood after branch-length optimisation. It looks like tree 1 (DNA) fits the data better than tree 2. The log-likelihood difference is Δℓ ≈ 5 (noticeable, but not necessarily significant). It looks like 83% of the bootstrap replicates were found to support tree 1 (17%, tree 2). We shall ignore columns p-KH and p-SH because these tests are not valid for our dataset. Column c-ELW summraises the expected likelihood weights, so it looks like tree 1 is favoured, but tree 2 still has ~18% of weight. The AU test does not reject any tree (both p-AU > 0.05). Nevertheless, tree 1 is preferred, but tree 2 is statistically compatible with the data.

> AA alignment

```text
Tree      logL    deltaL  bp-RELL    p-KH     p-SH       c-ELW       p-AU
-------------------------------------------------------------------------
  1 -9831.977979  3.9046   0.242 +  0.221 +  0.221 +     0.253 +    0.197 + 
  2  -9828.07342       0   0.758 +  0.779 +      1 +     0.747 +    0.803 + 
```

When comparing the tree topologies for the AA alignment, we see the opposite: tree 2 (generated with the AA alignment) is preferred, but tree 1 is statistically compatible with the data.

</details><br>

> [!TIP]
> These tests can be relevant when comparing tree topologies inferred using the same dataset!

### Partitioning DNA dataset in codon positions

We are going to see how partitioning the dataset in codon positions may affect the inferred tree topology. We will show you how you can specify a partitioning scheme in which first and second codon positions are kept in one alignment block (partition 1) and the third codon positions are kept in a second block (partition 2):

```sh
# Run from my_session/day2/iq-tree/dna
# Change directories if you are not
# there yet

# Create a new directory
mkdir part_CP12CP3
# Copy the NEXUS file that we have already
# prepared
cp ../info_parts/part12_3.nex part_CP12CP3/
##> NOTE: The ".\3" means that the program will read
##> until the end of the alignment in steps of 3

# Run `IQ-TREE`
## New options
## -st DNA   Sequence type is DNA 
## -m MFP    This option will allow different models for
##            each partition
cd part_CP12CP3
iqtree2 -s ../aln_dna.fasta -st DNA -p part12_3.nex -m MFP -B 1000 -pre dna_partCP12CP3

# Compare this tree against the best-scoring ML tree inferred
# for the concatenated sequence
## Tree 1: partitioned | Tree 2: concatenated
cat *treefile ../*treefile > ../../candidates_dna_concVSpart.tree
mkdir compare_partVSconc
cd compare_partVSconc
iqtree2 -s ../../aln_dna.fasta -m TIM3+F+G4 -z ../../../candidates_dna_concVSpart.tree -n 0 -zb 10000 -au -pre dna_conc_AUtest
iqtree2 -s ../../aln_dna.fasta -p ../dna_partCP12CP3.best_model.nex -z ../../../candidates_dna_concVSpart.tree -n 0 -zb 10000 -au -pre dna_part_AUtest
```

> [!IMPORTANT]
> Discuss the results you have obtained for both partitioned and concatenated DNA datasets -- you may want to check the `*iqtree` files! If you have time (or at home), you can also run different partitioning schemes to see their impact on phylogeny inference.

Once you have spent some time discussing the results, you may reveal the following section:

<details>
<summary><b>[ Click here to learn more about the results when comparing tree topologies inferred with a partitioned and a concatenated dataset ]</b></summary>
<br>

Inside the relevant `*iqtree` files, you shall see these tables:

> DNA concatenated alignment

```text
Tree      logL    deltaL  bp-RELL    p-KH     p-SH       c-ELW       p-AU
-------------------------------------------------------------------------
  1 -15216.00452 1.346e-05   0.498 +    0.5 +    0.5 +       0.5 +    0.514 + 
  2  -15216.0045       0   0.501 +    0.5 +      1 +       0.5 +    0.486 + 
```

> DNA partitioned alignment

```text
Tree      logL    deltaL  bp-RELL    p-KH     p-SH       c-ELW       p-AU
-------------------------------------------------------------------------
  1 -14971.32836       0   0.495 +  0.511 +      1 +       0.5 +    0.505 + 
  2  -14971.3284 3.5716e-05   0.505 +  0.489 +  0.489 +       0.5 +    0.495 + 
```

In both cases, the AU test does not reject any tree (both p-AU > 0.05). Both have similar maxim log-likelihood after branch-length optimisation and very close p-values. Nevertheless, tree 1 is narrowly preferred, but tree 2 is still statistically compatible with the data.

</details>

Now, you are ready to learn how to use a Bayesian program for phylogeny inference!

## Phylogeny reconstruction under a Bayesian approach

The second part of this tutorial consists of using `MrBayes` ([Ronquist et al., 2012](https://academic.oup.com/sysbio/article/61/3/539/1674894)) for phylogeny inference under a Bayesian approach. Given that AA models are not implemented in this program, we shall focus on the codon-aware alignment. Based on the best-fitting models found by `IQ-TREE`, the closest implemented nucleotide substitution model in `MrBayes` is the HKY+G4.

> [!IMPORTANT]
> Some phylogenetic software may not have implemented models of evolution available in other programs. Consequently, if the best-fitting model is not available in a program you want to use, you can always run the analysis with the closest implemented model.

We will divide this practical session into the following sections:

* Getting familiar with the commands in `MrBayes`.
* Running `MrBayes` with the codon-aware alignment under the HKY+G4 model.
* Learning how to use `Tracer` for MCMC diagnostics.
* Learning how to use `FigTree` to visualise inferred phylogenies.

Now... Let's get started!

### 1. Getting familiar with the NEXUS format

#### Input data

Firstly, we will be looking at the different commands needed to define a sequence alignment in NEXUS format. In [directory `mrbayes/inp_data`](mrbayes/inp_data), you will find a NEXUS file called [`aln_dna.nex`](inp_data/aln_dna.nex). You can visualise this file here on the GitHub repository or, if you have cloned this repository, you may open this file with your preferred text editor. The main commands to be highlighted are the following:

* `begin data`: this command specifies the beginning of the data block. More details [in the `MrBayes` documentation](https://github.com/NBISweden/MrBayes/blob/develop/doc/manual/Manual_MrBayes_v3.2.pdf).
* `dimensions`: after this command, you define the number of taxa (`ntax`) and the number of characters in your alignment (`nchar`). More details [in the `MrBayes` documentation](https://github.com/NBISweden/MrBayes/blob/develop/doc/manual/Manual_MrBayes_v3.2.pdf).
* `format`: this command specifies the type of data of your alignment (`datatype`), e.g., DNA, AA, RNA, etc. You can also indicate if it is interleaved or not. More details [in the `MrBayes` documentation](https://github.com/NBISweden/MrBayes/blob/develop/doc/manual/Manual_MrBayes_v3.2.pdf).
* `matrix`: after this command (normally next line), you need to include the alignment. Interleave format is accepted as aforementioned. More details [in the `MrBayes` documentation](https://github.com/NBISweden/MrBayes/blob/develop/doc/manual/Manual_MrBayes_v3.2.pdf).
* `end`: this command terminates the data block.

> [!TIP]
> Add a `;` to terminate a command (i.e., `begin data;`). If you do not add a semicolon, `MrBayes` will understand that all the information/arguments you add until the next `;` are part of the last settings block despite being in a new line! In addition, whatever information you add within square brackets is treated as a comment, and thus `MrBayes` does not execute it. Please use the square brackets to add useful comments in your NEXUS files, which will always help you remember what each settings block corresponds to.
> Example of a comment in the NEXUS control file to be read by `MrBayes`: `[ this is a comment and it is not run by MrBayes ]`.

#### `MrBayes` commands

Once we have defined our input alignment in NEXUS format, we can generate the control file that will execute `MrBayes`, which will also be in NEXUS format. You can find a template control file in directory `inp_data` called [ctl_mb_dna.nex](mrbayes/inp_data/ctl_mb_dna.nex). Below, you will find a summary of the commands that we will go through:

> **BLOCK 1: Start `MrBayes` and read input data**

* `begin mrbayes`: this command initiates the block with instructions to run `MrBayes`.
* `log`: this command specifies that you want a log file where all the screen messages printed out by the program will be saved. You can specify your preferred file name (`filename`) and whether you want to append new content or replace the old content with new (`append/replace`). You start or stop this command by including the option `start` or `stop`, respectively. More details [in the `MrBayes` documentation](https://github.com/NBISweden/MrBayes/blob/develop/doc/manual/Manual_MrBayes_v3.2.pdf).
* `execute`: this command is used to let `MrBayes` know that you want to read your input alignment file, which you have previously prepared.

> [!NOTE]
> If you wanted to add one outgroup to root the tree, you could also use command `outgroup` indicating which taxon is to be used as an outgroup.

> **BLOCK 2: Define analysis in `MrBayes`**:

* `charset`: this command is used to specify different datasets that you want to analyse based on your pre-defined input alignment. For instance, you may have two genes which could be labelled with a specific tag and which number of nucleotide bases is defined. This command can also be useful if you decide to partition your dataset and give a specific tag to each alignment block. The format followed is the following: `charset <name_set>=<start_pos>:<end_pos>`.
* `partition`: this command uses the information you have passed to `charset` and the names you gave to each character set. The format followed is the following: `partition <name_partition> = <num_partitions>:<name_charset1>, ...,<name_chraset_n`. E.g.:

  ```text
  # Example 1: partition a concatenated alignment by gene
  # The first 500 nucleotides correspond to the first gene, the
  # other 500 to a second gene
  charset gene1=1-500;
  charset gene2=501-1000;
  partition by-gene=2: gene1, gene2;
  ```

  ```text
  # Example 2: partition the first gene by codon positions 
  # but keep the second gene concatenated
  charset gene1_cp1=1-500\3;
  charset gene1_cp2=2-500\3;
  charset gene1_cp3=3-500\3;
  charset gene2=501-1000;
  partition by-codpos=4: gene1_cp1, gene1_cp2, gene1_cp3, gene2;   
  ```

  More details [in the `MrBayes` documentation](https://github.com/NBISweden/MrBayes/blob/develop/doc/manual/Manual_MrBayes_v3.2.pdf) and if you type `help charset` and `help partition` when running `MrBayes` interactively.

* `set`: this command is to be used alongside command `partition` because it "sets up" what you previously defined.
* `lset`: this command sets the parameters of the likelihood model. There are different options this command can take, but we will focus on `nst`, `applyto`, and `rates`. More details [in the `MrBayes` documentation](https://github.com/NBISweden/MrBayes/blob/develop/doc/manual/Manual_MrBayes_v3.2.pdf) and if you type `help lset` when running `MrBayes` interactively.
* `unlink`: this command is to be used when you enable the `partition` command. As the name of this command says, it will "unlink" model parameters across the data partitions you have defined. You can type `all` or the specific name of the partition/s for which you want to unlink the model parameters. By default, if the same parameter applies to different partitions and if this parameter has the same prior, `MrBayes` will use a single value for this parameter. If you want to use different parameter values for each partition you have established, then you need to use this command to "unlink" the model parameters, and specific parameter values will be inferred for each partition (see command below). If you use the command `link` instead of `unlink`, the opposite will occur. More details [in the `MrBayes` documentation](https://github.com/NBISweden/MrBayes/blob/develop/doc/manual/Manual_MrBayes_v3.2.pdf) and if you type `help unlink` when running `MrBayes` interactively.
* `prset`: use this command to set the priors for the phylogenetic model you want to use. This command enables various options, but we will focus on `statefreq` (stationary nucleotide frquencies), `shape` (shape parameter of the Gamma distribution of rate variation), and `tratio` (kappa parameter for HKY, the one we will be running our analyses under) when analysing the partitioned dataset as we are "unlinking" these model parameters for partitions 1 and 2. More details [in the `MrBayes` documentation](https://github.com/NBISweden/MrBayes/blob/develop/doc/manual/Manual_MrBayes_v3.2.pdf) and if you type `help prset` when running `MrBayes` interactively.
* `mcmc`: this command is used to set up and start the MCMC analysis. There are different options this command can activate, but we will focus on `seed`, `ngen`, `nruns`, `nchains` (default is 4, 1 cold chain and 3 heated chains), `printfreq`, `samplefreq`, `diagnfreq`, `diagnstat`, `savebrlens`, and `filename`. More details [in the `MrBayes` documentation](https://github.com/NBISweden/MrBayes/blob/develop/doc/manual/Manual_MrBayes_v3.2.pdf) and if you type `help mcmc` when running `MrBayes` interactively.

> **BLOCK 3: Summarise trees and other model parameters (as many blocks as analysis you want to perform!)**:

* `sumt`: this command produces summary statistics for the trees that have been sampled during the MCMC. You can specify the file name (`filename`) where you want the output to be written. By default, the burnin is established to be 25% of the samples collected (you could modify this, if required). More details [in the `MrBayes` documentation](https://github.com/NBISweden/MrBayes/blob/develop/doc/manual/Manual_MrBayes_v3.2.pdf) and if you type `help sumt` when running `MrBayes` interactively.
* `sump`: this command prints the values that have been sampled for the model parameters during the MCMC. You can specify the file name (`filename`) where you want the output to be written. By default, the burnin is established to be 25% of the samples collected (you could modify this if required). More details [in the `MrBayes` documentation](https://github.com/NBISweden/MrBayes/blob/develop/doc/manual/Manual_MrBayes_v3.2.pdf) and if you type `help sump` when running `MrBayes` interactively.

> **BLOCK 4: Stop `MrBayes` and end of nexus file**:

* `quit`: this command exits `MrBayes`.
* `end`: this command is used to indicate that the `MrBayes` block has come to an end.

### 2. Run `MrBayes`

Now, we are ready to run `MrBayes`!

```sh
# Run from my_session/day2/mrbayes
# Change directories if you are not
# there yet

# Go to `inp_dat` dir and run `MrBayes`
cd inp_data
mb ctl_mb_dna.nex

# Move results to another dir
mkdir ../{mb_dnaconc_out,mb_dnapart_out}
mv dna-conc* ../mb_dnaconc_out
mv dna-part* ../mb_dnapart_out
mv *log.txt ../
```

> [!IMPORTANT]
> If these analyses take too long, we will start analysing the results that have been generated prior to the practical session to save time. If that's the case, please run the following commands:
>
> ```sh
> # Run from my_session/day2/mrbayes
> cp -R out_of_time/mb_dnaconc_out .
> cp -R out_of_time/mb_dnapart_out .
> cp out_of_time/*log.txt .
> ```

We will go through the output files together and, subsequently, we will learn how to use `Tracer` and `FigTree` to analyse and visualise the results obtained!

### 3. Analysing the `MrBayes` MCMC output

To analyse the MCMC output, we are going to use the program `Tracer`, which you should have installed on your PC.

Please run `Tracer` by double clicking on the `Tracer` icon.

> [!TIP]
> If you want to run `Tracer` from the command line, you can either (i) set an alias in your `~/.bashrc` or `~/.bash_profile` file, (ii) execute the program from the directory where the `jar` file can be found, or (iii) execute this `jar` file using an absolute or a relative path to the directory where the terminal is running.

<details>
<summary><b>[ Click here to learn more about how you could run <code>Tracer</code> on the terminal ]</b></summary>
<br>

```sh
# Option 1: Run from the directory where 
# `tracer.jar` can be found
# For instance, if you open a terminal from
# `Tracer_v1.7.2/lib` (or the corresponding
# path on your PC), then type the following
java -jar tracer.jar

# Option 2: Use an absolute or a relative 
# path to execute the file
# An example of an absolute path is shown
# below -- please note this may differ from
# yours, so change the command below accordingly
java -jar Applications/Tracer_v1.7.2/lib/tracer.jar

# Option 3: Add an alias in your 
# `~/.bashrc` or `~/.bash_profile
# Please replace "<path_to_Tracer>" with the absolute path
# to the location where you have saved `Tracer`
# Modify the name of the directory that you have unzipped
# if needed too (e.g., `Tracer_v1.7.2`)
alias tracer1.7.2='java -jar <path_to_Tracer>/Tracer_v1.7.2/lib/tracer.jar'
# Now, open a terminal from any location on 
# your PC and type the following command to 
# execute `Tracer`
# Before you do this, please  make sure that you
# have X11 installed or, if you are on Windows,
# that you have the Xming 
# Server running on your PC
# If you are on Windows, you might need to run 
# the command `export DISPLAY=:0.0` before 
# you can execute `Tracer` from your WSL
tracer1.7.2
```

</details><br>

Regardless of the approach you have used to open `Trace`, please load the output files that have saved the samples collected during the MCMC for all model parameters (i.e., files which extension is `.p`). You can select the `Import Trace File...` option from the `File` menu. Then, select the file with extension `.p` that was output by `MrBayes` to load it onto `Tracer`. Alternatively, you can also drag the file onto the `Tracer` icon.

We will go through the most important features of `Tracer` together but, in general, we will focus on the effective sample sizes (ESSs) calculated for each of the model parameters, the frequency plot of the samples, and the trace plots.

> [!IMPORTANT]
> Please take sometime to answer the following questions:
>
> * Do you think we need to run the chains longer?
> * Is the ESS large enough for all model parameters?
> * How efficient is the chain?

### 4. Viewing the annotated tree

As you have already learnt when visualising the best-scoring ML trees, `FigTree` is a user-friendly, graphical program for viewing trees. You can run it either by double-clicking on the icon or via the command line (you will need to follow the same procedure as described above for Tracer).

<details>
<summary><b>[ Click here to learn more about how you could run <code>FigTree</code> on the terminal ]</b></summary>
<br>

```sh
# Option 1: Run from the directory where 
# `figtree.jar` can be found
# For instance,  if you open a terminal from
# `FigTree_v1.4.4/lib` (or corresponding path
# on your PC), please type the following
java -jar figtree.jar

# Option 2: Use an absolute or a relative 
# path to execute the file
# An example of an absolute path is shown below,
# but please update the command below to match
# your settings
java -jar Applications/FigTree_v1.4.4/lib/figtree.jar

# Option 3: Add an alias in your 
# `~/.bashrc` or `~/.bash_profile
# Please replace "<path_to_FigTree>" with the absolute path
# to the location where you have saved `FigTree`
# Modify the name of the directory that you have unzipped
# if needed too (e.g., `FigTree_v1.4.4`)
alias figtree1.4.4='java -jar <path_to_Tracer>/FigTree_v1.4.4/lib/figtree.jar'
# Now, open a terminal from any location on 
# your PC and type the following command to 
# execute `FigTree`
# Before you do this, please make sure that yo
# have X11 installed or, if you are on Windows,
# that you have the Xming Server running on your PC
# If you are on Windows, you might need to run 
# the command `export DISPLAY=:0.0` before 
# you can execute `FigTree` from your WSL
figtree1.4.4
```

</details><br>

Now, you can launch `FigTree` to view the file with the consensus tree that `MrBayes` has output (i.e., file names that end with `*con.tre`). The tree will be displayed in the `FigTree` window. On the left hand side, you can find the options and settings which control how the tree is displayed. We will see together the main options you can use to display the tree.


## Tree and Alignment Comparison with Phykit

There are various methods one can use to evaluate trees and alignments.  One package that contains a number of useful comparison metrics is [Phykit](https://jlsteenwyk.com/PhyKIT/).  

We will demonstrate how to use Phykit to compare Maximum Likelihood trees based on trimmed and untrimmed AA alignments between structural and sequence alignments.

First we need to create the trees.  

Earlier you created a tree for the structure-based AA alignment.  The alignment is saved in `~/my_session/day2/iq-tree/aa/aln_aa.fasta`.  We will now make another alignment based on the protein tree created by Muscle5 on day1 from `my_session/day1/aln/aln_prot_super5_muscle.fasta aa_seq/aln_aa_seq.fasta` 

```sh
# Run from my_session/day2/iq-tree

#Make a directory aa_seq and then copy the muscle alignment file into that folder
mkdir aa_seq
cp ../../day1/aln/aln_prot_super5_muscle.fasta aa_seq/aln_aa_seq.fasta

#run iqtree
cd aa_seq
iqtree2 -s aln_aa_seq.fasta -B 1000 -T AUTO

```

For this alignment notice that a different model of evolution was determined to be most appropriate to construct the tree than the structure-based alignment: `Q.plant+G4` rather than `Q.plant+I+R3`, which does not include an invariant class, and uses a discrete gamma distribution with 4 categories for rate heterogeneity rather than a free rate model with 3 categories.

Phykit has a number of functions to compare alignments, trees, or both.  Here are a few that might be useful to you and more can be found in the [phykit documentation](https://jlsteenwyk.com/PhyKIT/usage/index.html)

<table>
<!-- HEADER -->
<tr>
<th>Comparison Name</th>
<th>Data Used</th>
<th>Phykit function </th>
<th>Description</th>
</tr>

<!-- 1st ROW -->
<tr>
<td>Pairwise Identity</td>
<td>alignment</td>
<td>pairwise_identity</td>
<td>Average pairwise identity among sequences.  Defined as the number of identical columns (including gaps) between two aligned sequences divided by the number of columns in the alignment.  </td>
</tr>
<!-- 2nd ROW -->
<tr>
<td>Alignment Length</td>
<td>alignment</td>
<td>aln_len</td>
<td>Length of alignment</td>
</tr>
<!-- 3rd ROW -->
<tr>
<td>Relative composition variability</td>
<td>alignment</td>
<td>relative_composition_variability</td>
<td>Describes the average variability in sequence composition among taxa.</td>
</tr>
<!-- 4th ROW -->
<tr>
<td>Column score</td>
<td>alignment</td>
<td>column_score</td>
<td>Compare an alignment to a reference alignment.  Ratio of correctly aligned columns to total number of columns in an alignment.</td>
</tr>
<!-- 5th ROW -->
<tr>
<td>Internal Branch Statistics</td>
<td>tree</td>
<td>internal_branch_stats</td>
<td>Summary statistics for internal branch lengths.</td>
</tr>
<!-- 6th ROW -->
<tr>
<td>Terminal Branch Statistics</td>
<td>tree</td>
<td>terminal_branch_stats</td>
<td>Summary statistics for terminal branch lengths.</td>
</tr>
<!-- 8th ROW -->
<tr>
<td>Total Tree length</td>
<td>tree</td>
<td>total_tree_length</td>
<td>Sum of all branches.</td>
</tr>
<!-- 7th ROW -->
<tr>
<td>Treeness</td>
<td>tree</td>
<td>treeness</td>
<td>Sum of internal branch lengths divided by the total tree length.</td>
</tr>
<!-- 9th ROW -->
<tr>
<td>Robinson foulds distance</td>
<td>tree</td>
<td>robinson_foulds_distance</td>
<td>After pruning to include only shared tips, provides a metric to compare two phylogenies.  A low RF distance indicates the two trees are similar.  The function provides a plain RF value and a normalized RF value (divided by 2(N-3) where N is the number of tips). </td>
</tr>
<!-- 10th ROW -->
<tr>
<td>Treeness/RCV</td>
<td>tree and aligment</td>
<td>treeness_over_rcv</td>
<td>Divides the Treeness (a tree-based value) by the RCV (an alignment-based value).  High values indicate a high signal to noise ratio (higher treeness) and low composition bias (lower RCV)</td>
</tr>
<!-- 11th ROW -->
<tr>
<td>Saturation</td>
<td>tree and alignment</td>
<td>saturation</td>
<td>Assesses how much the sequences in the alignment have undergone numerous substitutions such that the distances between taxa are underestimated.  Compares the patristic distance (that is the distance between two sequences based on the branches of the phylogenetic tree which take multiple substitutions into account) and the raw uncorrected distance computed between the two sequences in the alignment.  See [Philippe et al. 2011](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1000602).</td>
</tr>
</table>

We will start by comparing both alignments based on column score and both trees based on Robinson-Foulds distance. 

For column score to work correctly, the sequences in the alignment should be in the same order.  We have provided a small python script that puts the sequences from aln_aa_seq.fasta into the same order as aln_aa.fasta. 

```sh
# Run from my_session/day2/iq-tree

cp ../../../biol0033-tutorial/day2/scripts ../
chmod 775 ../scripts/* 


# Phykit is installed in the virtual environment py-env along with biopython which is used in the python script.  You need to activate the environment.   
source /opt/py-env/bin/activate

../scripts/reorder_aln.py aa_seq/aln_aa_seq.fasta aa/aln_aa.fasta aa_seq

#Get Column Score with the sequence-based alignment as the reference
phykit column_score aa/aln_aa.fasta --reference aa_seq/aln_aa_seq_reordered.fasta

#Get Column Score with the structure-based alignment as the reference
phykit column_score aa_seq/aln_aa_seq_reordered.fasta --reference aa/aln_aa.fasta
```

You should get something like: 

```sh
(py-env) [<username>@rstudio-biol0033 iq-tree]$ phykit column_score aa/aln_aa.fasta --reference aa_seq/aln_aa_seq_reordered.fasta
0.3249
(py-env) [<username>@rstudio-biol0033 iq-tree]$ phykit column_score aa_seq/aln_aa_seq_reordered.fasta --reference aa/aln_aa.fasta
0.403
```

An identical alignment would have a column score of 1, but when we do this, we see a number less than 1. 


```sh
# Run from my_session/day2/iq-tree

#Get Column Score with the sequence-based alignment as the reference
phykit column_score aa/aln_aa.fasta --reference aa/aln_aa.fasta

#Get Column Score with the structure-based alignment as the reference
phykit column_score aa_seq/aln_aa_seq_reordered.fasta --reference aa_seq/aln_aa_seq_reordered.fasta

```

```sh
(py-env) [<username>@rstudio-biol0033 iq-tree]$ phykit column_score aa/aln_aa.fasta --reference aa/aln_aa.fasta
0.6973
(py-env) [<username>@rstudio-biol0033 iq-tree]$ phykit column_score aa_seq/aln_aa_seq_reordered.fasta --reference aa_seq/aln_aa_seq_reordered.fasta
0.6796
```
This is because columns with gaps incur a penalty in the column score.  For more information see [Thompson et al. 1999](https://academic.oup.com/nar/article/27/13/2682/2376831).


Now calculate the Robinson-Foulds Distance based on the trees.  You do not need to reorder the leaves or designate a reference to calculate the Robinson Foulds Distance. The output is two values, the raw RF distance and the normalized RF distance (See [Robinson & Foulds 1981](https://www.sciencedirect.com/science/article/pii/0025556481900432) for details)

```sh
# Run from my_session/day2/iq-tree
# Get Robinson Foulds distance between the structure and sequence generated trees
phykit robinson_foulds_distance aa/aln_aa.fasta.treefile aa_seq/aln_aa_seq.fasta.treefile

```

We get 0 for both values, as both of these trees are topologically equivlent (that is not to say their branch lengths are the same). 

```sh
(py-env) [<username>@rstudio-biol0033 iq-tree]$ phykit robinson_foulds_distance aa/aln_aa.fasta.treefile aa_seq/aln_aa_seq.fasta.treefile
0       0.0
```

Now we will gather some more metrics about each tree and alignment.  For this we will make a directory to save the data and use a small script, `my_session/day2/scripts/phykit_comparisons.py` that takes a file listing input alignments and trees that we want to analyze and another file listing the metrics we want to calculate and prints the output to a text file.  The inputs are provided in the day2 directory.   


```sh
# Run from my_session/day2
# If not already there change to that directory

cp /biol0033/day2/phykit_comparisons -R .

scripts/phykit_comparisons.py phykit_comparisons/seq_struct_to_compare.txt phykit_comparisons/phykit_metrics.csv phykit_comparisons/seq_struct.txt

```

This produces the file  `my_session/day2/phykit_comparisons/seq_struct.txt`

Which has a number of metrics including: 

```sh
Metric: Alignment Length
========================================

Input: aa_struct
----------------------------------------
991

Input: aa_seq
----------------------------------------
799
```

```sh
Metric: Terminal Branch Statistics
========================================

Input: aa_struct
----------------------------------------
mean: 0.6264
median: 0.2954
25th percentile: 0.0803
75th percentile: 0.4813
minimum: 0.0263
maximum: 4.1386
standard deviation: 1.1462
variance: 1.3137

Input: aa_seq
----------------------------------------
mean: 0.1634
median: 0.1235
25th percentile: 0.0339
75th percentile: 0.1938
minimum: 0.0098
maximum: 0.5198
standard deviation: 0.1748
variance: 0.0305
```

```sh
Metric: Tree length
========================================

Input: aa_struct
----------------------------------------
11.9309

Input: aa_seq
----------------------------------------
3.8441
```

```sh
Metric: Treeness/RCV
========================================

Input: aa_struct
----------------------------------------
4.7599	0.37	0.0777

Input: aa_seq
----------------------------------------
5.0645	0.49	0.0968
```

We see in particular that the alignment length for the structure-based alignment is much longer than the sequence based alignment.  It is likely that this is due to the large gaps in the beginning of the alignment which occur because of non-overlapping disordered N-terminal regions. 

### Trimming Alignments

We can try to tidy up the alignment by clipping out very gappy regions.  

To do this we will use the software, [ClipKit](https://jlsteenwyk.com/ClipKIT/) 

We will use the default mode which implements the smart-gap dynamic algorithm, but there are various other modes available. 

The flag `-l` ensures that a log is printed out that tracks which alignment columns were trimmed. 

The output is a multi-line fasta so we change it to a single line with our `one_line_fasta.pl` script

```sh
# Run from my_session/day2
# If not already there change to that directory
mkdir {aa,aa_seq}_trimmed

#trim structural alignment
clipkit aa/aln_aa.fasta -o aa_trimmed/aln_aa_trimmed.fasta -l

#Convert the fasta file to have one line per sequence
../../day1/scripts/one_line_fasta.pl aa_trimmed/aln_aa_trimmed.fasta
mv aa_trimmed/aln_aa_trimmed_one_line.fa aa_trimmed/aln_aa_trimmed.fasta

#trim sequence alignment
clipkit aa_seq/aln_aa_seq.fasta -o aa_seq_trimmed/aln_aa_seq_trimmed.fasta -l

#Convert the fasta file to have one line per sequence
../../day1/scripts/one_line_fasta.pl aa_seq_trimmed/aln_aa_seq_trimmed.fasta
mv aa_seq_trimmed/aln_aa_seq_trimmed_one_line.fa aa_seq_trimmed/aln_aa_seq_trimmed.fasta

```

You should get an output that provides details on the clipping, such as this one for the structural alignment

```sh
-------------
| Arguments |
-------------
Input file: aa/aln_aa.fasta (format: fasta)
Output file: aa_trimmed/aln_aa_trimmed.fasta (format: fasta)
Sequence type: Protein
Gaps threshold: 0.9167
Gap characters: ['-', '?', '*', 'X', 'x']
Trimming mode: smart-gap
Create complementary output: False
Process as codons: False
Trim ends only: False
Create log file: True


------------------------
| Writing output files |
------------------------
Trimmed alignment: aa_trimmed/aln_aa_trimmed.fasta
Complement file: False
Log file: aa_trimmed/aln_aa_trimmed.fasta.log


---------------------
| Output Statistics |
---------------------
Original length: 991
Number of sites kept: 740
Number of sites trimmed: 251
Percentage of alignment trimmed: 25.328%

Execution time: 0.046s
```

The trimming for the sequence alignment removed fewer sites, and ended up with fewer columns.  
```sh
---------------------
| Output Statistics |
---------------------
Original length: 799
Number of sites kept: 699
Number of sites trimmed: 100
Percentage of alignment trimmed: 12.516%
```

Now we build new trees with the Trimmed alignments: 


```sh
# Run from my_session/day2/iq-tree

#run iqtree
cd aa_trimmed
iqtree2 -s aln_aa_trimmed.fasta -B 1000 -T AUTO

cd ../aa_seq_trimmed
iqtree2 -s aln_aa_seq_trimmed.fasta -B 1000 -T AUTO

```

IQ-tree generated the same evolutionary model for trees generated with the trimmed alignments as it did for the untrimmed alignments:  `Q.plant+I+R3` for the structure-based tree and `Q.plant+G4` for the sequence-based tree. 

Now we can compare the various metrics between trimmed and untrimmed alignments


```sh
# Run from my_session/day2
# If not already there change to that directory

scripts/phykit_comparisons.py phykit_comparisons/seq_struct_trim_to_compare.txt phykit_comparisons/phykit_metrics.csv phykit_comparisons/seq_struct_trim.txt

```

We see that, besides alignment length, many of these metrics are similar between trimmed and untrimmed alignments.  Indeed, for the tree-based metrics the values are almost identical: 

```sh

Metric: Alignment Length
========================================

Input: aa_struct
----------------------------------------
991

Input: aa_seq
----------------------------------------
799

Input: aa_struct_trimmed
----------------------------------------
740

Input: aa_seq_trimmed
----------------------------------------
699


Metric: Internal Branch Statistics
========================================

Input: aa_struct
----------------------------------------
mean: 0.4905
median: 0.2775
25th percentile: 0.1503
75th percentile: 0.6303
minimum: 0.0538
maximum: 1.9568
standard deviation: 0.5952
variance: 0.3543

Input: aa_seq
----------------------------------------
mean: 0.2093
median: 0.1176
25th percentile: 0.0767
75th percentile: 0.2451
minimum: 0.0199
maximum: 0.8393
standard deviation: 0.2542
variance: 0.0646

Input: aa_struct_trimmed
----------------------------------------
mean: 0.4904
median: 0.2775
25th percentile: 0.1503
75th percentile: 0.6302
minimum: 0.0538
maximum: 1.9564
standard deviation: 0.5951
variance: 0.3541

Input: aa_seq_trimmed
----------------------------------------
mean: 0.2093
median: 0.1176
25th percentile: 0.0767
75th percentile: 0.2451
minimum: 0.0199
maximum: 0.8393
standard deviation: 0.2542
variance: 0.0646

Metric: Tree length
========================================

Input: aa_struct
----------------------------------------
11.9309

Input: aa_seq
----------------------------------------
3.8441

Input: aa_struct_trimmed
----------------------------------------
11.9285

Input: aa_seq_trimmed
----------------------------------------
3.844
```

What about the column score?  After reordering the trimmed sequence based tree we can recalculate the column score with the trimmed alignments.  


```sh
# Run from my_session/day2

# reorder sequence-based trimmed alignment to match structure based trimmed alignment so we can compare them with column score. 
scripts/reorder_aln.py iq-tree/aa_seq_trimmed/aln_aa_seq_trimmed.fasta iq-tree/aa_trimmed/aln_aa_trimmed.fasta iq-tree/aa_seq_trimmed/

#Get Column Score with the sequence-based alignment as the reference
phykit column_score iq-tree/aa_trimmed/aln_aa_trimmed.fasta --reference iq-tree/aa_seq_trimmed/aln_aa_seq_trimmed_reordered.fasta

#Get Column Score with the structure-based alignment as the reference
phykit column_score iq-tree/aa_seq_trimmed/aln_aa_seq_trimmed_reordered.fasta --reference iq-tree/aa_trimmed/aln_aa_trimmed.fasta

#Get the baseline Column Score for the structure-based alignment vs itself 
phykit column_score iq-tree/aa_trimmed/aln_aa_trimmed.fasta --reference iq-tree/aa_trimmed/aln_aa_trimmed.fasta

#Get the baseline Column Score with the sequence-based alignment vs itself 
phykit column_score iq-tree/aa_seq_trimmed/aln_aa_seq_trimmed_reordered.fasta --reference iq-tree/aa_seq_trimmed/aln_aa_seq_trimmed_reordered.fasta

```

We see that for the trimmed structure-based alignment, there seems to be an increased similarity with the sequence-based alignment.  Also the control baseline column score is much higher for the trimmed structure-based tree. 

<table>
<!-- HEADER -->
<tr>
<th>Comparison for Column Score</th>
<th>Untrimmed</th>
<th>Trimmed </th>
</tr>

<!-- 1st ROW -->
<tr>
<td>Structure-based vs Sequence-based (ref)</td>
<td>0.3249</td>
<td>0.377</td>
</tr>
<!-- 2nd ROW -->
<tr>
<td>Sequence-based vs Structure-based (ref)</td>
<td>0.403</td>
<td>0.3991</td>
</tr>
<!-- 3rd ROW -->
<tr>
<td>Structure-based vs Structure-based</td>
<td>0.6973</td>
<td>0.773</td>
</tr>
<!-- 4th ROW -->
<tr>
<td>Sequence-based vs Sequence-based</td>
<td>0.6796</td>
<td>0.7082</td>
</tr>
</table>


Which alignment / tree is better will depend on what you prioritize in your analysis as each metric measures different properties.  It appears that in this case, the structure-based tree gives longer branches and has higher saturation.  Treeness measures the internal branch lengths divided by the total tree length and the sequence-based trees had a higher value than the structure-based trees (0.49 vs 0.37) which would tend to indicate a higher signal to noise ratio.  Relative Compositional Variability measures bias in the composition of amino acids between sequences and high RCV indicates that some branches have specific amino acids overrepresented which could lead to systematic errors.  It is lower in structure-based alignments than in sequence-based alignments (0.0777/0.0782 vs 0.0968/0.0906).  Treeness/RCV is the ratio between these two metrics and was better for the sequence-based alignments. 



---

This is the end of day 2! Hope you have enjoyed reconstructing phylogenies under different approaches, models of evolution, and datasets! :smiley:
