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

The preferred model for the DNA alignment is `Q.plant+I+R3`. Let's evaluate each part:

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

Once we have defined our input alignment in NEXUS format, we can generate the control file that will execute `MrBayes`, which will also be in NEXUS format. You can find a template control file in [directory `inp_data`](mrbayes/inp_data/ctl_mb_dna.nex). Below, you will find a summary of the commands that we will go through:

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

---

This is the end of day 2! Hope you have enjoyed reconstructing phylogenies under different approaches, models of evolution, and datasets! :smiley:
