# Day 3 - testing for positive selection

> [!IMPORTANT]
> **Disclaimer**: this tutorial is based on the article **Beginner's Guide on the Use of PAML to Detect Positive Selection** ([Álvarez-Carretero et al. 2023](https://doi.org/10.1093/molbev/msad041)). Consequently, you may find similarities between the aforementioned paper and its [step-by-step protocol](https://github.com/abacus-gene/paml-tutorial/tree/main/positive-selection) when going through this practical session.

## Introduction

In the theoretical lectures, you have learnt that the ratio between nonsynonymous and synonymous substitution rates, $\omega=dN/dS$, has been widely used as a model parameter to measure the effect of natural selection on protein-coding genes ([Kimura 1977](https://www.nature.com/articles/267275a0), [Miyaga & Yasunaga 1980](https://link.springer.com/article/10.1007/BF01732067)). While synonymous mutations (also known as "silent" mutations) do not change the amino acid, nonsynonymous mutations (also known as replacement mutations) do. The latter mutations may often be deleterious and purged from the population by natural selection, resulting in reduced nonsynonymous substitution rate ($d_{N}$). **The $\omega$ ratio measures the direction and magnitude of such selection on nonsynonymous mutations**.

Once we calculate $\omega$, we may be able to suggest one of the following:

* The analysed protein-coding gene is **under positive (or diversifying) selection (or adaptive evolution)**, if $\omega>1$. Beneficial nonsynonymous mutations are favoured and quickly fix in the population.
* The analysed protein-coding gene shows **no evidence for selection**, if $\omega=1$. In such a case, most mutations are neither beneficial nor detrimental (i.e., the mutation makes no fitness difference), thus often related to neutral mutations.
* The analysed protein-coding gene may be **under negative (or purifying) selection**, if $\omega<1$.

> [!IMPORTANT]
> After running `CODEML`, we will run **Likelihood Ratio Tests** (LRTs) to assess the fitness of nested models. In some model comparisons, we will see that some models in which $\omega=1$ is fixed to account for sites that are nearly neutral or under a weak constraint become our null models.

As every functional protein has some structural constraints, the $\omega$ ratio average over the whole protein sequence is often less than 1, even if positive selection operates at some sites or over certain time intervals. Thus, statistical tests have been designed to detect positive selection that targets only a small subset of the amino acid residues or affects a limited time interval, leading to the development of **site** models and **branch** models, which subsequently led to the the **branch-site** models.

During this practical session, we will learn how to test for positive selection under the codon models implemented in the the maximum-likelihood program **`CODEML`**, part of the `PAML` software ([Yang 2007](https://doi.org/10.1093/molbev/msm088)). Particularly, we shall focus on the following analyses:

* Calculating $\omega$ (i.e., nonsynonymous/synonymous ratio $d_{N}/d_{S}$) as a measure of average selective pressure on the protein-coding gene we are analysing (i.e., Mx gene) under the **homogeneous model**, which assumes one $\omega$ for all sites and branches.
* Detecting positive selection affecting a **subset of sites** in the protein-coding gene evolving under positive selection (**site test**).
* Detecting a **specific branch or branches** of a phylogeny evolving under positive selection (**branch test**).
* Detecting a **subset of sites for particular branches** of a phylogeny evolving under positive selection (**branch-site test**).

> [!IMPORTANT]
> After running `CODEML`, we will run **Likelihood Ratio Tests** (LRTs) to assess the fitness of nested models. The homogeneous model will be our null hypothesis in some of our tests for positive selection.

In order to start our analyses with `CODEML`, we will need two important input files: a **fixed tree topology** and a **sequence alignment**. During days 1 and 2, we have learnt how to infer both so... We have everything we need to get started!

## Using `PAML` program `CODEML` to detect positive selection

### Understanding the control file

To run any `PAML` program, we need a control file. In other words, instead of having all options passed through a command, all variables and the corresponding options are specified in a text file that the program reads. Below, you can find an example of control file to run `CODEML` that you can always use as a template (options have been kept as "flags", which should be replaced with adequate values):

```txt
* [ INPUT/OUTPUT FILES ]
      seqfile = ALN            * Path to the sequence alignment file
     treefile = TREE           * Path to the tree file
      outfile = OUT            * Path to the output file
   
* [ SCREEN OUTPUT ]
        noisy = 3              * How much information you want to see on the screen
      verbose = 1              * Detailed output file

* [ TYPE OF SEQUENCE DATA ]
      seqtype = 1              * Data type
        ndata = NDAT           * Number of alignment blocks / partitions
        icode = 0              * Genetic code 
    cleandata = 0              * Remove sites with ambiguity data?

* [ SUBSTITUTION MODEL ]
        model = CODMOD         * Models for ω varying across lineages
      NSsites = NSSIT          * Models for ω varying across sites
    CodonFreq = CODFREQ        * Codon frequencies
      estFreq = ESTFREQ        * Use observed freqs or estimate freqs by ML
    fix_omega = FIXOME         * Estimate or fix omega
        omega = INITOME        * Initial or fixed omega
        clock = CLOCK          * Clock model
```

Let's explain each variable and their options!

* Variables related to input/output files:
  * `seqfile`: flag `ALN` needs to be replaced with the name of the input file with the sequence alignment (if available within the same directory where this control file is saved) or the path to such file. When the latter, you can use either absolute or relative paths to the input sequence file.
  * `treefile`: flag `TREE` needs to be replaced with the name of the input tree file (if available within the same directory where this control file is saved) or the path to such file. When the latter, you can use either absolute or relative paths to the input tree file.
  * `outfile`: flag `OUT` needs to be replaced with the name of the output file (if you want to save it in the same directory where this control file is saved) with relevant information about the analysis. If you want to save this file elsewhere, you can use either absolute or relative paths to do so.
* Variables related to screen and output details:
  * `noisy` and `verbose`: these are two variables that control how much information you want to see printed on the screen and on your output file, respectively. To see a moderate amount of detail printed on the screen, you may choose `noisy = 3`(possible options: `0` to `3`). If you want to obtain a very detailed output, you should choose `verbose = 1`. Choose `verbose = 0`  for a much more concise output file.
* Variables related to the input data type:
  * `seqtype`: this variable indicates the data type. When running `CODEML`, you should always use `seqtype = 1` as you must use codon data. The other options are available for other `PAML` programs that enable the usage of nucleotide data (`seqtype = 0`) or amino acid data (`seqtype = 2`).
  * `ndata`: flag `NDAT` should be replaced with the number of alignment blocks (or partitions) that you have in your input sequence alignment. We are only analysing one gene in this tutorial, so `NDAT` should be replaced with `1`.
  * `icode`: this variable related to the genetic code. When setting `icode = 0`, we are specifying the standard genetic code. All genetic codes that can be specified in `CODEML` are enabled using the following options (see [`transl_table` 1 to 11 of GENEBANK](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1) for more details on each of the below):
    * `0`: standard.
    * `1`: mammalian mt.
    * `2`: yeast mt.
    * `3`: mold mt.
    * `4`: invertebrate mt.
    * `5`: ciliate nuclear
    * `6`: echinoderm mt.
    * `7`: euplotid mt.
    * `8`: alternative yeast nu.
    * `9`: ascidian mt.
    * `10`: blepharisma nu.
  * `cleandata`: this variable can be used to decide whether sites with ambiguity data and alignment gaps should be kept (`cleandata = 0`) or removed (`cleandata = 1`). It is recommended that all alignment filtering and trimming takes place prior to running `CODEML`, thus encouraging users to specify `cleandata = 0`. More recommendations about this step can be found in sections `Alignment, sequence data file, and tree file" and "Gene tree versus species tree" in the [Supplementary Material of Álvarez-Carretero et al. 2023](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/mbe/40/4/10.1093_molbev_msad041/2/msad041_supplementary_data.pdf?Expires=1774651114&Signature=phKHDR~EitWyE~wSuNTBeR1jhwChp3q7E5SFRBv~u9Z1uVP~Q7gUJNs-4vRJxxYdgg7uBidCOFxKOotoLJeJMmFdqTs3AC3VLROgn~mSVgrICce7C6OOE-LqS4Dh3uM9NLhSKEOqtxKMfaTpcgTNXTHW2woP9WHwIP2i-QVV~W9z96ZuDFBg6HhqT3xFyyYmhI2JNDe~L8TGuSLDjI-k2kACz3GBt46-GYSS9dPJ75GexYd-hJ1G-s4Ravm329wSba4UmYdyeNHCV0QB4GKrAx1tIRoE6NOdWMS-4YDJuMdAyq2IdIvOK-xJF2MFFgx5vB5Ip-HJqAWAAlPIFO4RXw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA).
* Variables related to specifying the codon model:
  * `model` and `NSsites`: these two variables can jointly specify which model is being enabled:
    * homogeneous model: `model = 0` and `NSsites = 0`.
    * Site models: `model = 0` and `NSsites` can equal `1`, `2`, `7`, or `8` (depending on the site model under which the analysis takes place).
    * Branch model: `model = 2` and `NSsites = 0`.
    * Branch-site model: `model = 2` and `NSsites = 2`.
    <p align="center">
    <img width="350" height="400" src="../figs/CODEML_models.jpeg">
    </p>

    > Figure 1. Illustration of four different types of models implemented in `CODEML`, from Figure 3 in [Álvarez-Carretero et al. 2023](https://doi.org/10.1093/molbev/msad041). (A) homogeneous evolutionary pressure throughout the history of the gene (`M0`: one ratio, with one $\omega$ ratio for all sites and branches, specified as `model = 0` and `NSsites = 0`); (B) heterogeneous pressure across codons (site models: `model = 0` and `NSsites = 1, 2, 7, 8`, etc.); (C) heterogeneous pressure across branches of a tree but homogeneous across codons (branch model: `model = 2` and `NSsites = 0`); and (D) heterogeneous pressure across sites and branches (branch-site model: `model = 2` and `NSsites = 2`).
  * `CodonFreq`: the codon frequencies can be set to either equal or different to account for codon usage bias. If you choose `verbose = 1`, the program will print out the codon frequencies expected under the latter models, calculated using the observed nucleotide frequencies in the dataset. The options accepted for variable `CodonFreq` are as follows:
    * `0`: codon frequencies are assumed to be equal (i.e., 1/61 each). Number of free parameters for standard genetic code: 0.
    * `1`: "F1X4". Codon frequencies are calculated from the average of nucleotide frequencies. Number of free parameters for standard genetic code: 3.
    * `2`: "F3X4". Codon frequencies are calculated from the average nucleotide frequencies at the three codon positions. Number of free parameters for standard genetic code: 9.
    * `3`: codon frequencies are used as free parameters. Number of free parameters for standard genetic code: 60.
    * `4`: "F1x4MG". This option uses one set of base frequencies for all three codon positions just like F1x4 but is implemented in the style of [Muse and Gaut (1994)](https://doi.org/10.1093/oxfordjournals.molbev.a040152): the substitution rate from codon _i_ to codon _j_ is proportional to the frequency of the target nucleotide rather than the frequency of the target codon. Number of free parameters for standard genetic code: 3.
    * `5`: "F3X4MG". This option calculates codon frequencies from the average nucleotide frequencies at the three codon positions just like F3x4 but is implemented in the style of [Muse and Gaut (1994)](https://doi.org/10.1093/oxfordjournals.molbev.a040152): the substitution rate from codon _i_ to codon _j_ is proportional to the frequency of the target nucleotide rather than the frequency of the target codon. Number of free parameters for standard genetic code: 9.
    * `6`: "FMutSel0". This is a special case of "FMutSel" (see below) that assigns the same fitness value for all synonymous codons, so that only 19 (= 20 – 1) amino acid fitness parameters are used. This model assumes that the amino acid frequencies are determined by the functional requirements of the protein and, given the amino acid frequencies, the relative frequencies of synonymous codons are determined solely be the mutational-bias parameters. Number of free parameters for standard genetic code: $\pi_{TCA}$ (3) + 19 AA frequencies.
    * `7`: "FMutSel". This setting assigns a fitness to every codon with 60 (= 61 – 1) codon fitness parameters for the standard genetic code. Number of free parameters for standard genetic code: $\pi_{TCA}$ (3) + 60 codon frequencies.
  * `estFreq`: if the frequency/fitness parameters are estimated by maximum likelihood from the data, then `estFreq = 1`. Otherwise, `estFreq = 0` will enable their calculation using the observed frequencies in the data.
  * `fix_omega` and `omega`: depending on the mode assigned using variables `model`, `NSsites`, `CodonFreq`, and `estFreq`, we can decide whether fixing the $\omega$ ratio (`fix_omega = 1`) or estimating it (`fix_omega = 0`). If `fix_omega = 1`, the value of the $\omega$ ratio will fixed to the value passed to variable `omega`. If `fix_omega = 0`, the starting value of $\omega$ during ML estimation will be the value passed to variable `omega`.
  * `clock`: this variable will enable models concerning rate constancy or variation among lineages. When `clock = 0` is specified, there is no clock assumed, and thus rates are entirely free to vary from branch to branch (number of parameters: $2n-3$, where $n$ is the number of species in the sequence alignment and a binary tree is assumed). An unrooted tree should be used under this model. When `clock = 1`, the strict clock model is enabled, and thus all branches have the same rate (number of parameters: $n-1$, corresponding to the "n-1" internal nodes in the binary tree). A test of the molecular clock assumption comparing these two assumptions should have a degree of freedom $d.f = n-2$. You can enable other types of clocks, but they are out of the scope of this protocol (if interested, please consult the [`PAML` Wiki](https://github.com/abacus-gene/paml/wiki/) for more details).

### Input sequence alignment and tree files

Now that you are more familiar with the settings you need to specify in the control file to run `CODEML`, let's focus on our input files: the sequence alignment and the tree file.

During day 1 and day 2, we ran all necessary steps to obtain our codon-aware sequence alignments:

* `aln_nuc_against_protsuper5.fasta`: unaligned nucleotide sequences were aligned against a protein alignment.
* `aln_nuc_against_protusalign.fasta`: unaligned nucleotide sequences were aligned against a structural alignment.

To complete this practical session on time, we shall use only one of those alignments (`aln_nuc_against_protsuper5.fasta`). Nevertheless, you are more than free to repeat this tutorial using the other alignment too!

Nevertheless, we will need to convert the FASTA-formatted sequence files to PHYLIP. There are many tools that can help you convert these alignments, so please feel free to choose the tool you want for this purpose (you can even write your own script if you want to!). To speed things up, you will find the sequence alignment in PHYLIP format inside `day3/inp_data` (see below when we start with the analyses!).

During day 2, we learnt how to infer the best-scoring maximum-likelihood tree wit `IQ-TREE`. Today, we have learnt how to infer a phylogeny under a Bayesian approach with `MrBayes`. When inferring the trees with alignment `aln_nuc_against_protsuper5.fasta`, we obtain the same **unrooted tree topology** (branch lengths may differ):

```txt
((((((Pig_Mx, (Sheep_Mx, Cow_Mx)), Dog_Mx), (Rhesus_macaque_Mx, (Orangutan_Mx, Human_Mx))), (Rat_Mx, Mouse_Mx)), Chimpanzee_Mx), Duck_Mx, Chicken_Mx);
```

If we were to root it by separating birds (Duck and Chicken) from the mammals, we would obtain the following **rooted tree topology** in Newick format:

```txt
((((((Pig_Mx, (Sheep_Mx, Cow_Mx)), Dog_Mx), (Rhesus_macaque_Mx, (Orangutan_Mx, Human_Mx))), (Rat_Mx, Mouse_Mx)), Chimpanzee_Mx), (Duck_Mx, Chicken_Mx));
```

You can clearly see that the inferred topology of this **gene tree** is in conflict with the expected topology of the **species tree** (i.e., "Chimpanzee_Mx" should cluster with the rest of sequences from primates):

<p align="center">
<img width="400" height="350" src="../figs/CODEML_PhyloTree.jpeg">
</p>

> Figure 2. Phylogenetic tree for ten mammals and two bird species, , from Figure 1 in [Álvarez-Carretero et al. 2023](https://doi.org/10.1093/molbev/msad041). All silhouettes were obtained from images dedicated to the public domain, accessible via [Phylopic](https://www.phylopic.org/images/).

Mismatches between gene trees and species trees are not rare, and can be due to different scenarios:

* Compositional biases.
* Long-branch attraction.
* Model assumptions being violated.
* Systematic errors.
* Data filtering.
* And many others!

Under such a conflict, deciding whether fixing the species tree or the gene tree is not easy. The models used in tests for positive selection assume that the fixed tree topology represent the "true" evolutionary relationships of the sequences, and so one should use whichever tree is most likely to be correct. As discussed in section "Gene tree versus species tree" in the [Supplementary Material of Álvarez-Carretero et al. 2023](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/mbe/40/4/10.1093_molbev_msad041/2/msad041_supplementary_data.pdf?Expires=1774651114&Signature=phKHDR~EitWyE~wSuNTBeR1jhwChp3q7E5SFRBv~u9Z1uVP~Q7gUJNs-4vRJxxYdgg7uBidCOFxKOotoLJeJMmFdqTs3AC3VLROgn~mSVgrICce7C6OOE-LqS4Dh3uM9NLhSKEOqtxKMfaTpcgTNXTHW2woP9WHwIP2i-QVV~W9z96ZuDFBg6HhqT3xFyyYmhI2JNDe~L8TGuSLDjI-k2kACz3GBt46-GYSS9dPJ75GexYd-hJ1G-s4Ravm329wSba4UmYdyeNHCV0QB4GKrAx1tIRoE6NOdWMS-4YDJuMdAyq2IdIvOK-xJF2MFFgx5vB5Ip-HJqAWAAlPIFO4RXw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA):

* In analysis of duplicated genes with orthologs and paralogs, the species tree may not be applicable so that the inferred gene tree is the only choice.
* In analysis of viral sequences, a species tree does not exist, and hence the gene tree is the only choice.
* If the gene sequences are short or otherwise do not contain much phylogenetic information, the inferred gene tree may be unresolved or incorrect, and the species tree will be preferable.
* If convergent evolution is likely to have misled gene tree reconstruction, the species tree will be preferable.

In our case, we are analysing a relatively short sequence, which could explain the mismatch we observed when comparing the gene tree with the species tree. This scenario does not necessarily mean that the phylogeny we inferred is wrong: structurally-speaking, the tree may be sensible and may not align with what is expected from the relationships depicted in the species tree. Nevertheless, if further research does not support the latter (e.g., structural trees may can reflect better long-term evolution [deep phylogenies] than short-term evolution [shallow trees]), using the sequence-based species tree could be a better approach.

Given all the discussion above, we will be using the species tree in this tutorial:

> Unrooted, PHYLIP format (header: 12 species, 1 tree in Newick format)

```txt
12  1
((((((Chimpanzee_Mx,Human_Mx),Orangutan_Mx),Rhesus_macaque_Mx),(((Sheep_Mx,Cow_Mx),Pig_Mx),Dog_Mx)),(Mouse_Mx,Rat_Mx)),Duck_Mx,Chicken_Mx);
```

> Rooted, PHYLIP format (header: 12 species, 1 tree in Newick format)

```txt
12  1
((((((Chimpanzee_Mx,Human_Mx),Orangutan_Mx),Rhesus_macaque_Mx),(((Sheep_Mx,Cow_Mx),Pig_Mx),Dog_Mx)),(Mouse_Mx,Rat_Mx)),(Duck_Mx,Chicken_Mx));
```

You will find a copy of these input files (and relevant control files to run `CODEML`) in the `day3` directory that you will now copy to your main working directory!

### Running `CODEML`

In order to start the analyses with `CODEML`, please access the `biol0033` server and type the following from your main working directory `my_session`:

```sh
# Run from "my_session"
# Change directories if you are not
# there yet

# Copy data and scripts for day 3
# and go to dir "day3"
cp -R ~/biol0033-tutorial/day3 .
cd day3
```

If you explore the different directories inside `day3` (either via the command line or the `Files` section on the bottom-right panel), you will see that we have different tree files (you will see which ones we use depending on the `CODEML` analysis) and the PHYLIP-formatted sequence alignment inside `inp_data` and then one directory for analysis: `00_homogeneous_model`, `01_site_models`, `02_branch_models`, `03_branchsite_models`. Let's see how we can run each of them!

#### Homogeneous model

The simplest codon model implemented in `CODEML` is the so-called `M0`, which assumes that  $\omega$ does not vary across sites or across lineages. All alignment sites of a gene are assumed to evolve under the same evolutionary pressure in all taxa.

The control file is as follows:

```txt
      seqfile = ../../inp_data/aln_codawdna.phy   * Path to input alignment file
     treefile = ../../inp_data/Mx_unroot.tree     * Path to input tree file
      outfile = out_M0.txt                        * Path to output file

        noisy = 3     * Moderate information printed on the screen
      verbose = 1     * Detailed report

      seqtype = 1     * Codon data
        ndata = 1     * One alignment block
        icode = 0     * Standard genetic code 
    cleandata = 0     * Do not remove sites with ambiguity data

        model = 0     * ω varying across lineages? No
      NSsites = 0     * ω varying across sites? No
    CodonFreq = 7     * FMutSel model for codon frequencies
      estFreq = 0     * Estimate codon frequencies by ML
    fix_omega = 0     * Estimate omega
        omega = 0.5   * Initial value of omega
        clock = 0     * No clock model
```

You can see that we are using `model = 0` and `NSsites = 0` to enable model `M0`. We want to estimate $\omega$, and so we use `fix_omega = 0` with a starting value of `omega = 0.5`. We will not assume the molecular clock (`clock = 0`) because the specified codon model is time-reversible. We will therefore use an unrooted tree (see `treefile`). All paths to input/output files are either relative paths (`seqfile` and `treefile`) or will be saved in the same directory as the control file (`outfile`).

To run `CODEML` with these settings, please do the following:

```sh
# Run from "my_session/day3"
# Change directories if you are not
# there yet

# Go to directory with control file to
# run M0 with CODEML
cd 00_homogeneous_model/Model_M0
# Run CODEML but also save a copy of the 
# screen output in a log file
codeml codeml-M0.ctl >(tee logfile_codemlM0.txt >&2)
```

As soon as you start the run, you will see that various intermediate files alongside output file `out_M0.txt` start being generated inside directory `Model_M0` while lots of information about this analysis is printed on the screen. Once `CODEML` finishes, you are ready to look at output file `out_M0.txt`:

> **Summary of site patterns in the input sequence alignment**

The input sequence alignment and its compressed version showing only site patterns are printed at the top of the output file. The latter will look like the following:

```txt
Printing out site pattern counts


        12       2229  P

Duck_Mx                          --- --- --- --- --- --- --- --- --- --- --- [...]
Chimpanzee_Mx                    --- --- --- --- --- --- --- --- --- --- --- [...]
Chicken_Mx                       --- --- --- --- --- --- --- AAT AGC ATG CAG [...]
Orangutan_Mx                     --- --- --- --- --- --- TCG --- --- --- --- [...]
Human_Mx                         --- --- --- --- --- --- TCG --- --- --- --- [...]
Pig_Mx                           --- --- --- --- --- GAA --- --- --- --- --- [...]
Rhesus_macaque_Mx                --- --- --- --- --- --- --- --- --- --- --- [...]
Dog_Mx                           --- --- --- --- --- --- --- --- --- --- --- [...]
Sheep_Mx                         --- --- --- --- --- --- --- --- --- --- --- [...]
Cow_Mx                           --- --- --- --- --- --- --- --- --- --- --- [...]
Rat_Mx                           --- --- CCT CTA GGT --- TCT --- --- --- --- [...]
Mouse_Mx                         CAG TTT CCT --- --- --- TCT --- --- --- --- [...]



    1    1    1    1    1    1    1    1    1    1    1    1    1    2    1
    1    1    1    1    1    1    1    1    1    1    1    1    1    1    2
    1    1    1    1    1    1    1    1    1    1    2    1    1    3    2
    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
    1    1    1    1    1    1    1    1    1    1    1    1    4    1    1
    1    1    1    1    1    1    3    1    1    1    1    1    1    1    1
    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
```

You can see that, in the snapshot above taken from the `out_M0.txt` file, most site patterns are repeated once, but there are a few site patterns that are repeated twice, three times, and one site pattern is repeated four times. You can explore the full output file to check the rest of site patterns.

> **Summary of the input sequence alignment and the model selected**

You will then see the version of `PAML` you are using followed by the path to the sequence alignment file. Then, you will see details about the model specified in the control file: model `M0` (one dN/dS ratio) and "FMutSel" to estimate codon frequencies. Then, you will see the number of species (`ns = 12`) and the number of codons (`ls = 799`; number of base pairs in the alignment [2,397] divided into 3):

```txt
CODONML (in paml version 4.10.10, 27 Jan 2026)  ../../inp_data/aln_codawdna.phy
Model: One dN/dS ratio, 
Codon frequency model: FMutSel
ns =  12  ls = 799
```

> **Summary of nucleotide and codon frequencies**

The next blocks correspond to observed nucleotide and codon frequencies for each sequence and their average over all sequences followed by the codon frequencies under the model. The latter can be used for simulation purposes (e.g., you can use the `PAML` program `evolver`).

> **Summary of tree scores and estimated model parameter values**

Next, you will see a matrix with estimates of dS and dN using the method of [Nei and Gojobori (1986)](https://doi.org/10.1093/oxfordjournals.molbev.a040410). According to the `PAML` documentation, these values are used to construct initial estimates of branch lengths for the subsequent likelihood-based analysis, but they are not MLEs themselves. Consequently, this matrix is printed just for reference (this matrix can be ignored).

The next section shows a relevant summary of our analysis:

```txt
TREE #  1:  ((((((2, 5), 4), 7), (((9, 10), 6), 8)), (12, 11)), 1, 3);   MP score: -1
lnL(ntime: 21  np: 26): -14896.436992      +0.000000
  13..14   14..15   15..16   16..17   17..18   18..2    18..5    17..4    16..7    15..19   19..20   20..21   21..9    21..10   20..6    19..8    14..22   22..12   22..11   13..1    13..3  
 1.498321 0.258311 0.265745 0.051382 0.009002 1.465875 0.031683 0.025348 0.062455 0.062712 0.186517 0.276587 0.064576 0.053741 0.237447 0.411994 0.456193 0.250451 0.241940 0.317898 0.536876 2.193232 0.697237 1.364936 1.729643 0.329189

Note: Branch length is defined as number of nucleotide substitutions per codon (not per nucleotide site).

tree length =  6.765056

((((((2: 1.465875, 5: 0.031683): 0.009002, 4: 0.025348): 0.051382, 7: 0.062455): 0.265745, (((9: 0.064576, 10: 0.053741): 0.276587, 6: 0.237447): 0.186517, 8: 0.411994): 0.062712): 0.258311, (12: 0.250451, 11: 0.241940): 0.456193): 1.498321, 1: 0.317898, 3: 0.536876);

((((((Chimpanzee_Mx: 1.465875, Human_Mx: 0.031683): 0.009002, Orangutan_Mx: 0.025348): 0.051382, Rhesus_macaque_Mx: 0.062455): 0.265745, (((Sheep_Mx: 0.064576, Cow_Mx: 0.053741): 0.276587, Pig_Mx: 0.237447): 0.186517, Dog_Mx: 0.411994): 0.062712): 0.258311, (Mouse_Mx: 0.250451, Rat_Mx: 0.241940): 0.456193): 1.498321, Duck_Mx: 0.317898, Chicken_Mx: 0.536876);

Detailed output identifying parameters

kappa (ts/tv) =  2.19323

Frequency parameters:
   0.14551 (T)   0.28485 (C)   0.36096 (A)   0.20869 (G)
omega (dN/dS) =  0.32919

dN & dS for each branch

 branch          t       N       S   dN/dS      dN      dS  N*dN  S*dS

  13..14     1.498  1784.8   612.2  0.3292  0.3285  0.9979 586.3 610.9
  14..15     0.258  1784.8   612.2  0.3292  0.0566  0.1720 101.1 105.3
  15..16     0.266  1784.8   612.2  0.3292  0.0583  0.1770 104.0 108.3
  16..17     0.051  1784.8   612.2  0.3292  0.0113  0.0342  20.1  20.9
  17..18     0.009  1784.8   612.2  0.3292  0.0020  0.0060   3.5   3.7
  18..2      1.466  1784.8   612.2  0.3292  0.3214  0.9763 573.6 597.6
  18..5      0.032  1784.8   612.2  0.3292  0.0069  0.0211  12.4  12.9

[...]

tree length for dN:       1.4831
tree length for dS:       4.5054
```

Remember that we used an unrooted tree with $n = 12$ species, which has $2\times n-3$ branch lengths (i.e., $2\times 12-3= 21$ branch lengths). Free parameters include branch lengths, the equilibrium frequencies, the transition/transversion rate ratio ($\kappa$), and the parameters in the omega distribution. It may be difficult to match the MLEs in line 4 to the corresponding branch lengths in the phylogeny or to the model parameters, but the branch lengths are printed again in the Newick tree. The last block shows the estimates of $\kappa$ (transition/transversion rate ratio) and $\omega$, as well as the mutation bias (nucleotide frequency) parameters (i.e., line `0.14551 (T)   0.28485 (C)   0.36096 (A)   0.20869 (G)`). The mutation-selection model ("FMutSel") accommodates different codon frequencies by modelling mutational biases and fixations of mutations under selection ([Yang and Nielsen 2008](https://doi.org/10.1093/molbev/msm284), eqs. 1 and 4). In this analysis, it looks like the estimated frequency for "A" is higher than the rest, and thus the mutation process is biased toward A. Estimated values of $t$ (branch lengths), $d_{N}$ (nonsynonymous rate), and $d_{S}$ (synonymous rate) for each branch follow.

Since the homogenous model (`M0`) is specified (`model = 0` and `NSsites = 0`), the same $\omega$ value (**0.3292**) is reported for each branch under the column `dN/dS`. The tree length is $d_{N} = 1.4831$ at the nonsynonymous sites and $d_{S} = 4.5054$ at the synonymous sites, with $\omega = d_{N}/d_{S} = 0.3292$, suggesting that the myxovirus gene is overall under purifying selection, with a nonsynonymous mutation having on average a third as large a chance (i.e., $\omega = d_{N}/d_{S} = 0.3292$) of going to fixation as a synonymous mutation.

This is all the information you will find in the output file! Now, let's remove unnecessary intermediate files that we do not need to use and extract the estimated $\omega$ ratio values in an output file for subsequent analyses when comparing other models for tests of positive selection:

```sh
# Run from "00_homogeneous_model/Model_M0"
# Change directories if you are not
# there yet

# Delete unnecessary files
rm 2N*

# Extract omega vals
printf "omega\n" > omega_est.txt
grep 'omega ' out_M0.txt | sed 's/..*= *//' >> omega_est.txt
```

We will now analyse our data under other codon models useful for detecting positive selection: site models, branch models, and branch-site models! :muscle:

#### Site models

Now, we will run `CODEML` under various **site models**, which allow $\omega$ to vary across codons. Instead of running `CODEML` once for every different site model, there is a feature that you can enable through variable `NSsites` that will allow you to run a "batch" analysis. In other words, you can instruct `CODEML` to run more than one site model by typing the options that enable the corresponding site model to be enabled (i.e., `CODEML` will run the second analysis once after the first one has finished, the third after the second has finished, etc.).

We will be running a batch analysis under the `M0` model (this is the homogeneous model with one $\omega$ ratio, we have already run an analysis under this model!) and then the following site models: `M1a` (nearly neutral), `M2a` (positive selection), `M7` (beta) and `M8` (beta&$\omega$). You can find a summary of these site models in the table below, including the model comparisons you need to carry out with a **LRT** (the models being compared are **nested**!):

> Table 1. Site models for testing positive selection affecting amino acid residues in a protein (from Table 1 in in [Álvarez-Carretero et al. 2023](https://doi.org/10.1093/molbev/msad041))

<table>

<!-- HEADER -->
<tr>
<th>Site model</th>
<th>Free parameters</th>
<th>Number of site classes</th>
<th>Model parameters</th>
<th>Model comparison (LRT)</th>
</tr>

<!-- FIRST ROW -->
<tr>
<td>M1a</td>
<td>2 (p<sub>0</sub>, &omega;<sub>0</sub>)</td>
<td>2</td>
<td>p<sub>0</sub> (p<sub>1</sub> = 1 − p<sub>0</sub>),
&omega;<sub>0</sub> &lt; 1,
&omega;<sub>1</sub> = 1</td>
<td  rowspan="2">M1a vs. M2a: df = 2</td>
</tr>
<!-- SECOND ROW -->
<tr>
<td>M2a</td>
<td>4 (p<sub>0</sub>, p<sub>1</sub>, &omega;<sub>0</sub>, &omega;<sub>2</sub>)</td>
<td>3</td>
<td>p<sub>0</sub>,
p<sub>1</sub> (p<sub>2</sub> = 1 − p<sub>0</sub> − p<sub>1</sub>),
&omega;<sub>0</sub> &lt; 1,
&omega;<sub>1</sub> = 1,
&omega;<sub>2</sub> &gt; 1</td>
</tr>
<!-- THIRD ROW -->
<tr>
<td>M7</td>
<td>2 (p, q)</td>
<td>10</td>
<td>beta(p, q)</td>
<td  rowspan="2">M7 vs. M8: df = 2</td>
</tr>
<!-- FOURTH ROW -->
<tr>
<td>M7</td>
<td>4 (p<sub>0</sub>, p, q, &omega;<sub>s</sub>)</td>
<td>11</td>
<td>p<sub>0</sub>,
(p<sub>1</sub> = 1 − p<sub>0</sub>),
beta(p, q),
&omega;<sub>s</sub> &gt; 1</td>
</tr>
</table>

> [!IMPORTANT]
> The number of free parameters in the $\omega$ distribution will vary depending on the site model under which data are being analysed:
>
> * `M1a`: $p_{0}$ is the proportion of sites with $\omega_{0} < 1$, while $p_{1} = 1 − p_{0}$ is the proportion of sites with $\omega_{1} = 1$.
> * `M2a`: same as `M1a` except for the fact that `M2a` includes an additional class of sites with $\omega_{2} > 1$ in proportion $p_{2}$, with $p_{0} + p_{1} + p_{2} = 1$. The first test for positive selection that you will run compares model `M1a` (nearly neutral) against the model that accounts for positive selection include an extra class to allow for these sites being considered (i.e.,  $\omega_{2} > 1$ in proportion $p_{2}$): `M2a`.
> `M7`: this model uses a beta distribution with parameters $p$ and $q$ to describe variable $\omega$ for sites in the range $0 ≤ \omega ≤ 1$.
> `M8`: $p_{0}$ is the proportion of sites with $\omega$ from $beta(p, q)$ as in `M7`, but now an additional class is added (with proportion $p_{1}$) with $\omega_{s} > 1$. By comparing models `M7` (null) and `M8` (alternative), you can run another test for positive selection. Nevertheless, please note that this test is less stringent than the "M1a-M2a test" described above.

<!-- Table in markdown format, left here to make LaTeX format has been properly translated into HTML above

| Site model | Free parameters   | Number of site classes | Model parameters                                  | Model comparison (LRT) |
|------------|--------------------|------------------------|---------------------------------------------------|------------------------|
| M1a        | 2 ($p_{0}$, $\omega_{0}$)         | 2                      | $p_{0}$ $(p_{1} = 1 − p_{0})$, $\omega_{0} < 1$, $\omega_{1} = 1$                  | M1a vs. M2a: df = 2    |
| M2a        | 4 ($p_{0}$, $p_{1}$, $\omega_{0}$, $\omega_{2}$) | 3                      | $p_{0}$, $p_{1} (p_{2} = 1 − p0 − p1)$, $\omega_{0} < 1$, $\omega_{1} = 1$, $\omega_{2} > 1$ |                        |
| M7         | 2 ($p$, $q$)           | 10                     | $beta(p, q)$                                        | M7 vs. M8: df = 2      |
| M8         | 4 ($p_{0}$, $p$, $q$, $\omega_{s}$)   | 11                     | $p_{0}$, $(p_{1} = 1 − p0)$, $beta(p, q)$, $\omega_{s} > 1$             |  M7 vs. M8: df = 2                        |

-->

You will see that, under `day3/01_site_models/Site_models`, you have the control file to launch a batch analysis that will analyse your data under the homogeneous (`M0` model) and four additional heterogeneous models (site models, `M1a`, `M2a`, `M7`, and `M8`). The differences between this control file and the control file to run the `M0` model are the following:

```txt
outfile = out_sites.txt  * Path to output file
NSsites = 0 1 2 7 8      * ω varying across sites? Yes
```

Note that the name of the output file is now `out_sites.txt` and we are enabling the batch model via variable `NSsites`. Run the commands in the code snippet below to launch `CODEML`:

```sh
# Run from "01_site_models/Site_models"
# Change directories if you are not
# there yet

# Run CODEML but also save a copy of the 
# screen output in a log file
codeml codeml-sites.ctl >(tee logfile_codeml-sites.txt >&2)
```

The output file follows the same format as described above for the homogeneous model. The only difference is that now we will have one block of summarising tree scores and model parameters for each of the models under which the data were analysed (i.e., 5 blocks).

As you learnt during the theoretical session focused on phylogenetic reconstruction under maximum likelihood (and as you revised in the latest session on molecular adaptation), we can use the **LRT statistic** to compare **nested models**. The LRT statistic is defined as **twice the difference in log-likelihood between the null and alternative hypotheses**, $2\Deltaℓ = 2(ℓ_{1} − ℓ_{0})$, where $ℓ_{0}$ is the log-likelihood score for the null model and $ℓ_{1}$ is the log-likelihood under the alternative model. The LRT statistic is compared with the **$\chi_{2}$ distribution** with the **degree of freedom (d.f.)** equal to the **difference in the number of free parameters** (which `CODEML` prints in the output file for each model) between the two models being compared. The model comparisons we will carry out are the following:

* `M0` vs. `M1a` (one-ratio vs. nearly neutral): this is a **test for variable of selective pressure among amino acid sites** rather than a test of positive selection. We will only proceed to carry out the next test to detect for positive selection (i.e., `M1a` vs `M2a`) if we can reject the null hypothesis (`M0`) that the one-ratio model fits the data better than a site model allowing $\omega$ to vary among sites (i.e., `M1a`, nearly-neutral model).
* `M1a` vs. `M2a` (nearly neutral vs. positive selection): this is a **test for positive selection** given that `M2a` adds a class to consider sites under positive selection: $\omega_{2} > 1$, with a proportion of sites of $p_{2}$.
* `M7` vs. `M8` (beta vs. beta&$\omega$): this another **test for positive selection** that we will run (null: `M7`, alternative: `M8`), but it is not as stringent as the "M1a-M2a" test described above.

In order to carry out the LRTs for each model comparison, we will need the log-likelihood values. Once our batch analysis finishes, you can extract these values from the output file using the code snippet below:

```sh
# Run from "01_site_models/Site_models"
# Change directories if you are not
# there yet

# First, remove unnecessary files 
rm 2N*

# Then, extract the likelihood values to then
# calculate LRT statistic
lnL_vals=$( grep 'lnL' out_sites.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' )
np_vals=$( grep 'lnL' out_sites.txt | sed 's/..*np\:\ //' | sed 's/)..*//' )
header=$( grep 'Model ' out_sites.txt | sed 's/\:..*//' | sed 's/ Model /_/' )
echo $header > lnL_sites.txt
echo $lnL_vals >> lnL_sites.txt
echo $np_vals >> lnL_sites.txt
```

Now, you are ready to run the R script that we have prepared for you to run all model comparisons! In order to run `Find_best_site_model.R`, you will need to do the following:

* Go to the bottom right panel and navigate to `my_session/day3/01_site_models`, then click `Find_best_site_model.R` to open the R script on your top left panel.
* Go through the R script and ask questions about commands you do not understand. Basically, this script will do the following general tasks:
  * Set your working directory
  * Read the `lnL_sites.txt` file that you have just created with the log-likelihood values and model parameters for each of the models you have just run in this batch analysis.
  * Calculate the LRT statistic and the alpha-critical values at 5% and 1% for each test considering the degree of freedom (i.e., this is calculated based on the difference of model parameters!).
  * Plot the LRT statistic, the alpha-critical values at 5% and 1%, and the p-value for the $\chi^2$ test with the calculated degree of freedom -- there is one plot per model comparison.

When you obtain the final plot, you can click `PLots > Export > Save as PDF...` on the bottom right panel, and then choose "Landscape" and directory `~/my_session/day3/01_site_models` in the last pop-up window. You can save the output file as `LRT_site_models.pdf`. If you want to see the plot after saving it, then tick the box `View plot after saving`:

<p align="center">
<img width="400" height="300" src="../figs/CODEML_LRT_sitemods_plot.png">
</p>

> **[ EXERCISE ]**
> Open a text editor or any other program you feel comfortable with and write two tables summarising the following:
>
> * Table 1: table summarising the log-likelihood values for each of the models being compared, the corresponding free parameters, the degree of freedom, and the LRT statistic. The header could be the following: "Model comparison", "Log-likelihood values, $ℓ_{0}$ and $ℓ_{1}$", "Free parameters", "d.f", "LRT statistic".
> * Table 2: table summarising the log-likelihood values and parameter estimates for each model under which data were analysed in this batch analysis. The header could be the following: "Model", "Log-likelihood, $ℓ$", "$d_{N}/d_{s}", and "Estimates of model parameters".
>
> Then, answer the following questions:
>
> * Which is/are the best-fitting model/s according to the LRTs you have carried out?
> * How would you interpret the model parameters estimated under each model?

> [!TIP]
> You can use bash scripting to quickly parse the output files when searching for specific lines that would have information about estimates of model parameters:
> Example 1: which specific line under the output block for each model can we use to extract information about the $\omega$ ratio? Well, the line starting with " branch" is then followed by two lines that allow us to see the first value of estimated $\omega$ ratio. The estimated $d_{N}/d_{S} value is the same for each branch, thus we only need to see the first line. We can use the following command to find the five instances that match this pattern using command `grep` and option `-A2` (prints two lines after the line that matches the pattern to be found in the file):
`grep '^ branch' -A2 out_sites.txt`:
> 
> ```txt
>  branch          t       N       S   dN/dS      dN      dS  N*dN  S*dS
>
>  13..14     1.498  1784.8   612.2  0.3292  0.3285  0.9979 586.3 610.9
>--
> branch          t       N       S   dN/dS      dN      dS  N*dN  S*dS
>
>  13..14     2.235  1753.6   643.4  0.5038  0.5892  1.1696 1033.3 752.5
>--
> branch          t       N       S   dN/dS      dN      dS  N*dN  S*dS
>
>  13..14     2.235  1753.6   643.4  0.5038  0.5892  1.1696 1033.3 752.5
>--
> branch          t       N       S   dN/dS      dN      dS  N*dN  S*dS
>
>  13..14     2.308  1767.3   629.7  0.3924  0.5468  1.3934 966.3 877.4
>--
> branch          t       N       S   dN/dS      dN      dS  N*dN  S*dS
>
>  13..14     2.505  1759.0   638.0  0.4858  0.6515  1.3412 1146.0 855.6
> ```
>
> Each time the pattern has been found is separated with `--`!
>
> Now, your turn! How would you find model parameters for the site proportions and the estimated $\omega$ ratios for each site class?

<br>
<details>
<summary><b>[ Click here to find the command that you can use to find model parameters ]</b></summary>
<br>
The solution is <code>grep 'Frequency parameters' -A6 out_sites.txt</code>. While some models will have blank lines or lines with other info after printing the third or fourth lines after the pattern is found (<code>Frequency parameters</code>), models <code>M1a</code> and <code>M2a</code> require 6 lines to be printed after the pattern to see all parameter estimates:

```txt
Frequency parameters:
   0.14551 (T)   0.28485 (C)   0.36096 (A)   0.20869 (G)
omega (dN/dS) =  0.32919

dN & dS for each branch

 branch          t       N       S   dN/dS      dN      dS  N*dN  S*dS
--
Frequency parameters:
   0.14647 (T)   0.33716 (C)   0.28715 (A)   0.22922 (G)

MLEs of dN/dS (w) for site classes (K=2)

p:   0.54306  0.45694
w:   0.08631  1.00000
--
Frequency parameters:
   0.14647 (T)   0.33716 (C)   0.28715 (A)   0.22922 (G)

MLEs of dN/dS (w) for site classes (K=3)

p:   0.54306  0.19572  0.26122
w:   0.08631  1.00000  1.00000
--
Frequency parameters:
   0.14791 (T)   0.30733 (C)   0.31722 (A)   0.22753 (G)
Parameters in M7 (beta):
 p =   0.33803  q =   0.52339


MLEs of dN/dS (w) for site classes (K=10)
--
Frequency parameters:
   0.14578 (T)   0.31510 (C)   0.30783 (A)   0.23128 (G)
Parameters in M8 (beta&w>1):
  p0 =   0.91379  p =   0.39466 q =   0.77341
 (p1 =   0.08621) w =   2.05744
```

</details><br>

<details>
<summary><b>[ Click here only when you have finished compiling your first summary table to check its content ]</b></summary>
<br>

<table>

<!-- HEADER -->
<tr>
<th>Model comparison</th>
<th>Log-likelihood values (ℓ<sub>0</sub> and ℓ<sub>1)</sub></th>
<th>Free parameters</th>
<th>d.f</th>
<th>LRT statistic<br>
(2&Delta;ℓ)</th>
</tr>

<!-- FIRST ROW -->
<tr>
<td>M0 vs. M1a (one-ratio vs. nearly neutral)</td>
<td>ℓ<sub>0</sub> = −14,896.44<br>
ℓ<sub>1</sub> = −14,440.41</td>
<td>26 vs. 27</td>
<td>1</td>
<td>912.06</td>
</tr>
<!-- SECOND ROW -->
<tr>
<td>M1a vs. M2a (nearly neutral vs. positive selection)</td>
<td>ℓ<sub>0</sub> = −14,440.41<br>
ℓ<sub>1</sub> = −14,440.41</td>
<td>27 vs. 29</td>
<td>2</td>
<td>0</td>
</tr>
<!-- THIRD ROW -->
<tr>
<td>M7 vs. M8 (beta vs. beta&&omega;)</td>
<td>ℓ<sub>0</sub> = −14,400.51<br>
ℓ<sub>1</sub> = −14,380.93</td>
<td>27 vs. 29</td>
<td>2</td>
<td>39.15</td>
</tr>

</table>

</details><br>

<details>
<summary><b>[ Click here only when you have finished compiling your second summary table to check its content ]</b></summary>
<br>

<table>

<!-- HEADER -->
<tr>
<th>Model</th>
<th>Log-likelihood, ℓ</th>
<th>d<sub>N</sub>/d<sub>s</sub></th>
<th>Estimates of model parameters</th>
</tr>

<!-- FIRST ROW -->
<tr>
<td>M0 (one-ratio)</td>
<td>ℓ = −14,896.44</td>
<td>0.3292</td>
<td>&omega; = 0.32919</td>
</tr>
<!-- SECOND ROW -->
<tr>
<td>M2a (positive selection)</td>
<td>ℓ = −14,440.41</td>
<td>0.5038</td>
<td>p<sub>0</sub> = 0.54306 (p<sub>1</sub> = 0.45694)<br>
&omega;<sub>0</sub> = 0.08631</td>
</tr>
<!-- THIRD ROW -->
<tr>
<td>M2a (positive selection)</td>
<td>ℓ = −14,440.41</td>
<td>0.5038</td>
<td>p<sub>0</sub> = 0.54306, p<sub>1</sub> = 0.19572 (p<sub>2</sub> = 0.26122)<br>
&omega;<sub>0</sub> = 0.08631</td>
</tr>
<!-- FOURTH ROW -->
<tr>
<td>M7 (beta )</td>
<td>ℓ = −14,400.51</td>
<td>0.3924</td>
<td>p = 0.33803, q = 0.52339</td>
</tr>
<!-- FIFTH ROW -->
<tr>
<td>M8 (beta&&omega;)</td>
<td>ℓ = −14,380.93</td>
<td>0.4858</td>
<td>p<sub>0</sub> = 0.91379 (p<sub>1</sub> = 0.08621)<br>
p = 0.39466, q = 0.77341<br>
&omega; = 2.05744</td>
</tr>

</table>

</details><br>

<details>
<summary><b>[ Click here only when you have finished thinking about the results and proposed an interpretation to those ]</b></summary>
<br>

</details><br>

The last output we will pay attention is the **Bayes Empirical Bayes (BEB) method**, used to calculate the posterior probability for each site coming from the different site classes. This approach is printed on the screen when analysing the dataset under `M2a` and `M8` models. Sites with high posterior probabilities for the positively selected class are likely to be under positive selection.

> [!IMPORTANT]
> **The results for the BEB method will be always printed on the screen and calculated only under models of positive selection, not under the null models**. `CODEML` does not compute any LRT and does not know whether `M2a` or `M8` fit the data better than the null models they are being compared to. Once you have carried out your LRTs, if you cannot reject the null hypotheses (i.e., `M1a` in "M1-M2a test" or `M7` in "M7-M8 test"), **the BEB results should be ignored as they refer to a model that does not fit the data better than the null model**. In addition, please ignore all sections in the output file referring to the Naïve Empirical Bayes (NEB) method as it is an older method that does not accommodate the uncertainties in the MLEs that the BEB method does ([Yang et al. 2005](https://doi.org/10.1093/molbev/msi097)).

In our site model analysis, we cannot look at the BEB results printed under the output block for model `M2a` because the null model `M1a` cannot be rejected. Nevertheless, we can look at the results under `M8` as we can reject the null model `M7`:

```txt
Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)
Positively selected sites (*: P>95%; **: P>99%)
(amino acids refer to 1st sequence: Duck_Mx)

            Pr(w>1)     post mean +- SE for w

     3 T      0.946         2.430 +- 0.486
     5 R      0.725         2.063 +- 0.779
     6 N      0.968*        2.465 +- 0.429
     7 T      0.786         2.165 +- 0.727
    64 D      0.624         1.891 +- 0.850
    66 P      0.634         1.906 +- 0.839
    70 P      0.567         1.790 +- 0.870
    78 N      0.625         1.894 +- 0.842
    79 M      0.528         1.725 +- 0.831
    83 N      0.794         2.178 +- 0.724
    84 P      0.532         1.733 +- 0.854
   156 V      0.509         1.689 +- 0.871
   237 P      0.796         2.181 +- 0.716
   252 I      0.793         2.178 +- 0.722
   315 Q      0.596         1.844 +- 0.830
   452 N      0.926         2.397 +- 0.529
   511 R      0.708         2.029 +- 0.779
   518 H      0.852         2.274 +- 0.651
   559 W      0.638         1.913 +- 0.824
   607 E      0.680         1.985 +- 0.801
   692 A      0.772         2.142 +- 0.752
   693 S      0.878         2.319 +- 0.617
   696 S      0.924         2.393 +- 0.535
   697 D      0.963*        2.456 +- 0.444
   702 K      0.887         2.333 +- 0.599
   708 Q      0.816         2.214 +- 0.694
```

In each line, the first column shows the site position (e.g., 3, 5, 6, 7, etc.), which is followed by the amino acid at this site in the first sequence (this is for identification of the site in the sequence which, in this case, is "Duck_Mx"). The third column (`Pr (w > 1)`) shows the posterior probability for the site to be from the positive-selection class (i.e., with $\omega > 1$). The last columns show the posterior mean of the $\omega$ ratio and the standard deviation in the $\omega$ distribution for the site. When analysing our dataset, we see 26 sites that have a probability >50% for the positive-selection class with $\omega > 1$. For instance, site 3 (with amino acid T (Threonine) in the first sequence) has probability 0.946 of coming from the positive-selection class, and the posterior distribution of $\omega$ for the site has mean 2.430 and SD 0.486. Note that there are two sites which probability of coming from the positive-selection class is larger than 0.95, which is indicated with an asterisk after the probability in the second column:

```txt
     6 N      0.968*        2.465 +- 0.429
   697 D      0.963*        2.456 +- 0.444
```

Site 6 in the alignment has a posterior probability of 96.8% of coming from the positive-selection class with $\omega > 1$ (site 697 of 96.3%). The approximate posterior distribution of $\omega$ for site 6 has mean 2.465 and SD 0.429 (2.456 +- 0.444 for site 697).

> [!NOTE]
> In the `CODEML` output, posterior probabilities $P > 0.95$ are indicated with an asterisk, while those with $P > 0.99$ are indicated with two asterisks. You can look at section [“BEB Analysis” Supplementary Material of Álvarez-Carretero et al. 2023](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/mbe/40/4/10.1093_molbev_msad041/2/msad041_supplementary_data.pdf?Expires=1774651114&Signature=phKHDR~EitWyE~wSuNTBeR1jhwChp3q7E5SFRBv~u9Z1uVP~Q7gUJNs-4vRJxxYdgg7uBidCOFxKOotoLJeJMmFdqTs3AC3VLROgn~mSVgrICce7C6OOE-LqS4Dh3uM9NLhSKEOqtxKMfaTpcgTNXTHW2woP9WHwIP2i-QVV~W9z96ZuDFBg6HhqT3xFyyYmhI2JNDe~L8TGuSLDjI-k2kACz3GBt46-GYSS9dPJ75GexYd-hJ1G-s4Ravm329wSba4UmYdyeNHCV0QB4GKrAx1tIRoE6NOdWMS-4YDJuMdAyq2IdIvOK-xJF2MFFgx5vB5Ip-HJqAWAAlPIFO4RXw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA) for another example.

Now, we need to remember what happened during the LRTs. The test for positive selection when comparing "M1a-M2a" Was not significant, while the "M7-M8" comparison resulted in the null hypothesis (`M7`) being rejected and model `M8` fitting the data better at both 5% and 1% significant levels. Considering that the first test is more stringent, if only a few sites are under weak positive selection, then the "M1a-M2a" comparison might fail to reach significance, while the "M7-M8" comparison may detect it. In our case, the estimated $\omega$ ratio under `M8` is quite low ($\omega = 2.05744$), and so that could explain the abovementioned. Nevertheless, we can verify this evidence by further doing some of these checks (there are more!):

* Running `CODEML` with different starting values for $\omega$ to make sure that the MLEs converge to the same estimated values.
* Check the sites that have been identified with a probability larger than 95% under the BEB method of coming from the positive-selection class with $\omega > 1$. Structural analyses may reveal whether there is some structural functionality related to these sites that must have been favoured by a selective pressure.

Overall, we could say that there seems to be some evidence for sites under positive selection in the Mx gene, although this evidence could be further verified.

#### Branch models

#### Branch-site models
