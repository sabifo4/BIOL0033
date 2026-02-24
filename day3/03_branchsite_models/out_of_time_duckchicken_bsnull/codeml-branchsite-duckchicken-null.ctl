      seqfile = ../../inp_data/aln_codawdna.phy            * Path to input alignment file
     treefile = ../../inp_data/Mx_branch_duckchicken.tree  * Path to input tree file
      outfile = out_branchsite_duckchicken_null.txt        * Path to output file

        noisy = 3     * Moderate information printed on the screen
      verbose = 1     * Detailed report

      seqtype = 1     * Codon data
        ndata = 1     * One alignment block
        icode = 0     * Standard genetic code 
    cleandata = 0     * Do not remove sites with ambiguity data

        model = 2     * ω varying across lineages? Yes
      NSsites = 2     * ω varying across sites? Yes
    CodonFreq = 7     * FMutSel model for codon frequencies
      estFreq = 0     * Estimate codon frequencies by ML
    fix_omega = 1     * Estimate omega
        omega = 1     * Initial value of omega
        clock = 0     * No clock model