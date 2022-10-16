# uniprot_description_analysis
A tool to analyze significant differences in the frequencies of words found in the UniProt descriptions of set of proteins

The aim of this project is to characterize subsets of proteins identified through proteomics experiments, for example, a set of proteins that are enriched by a chemical probe or are significantly dysregulated upon treatment. Using the UniProt functional descriptions may capture details from the literature about these proteins that are not reflected in e.g. GO annotations. This script therefore takes a list of UniProt accessions for the protein subset and for a 'control' set. This 'control' set could be the proteome as a whole, the set of proteins in the proteome that are not in the subset of interest, or all detected proteins with the exception of the proteins found in the subset of interest. For each of these sets, the script extracts descriptions from a UniProt database, strips them of extraneous annotations. The descriptions are then tokenized and the frequencies of each token recorded. The significance of the difference in word frequencies between the two sets of proteins are then calculated using a log-likelihood approach. A usage example mockup is shown below in which apoptosis-related proteins were used as the test set, and a random selection of one thousand proteins was used as a control:

import uniprot_description_analysis as uda
analyzer = uda.Analyzer('UniProt 20220516.dat')
subset = ['P42575', 'Q14790', 'P55211', 'Q92851', 'P42574', 'P55212', 'P55212', 'P55210', 'O14727', 'P05067', 'P10415', 'Q9UMX3', 'Q9HD36', 'Q07820', 'P04637', 'O15350', 'P17066', 'P0DMV8', 'P07900']
control = np.random.choice(list(analyzer._uniprot_descriptions.keys()), 1000)
result = analyzer.generate_significance_table(subset, control)
result.sort_values('G2', ascending=False)

This script produces the output:

            Subset occurances  Control occurances         G2             p
cleaves                    17                  16  79.770437  0.000000e+00
apoptosis                  25                  99  63.797983  1.332268e-15
programmed                  7                   2  41.899085  9.610757e-11
chaperone                  10                  15  40.455263  2.011697e-10
cycles                      6                   2  35.103549  3.126305e-09
                      ...                 ...        ...           ...

And so this approach appears to work, at least for this contrived example.
