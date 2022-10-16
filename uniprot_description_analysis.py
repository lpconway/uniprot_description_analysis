# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 11:53:20 2022

@author: Louis
"""

import re
import pandas as pd
from nltk.tokenize import word_tokenize
import numpy as np
from scipy.stats import chi2

class Analyzer:    
    def _remove_pmed_refs(self, description):
        # Remove pubmed references from the descriptions
        result = re.search('\(PubMed:\d+(, PubMed:\d+)*\)', description)
        if type(result) == re.Match:
            limits = result.span()
            description = description[:limits[0]-1] + description[limits[1]:]
            description = self._remove_pmed_refs(description)
        return(description)
    
    def _process_description(self, description):
        # Processes the description strings to remove any unnecessary annotations
        # that might interfere with tokenization, i.e. the references
        description = self._remove_pmed_refs(description)
        result = re.search('{.*}\.', description)
        if type(result) == re.Match:
            limits = result.span()
            description = description[:limits[0]-1] + description[limits[1]:]
            description = self._process_description(description)
        return(description)
    
    def _dehyphenate_token(self, token):
        # Removes common hypenations that add little information so that the 
        # stem is added to the count of the parent word
        suffixes = [
            '-bound',
            '-dependent',
            '-mediated',
            '-associcated',
            '-regulated',
            '-induced',
            '-binding',
            '-specific',
            '-like',
            '-linked',
            '-containing',
            '-activated',
            '-regulation',
            '-related',
            '-coupled',
            '-stimulated',
            '-regulates',
            '-activating',
            '-derived',
            '-regulated',
            '-apoptotic', # Note
            '-based',
            '-positive',
            '-forming',
            '-selective',
            '-conjugating',
            '-recognition',
            '-inflammatory', #
            '-recognition',
            '-linking',
            '-helper',
            '-coding',
            '-signaling',
            '-responsive',
            '-sensitive',
            '-catalytic',
            '-promoting',
            '-homologous',
            '-producing',
            '-regulate',
            '-presenting',
            '-processing',
            '-inducing'
            ]   
        for suffix in suffixes:
            x = len(suffix)
            if token[-x:] == suffix:
                token = token[:-x]
        return(token)
        
    def _tokenize_description(self, description):
        # Takes the description and tokenizes it, returning a list of tokens
        tokens = word_tokenize(description)
        # Remove punctuation
        tokens = [x.lower() for x in tokens if type(re.search(r'\w+', x)) == re.Match]
        tokens = [self._dehyphenate_token(token) for token in tokens]
        return(tokens)
    
    def _get_description_dict(self, path):
        # Takes a UniProt search file, and isolates protein descriptions from
        # it, returning it as a dictionary of accession / description key / values
        uniprot_file = open(path)
        
        uniprot_descriptions = {}
            
        previous_AC = False
        cc_on = False
        ac_split = []
        
        lines = uniprot_file.readlines()
        for line in lines:
            line = line.strip()
            if line[0:2] == "AC":
                if previous_AC:
                    ac_split = ac_split + line[5:-2].split('; ')
                 
                else:# its the first line of entries
                    previous_AC = True
                    ac_split = line[5:-2].split('; ')
                    function_string = ''
            
            if line[0:2] == "CC":
                if line[0:18] == 'CC   -!- FUNCTION:':
                    cc_on = True
                    function_string = function_string + line[19:]
                elif line[0:8] == 'CC   -!-':
                    if cc_on == True:
                        for ac in ac_split:
                            uniprot_descriptions[ac] = function_string
                    cc_on = False
                elif cc_on == True:
                    if line[-1] == '-':
                        function_string = function_string + line[8:-1]
                    else:
                        function_string = function_string + line[8:]
            
            else:
                previous_AC = False
                if cc_on == True:
                    cc_on = False
                    for ac in ac_split:
                        uniprot_descriptions[ac] = function_string
        return(uniprot_descriptions)
    
    def _generate_freq_dict(self, protein_list):
        # Generates a dictionary of word occurances from a list of proteins
        freq_dict = {}
        for protein in protein_list:
            if protein in self._uniprot_descriptions:
                description = self._uniprot_descriptions[protein]
                description = self._process_description(description)
                tokens = self._tokenize_description(description)
                for token in tokens:
                    if token in freq_dict:
                        freq_dict[token] = freq_dict[token] + 1
                    else:
                        freq_dict[token] = 1
        return(freq_dict)
    
    def _calculate_G2(self, a, b, c, d):
        # Calculates the G squared value for a word from the count of that word
        # in the subset (a), the count of that word in the control set (b), the
        # total word count of the subset (c), and the total word count of the
        # control set (d)
        E1 = c*(a+b)/(c+d)
        E2 = d*(a+b)/(c+d)
        G2 = 2*((a*np.log(a/E1)) + (b*np.log(b/E2)))
        return(G2)
    
    def generate_significance_table(self, subset, control):
        """
        Calculates the significance of the differences in the occurances of
        words in the UniProt functional descriptions for a selected subset
        of proteins and a control set - e.g. all proteins or all proteins that
        do not appear in the subset.

        Parameters
        ----------
        subset : 1-D array-like
            The UniProt accessions of proteins in the test set
        control : 1-D array-like
            The UniProt accessions of proteins in the control set

        Returns
        -------
        A pandas DataFrame with the p-values calculated for each of the
        tokenized words in the descriptions. p-Values can only be caulcualted
        when the word appears in the descriptions of both sets at least once
        

        """
        # Generates a pandas DataFrame with the p-values for each word, if it
        # exists in both datasets
        subset_dict = self._generate_freq_dict(subset)
        control_dict = self._generate_freq_dict(control)
        subset_freq_table = pd.DataFrame.from_dict(subset_dict, orient='index', columns=['Subset occurances'])
        control_freq_table = pd.DataFrame.from_dict(control_dict, orient='index', columns=['Control occurances'])
        freq_table = subset_freq_table.join(control_freq_table, how='outer')
        freq_table = freq_table.fillna(0)
        freq_table = freq_table.astype(np.int64)
        subset_total, control_total = np.sum(freq_table)
        freq_table = freq_table.assign(G2=lambda x: self._calculate_G2(x['Subset occurances'], x['Control occurances'], subset_total, control_total))
        freq_table['p'] = 1 - freq_table['G2'].apply(chi2.cdf, args=[1])
        return(freq_table)

    def __init__(self, path):
        """
        
        Initializes with a UniProt search result file so that the functional
        descriptions for all proteins can be extracted.

        Parameters
        ----------
        path : str
            Path to a UniProt search result text file

        Returns
        -------
        None.

        """
        self._uniprot_descriptions = self._get_description_dict(path)