# Motivation
while preparing for the discussion with Mike Springer last night, I asked if it will a general pattern that the number of Pho4 consensus motifs is smaller in lineages with Pho2-independent Pho4. The logic is that in those lineages, Pho4 alone is able to induce binding and therefore there should be stronger purifying selection against spurious sites not meant to be regulated.

# Previous analysis
In the analysis I did last September, I used RSAT (originally at http://rsat.ulb.ac.be/rsat/, now at http://fungi.rsat.eu/) DNA-pattern tool (genome-wide version) to count the number of occurances of CACGTG motifs in S. cerevisiae and C. glabrata genomes, specifcally in the 800bp upstream of any protein-coding genes and also within the CDS. This led to the result then (22 sep 2014).

However, when I tried to repeat the analysis, with the intention of extending it to more species, I found the RSAT website has changed a little. What’s more, I realize that I didn’t record all the details / parameters I used during my last analysis, and thus couldn’t perfectly repeat the analysis. 

# Goal
To answer the question (underlined) above, I will first use a coarse approach by counting the number of Pho4 consensus motifs in different species. 

# Notes
The approach taken here relies on an important assumption, namely Pho4 binding preference hasn’t changed. Let me remember this issue, but put it aside for the moment. Another practical question is where should I count? In the whole genome? In genes? In their promoters? How to define their promoters?

After talking to Elmar and some research on my own, I decide to focus on the 800bp upstream region, while ignoring the CDS for the moment. Elmar mentioned that in his ChIP data, he rarely see those sites occupied. Besides, I notice that some studies use 1000 bp upstream, and don’t account for cases where this region overlaps adjacent ORFs. In my work here, I use 800 bp and only extend to the boundary of adjacent ORF.
