<div align="center">
<img src="logo.png" alt="QPromoter logo" width="15%" />
</div>

# QPromoters
Your go-to tool for predicting the strength of Saccharomyces cerevisiae promoters! 

# Website
[QPromoters.com](https://www.qpromoters.com/) (designed by @Mukulikaa)

# About
We present a theoretical model to describe the relationship between promoter strength and nucleotide sequence in *Saccharomyces cerevisiae*. We infer from our analysis that the -49 to 10 sequence with respect to the Transcription Start Site represents the minimal region that can be used to predict the promoter strength. This tool takes advantage of this fact to quickly quantify the strength of the promoters. See [our publication](https://doi.org/10.1101/2021.04.27.441621)  for more information.

# Usage
Download this repository and run the program ```terminal_app.py```  
Replace the file ```EPDPromoters/all_promters_49_10.fa``` to add or remove promoters use to generate [PSSM](https://en.wikipedia.org/wiki/Position_weight_matrix). Sequences for named promoters used in the program are also retrived from this file.

# Dependencies
1. Biopython (Installation instructions: https://biopython.org/wiki/Download)
2. Numpy (Installation instructions: https://numpy.org/install/)
3. Matplotlib (Installation instructions: https://matplotlib.org/stable/users/installing.html)

# Cite
QPromoters: Sequence based prediction of promoter strength in Saccharomyces cerevisiae
Devang Haresh Liya, Mirudula Elanchezhian, Mukulika Pahari, Nithishwer Mouroug Anand, Shivani Suresh, Nivedha Balaji, Ashwin Kumar Jainarayanan
bioRxiv 2021.04.27.441621; doi: https://doi.org/10.1101/2021.04.27.441621 
