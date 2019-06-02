# Specific_Projects
Project-Specific files for automating data analysis for paired patch-clamp recordings from electrically-coupled neurons. Processing, descriptive statistics, ANOVA, and plotting are all coded in Python (3.73). Dunnett's posthoc test is performed in R. 

## For the 6-2-19 China trip:
### Paired Recording data:
1. empty_paired_data.csv will contain cleaned raw data
2. Run paired_data.csv through py_calc to calculate the %cc etc...
3. Run the calculated files ('_attempted.csv', '_connected.csv' etc through py_scatter.py to generate desc stats, One-Way ANOVA, Dunnett's Posthoc Test, and plots for each measurement


