################################
# Taxonomic assignment table ####
##################################
# Using tax table only, grouped at each taxrank

# By taxrank, number of taxa found by tool, shown independently by approach
# and within approach, intersect size
# One table per dataset
# Is there a clear break across taxranks ?

## Alternatively : series of Venn diagrams ?

####################
# Alpha diversity ###
######################
# Sample_id, Dataset, idx, idx_value, CCE, Approach
# Import from previous work or re-generate?

# Tests at lowest (mostly) coherent taxonomic rank, than one step lower
# Or straight up at species level ? because that's one reason to do shotgun

### 1. TECHNICAL COMPARISON
# One plot per index (Richness and Shannon, to begin with)
# facet dataset ~ approach; x = tool and y = idx_value
# dots + grouped lines
# Compute mean change? 

### 2. HYPOTHESIS COMPARISON
# Test between groups
# Simple boxplots with pvalues

###################
# Beta diversity ###
#####################

# using collapsed pairwise matrices (BC and rAitchison):
# Sample_pair, Dataset, idx, idx_value, CCE, Approach

### 1. TECHNICAL COMPARISON
# Essentially the same plot as alphadiv

### 2. HYPOTHESIS COMPARISONS
# 2.1. PCoA comparison
# Procruste comparison of pcoas ?

# 2.2. perMANOVA 
# Presented as in poster? 



