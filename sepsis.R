library(readxl)
library(purrr)
library(tidyverse)

# Convert original data to col-based (columns are variables)
col_names <- read_excel("Data/sepsis.xlsx", 
                        sheet = 1, 
                        range = "A1:A171",
                        col_names = FALSE) %>% pull(1)

mb <- read_excel("Data/sepsis.xlsx",
                 sheet = 1,
                 range = "B1:IV171",
                 col_names = FALSE) %>%
    t() %>%
    as.data.frame() %>%
    type.convert(as.is = TRUE)

colnames(mb) <- col_names
rownames(mb) <- NULL

mb %>% count(Timepoint)
mb %>% count(SepsisControl)

# phyloseq Package --------------------------------------------------------
# https://joey711.github.io/phyloseq/import-data.html

library(phyloseq)

# OTU table
otu_mat <- mb %>% select(contains('D_0')) %>% as.matrix() %>% t()

# https://stackoverflow.com/questions/5006716/getting-the-text-that-follows-after-the-regex-match
taxonomy <- stringr::str_split(rownames(otu_mat), ";") %>%
    map(~ stringr::str_extract(., "(?<=__).*")) # will match everything after "__"

colnames(otu_mat) <- mb %>% pull(index)
rownames(otu_mat) <- paste0("OTU", 1:nrow(otu_mat))

# Taxonomy table
tax_mat <- taxonomy %>%
    tibble::enframe() %>%
    mutate(name = paste0("OTU", name)) %>%
    mutate(value = map(value, ~ setNames(., c("Domain", "Phylum", "Class", 
                                              "Order", "Family", "Genus", "Species")))
    ) %>% 
    unnest_wider(value) %>%
    column_to_rownames("name") %>%
    as.matrix()

# Turn into phyloseq objects
physeq <- phyloseq(otu_table(otu_mat, taxa_are_rows = TRUE), 
                   tax_table(tax_mat))
physeq

# Add sample data
sampledata <- mb %>%
    select(Timepoint, SepsisControl, index) %>%
    column_to_rownames("index") %>%
    sample_data()
physeq <- merge_phyloseq(physeq, sampledata)
sample_variables(physeq)

otu_table(physeq)[1:5, 1:5]
tax_table(physeq)[1:5, 1:7]

# Plot
plot_bar(physeq, x = "Family", y = "Abundance", facet_grid = "Timepoint~.") 
plot_bar(physeq, x = "Family", y = "Abundance", facet_grid = "Timepoint~SepsisControl") 
plot_bar(physeq %>% 
             filter_taxa(function(x) sum(x > 0) > (0.1*length(x)), prune = TRUE) %>%
             transform_sample_counts(function(x) x / sum(x)),  # divide by library size 
         x = "Family", y = "Abundance", facet_grid = "Timepoint~SepsisControl",
         fill = "Family") +
    geom_bar(aes(fill=Family), stat="identity", position="stack") + # get rid of the dividing lines
    scale_fill_viridis_d() 
    
plot_bar(physeq %>% 
             filter_taxa(function(x) sum(x > 0) > (0.2*length(x)), prune = TRUE),
         x = "Timepoint", y = "Abundance", fill = "Family") +
    facet_grid(vars(Family), vars(SepsisControl), scales = "free_y") +
    theme(strip.text.y = element_blank())

# Ordination
# physeq_sepsis <- prune_samples(get_variable(physeq, "SepsisControl") == "Sepsis", physeq)
# get_variable(physeq_sepsis, "SepsisControl") %>% unique()

# physeq_sepsis_ord <- ordinate(physeq_sepsis, method = "MDS", distance = "bray")
# physeq_sepsis %>% plot_ordination(physeq_sepsis_ord, 
#                                   type = "samples", 
#                                   color = "Timepoint") +
#     geom_point(size=3)

# MDS
dist_method <- "canberra"
#dist_method <- "bray"
physeq_ord <- ordinate(physeq, method = "MDS", distance = dist_method)
physeq %>% plot_ordination(physeq_ord, type = "samples", 
                           color = "SepsisControl",
                           title = "Sample Ordination") +
    scale_color_viridis_d(option = "plasma") +
    theme_bw()

# Shepard diagram
physeq_sh <- MASS::Shepard(distance(physeq, method = dist_method), physeq_ord$vectors)
plot(physeq_sh, pch = '.')
lines(physeq_sh$x, physeq_sh$yf, type = 'S', col = 'red')

physeq %>% plot_ordination(physeq_ord, type = "samples", 
                           color = "Timepoint", shape = "SepsisControl",
                           title = "Sample Ordination") +
    facet_wrap(~SepsisControl) +
    geom_point(size=3) +
    scale_color_viridis_d(option = "plasma") +
    theme_bw()

# CCA
physeq_ord <- ordinate(physeq, method = "CCA") # distance method not required
physeq %>% plot_ordination(physeq_ord, type = "split", color = "Timepoint", label = "Family") 
physeq %>% plot_ordination(physeq_ord, type = "taxa", label = "Family") +
    geom_text(aes(label = Family), size=3, color="blue")

physeq_by_family <- physeq %>% tax_glom("Family") 
tax_table(physeq_by_family)[, "Family"] %>% as.character() %>% stringr::str_subset("P.*")
# Rare OTU found in 3 samples, corresponding to the bottome left corner of the biplot
subset_taxa(physeq_by_family, Family == "Paenibacillaceae") %>% otu_table()



# vegan::cca(otu_table(physeq)@.Data %>% t()) -> ccaout
# plot(ccaout)


# ML with bartMachine -----------------------------------------------------

keep_otu <- function(x, freq_cut = 99/1){
    if(length(unique(x)) == 1)
        return(FALSE)
    # Ratio of freq of most common to the second most common value
    # as in recipe step_nzv
    s <- sort(table(x), decreasing = TRUE)
    return (isTRUE(s[1] / s[2] < freq_cut))
}

physeq_pruned <- physeq %>% 
    filter_taxa(keep_otu, prune = TRUE)


options(java.parameters = "-Xmx4000m") 
library(bartMachine)
set_bart_machine_num_cores(4)

X <- otu_table(physeq_pruned)@.Data %>% t() %>% as.data.frame()
y <- get_variable(physeq, "SepsisControl")
y <- as.factor(y)

# Fit with default params
bart_machine <- bartMachine(X, y)
bart_machine

# Check out-of-sample performance
oos_stats <- k_fold_cv(X, y, k_folds = 10)
oos_stats$confusion_matrix

# Variable importance
investigate_var_importance(bart_machine, num_replicates_for_avg = 20)
tax_table(physeq)[c("OTU61", "OTU156", "OTU140")]

# Variable selection
vs <- var_selection_by_permute(bart_machine, bottom_margin = 5, 
                               num_permute_samples = 100,
                               num_reps_for_avg = 20)
vs$important_vars_local_names
vs$important_vars_global_max_names
vs$important_vars_global_se_names

# Interaction effect
interaction_investigator(bart_machine, num_replicates_for_avg = 20, 
                         num_var_plot = 6) # no interaction detected

# Partial dependence
pd_plot(bart_machine, j = "OTU61")
pd_plot(bart_machine, j = "OTU156", levs = c(0.05, seq(0.1, 0.8, 0.05), 0.85))
pd_plot(bart_machine, j = "OTU140", levs = c(0.05, seq(0.1, 0.8, 0.05), 0.85))
