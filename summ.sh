Rscript -e '
library(data.table)
odir="bigdata/2025-08-01/";
tt=fread(paste0(odir,"tt_3rep_3x_per_group.csv.gz"))


CpG = colSums(tt[, grep(".CpG",names(tt),value=T), with=F],  na.rm = TRUE)
uCpG = colSums(tt[, grep(".uCpG",names(tt),value=T), with=F],  na.rm = TRUE)
group = factor(substr(names(CpG),1,1))

df <- data.frame(CpG=CpG, uCpG=uCpG, group=group)
# Binomial GLM with E as reference
glm_fit <- glm(cbind(CpG, uCpG) ~ group, data=df, family=binomial)
summary(glm_fit)

# Pairwise contrasts (including Y vs W)
emm <- emmeans(glm_fit, ~ group)
contrast(emm, method = "pairwise", adjust = "none")  # or use adjust="tukey" if multiple comparisons correction is desired

# Mixed-effects example: if you had a higher-level nesting (e.g., batch or sample grouping),
# here is a placeholder with a random intercept per group (not very useful with only group levels):
glmer_fit <- glmer(cbind(CpG, uCpG) ~ group + (1|group), data=df, family=binomial)
summary(glmer_fit)

'

