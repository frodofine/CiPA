vroom::vroom("data/drug_block.csv") -> drug_block
source("newIC50_mcmc.R")
Hill_fitting("astemizole", drug_block, desiredChannels = "INa") -> astemizole_INa

plot_sensitivity(astemizole_INa$INa, "Astemizole INa")
plot_MCMC(astemizole_INa$INa, "Astemizole INa")
