# main.R
library(yaml)
library(optparse)

option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config/default_config.yaml",
              help = "Path to configuration file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

config <- read_yaml(opt$config)
message(sprintf("Starting pipeline for project: %s", config$project_name))
