# Additional functions for snp index calculation and snp-wise fisher exact test.
# No additional library requirement, simplely source this file if you need.
# panyq 2022-05-30


# Calculate snp index from AD string
# Input: "10,20"
# Output: 2/3
#
# Only first two alleles will be considered
get_snp_index <- function(AD) {
    single_snp_index <- function(x) {
        x <- as.numeric(strsplit(x, ",")[[1]])
        return(x[2] / (x[1] + x[2]))
    }
    sapply(AD, single_snp_index, USE.NAMES=F)
}

# Do fisher exact test to ADs.
# Input: "1,2", "3,4", "5,6" ...
# Output: a p- value, return NA if encounter error.
#
# Inputs can also be columns like above, this function have been vectorized with mapply
AD_fisher_test <- function(..., mc.cores = 2, workspace = 2e9, simulate.p.value = T) {
    single_fisher_test <- function(..., workspace, simulate.p.value) {
        ADs <- list(...)
        ADs <- as.data.frame(lapply(ADs, function(x) as.numeric(strsplit(x, ",")[[1]][1:2])))
        p_value <- NA
        tryCatch(
            {p_value <- fisher.test(ADs, workspace = workspace, simulate.p.value = simulate.p.value)$p.value},
            warning = function(war){}, error = function(err){}
        )
        return(p_value)
    }
    parallel::mcmapply(single_fisher_test, ..., workspace = workspace, simulate.p.value = simulate.p.value,
                       USE.NAMES=F, mc.cores = mc.cores)
}

# Sliding window calculate average snp index
# Input: a data.frame, containing snp_index columns
# Output: a data.frame, containing per window snp_index
#
# Note: this function is not efficient at all, and it use parallel::mcmapply to speed up
# calculation process, provide more threads using mc.cores argument.
get_sliding_window_snp_index <- function(df, chr_name = "CHROM", pos_name = "POS",
                                         snp_index_name = "snp_index", window = 5e6,
                                         shift = 1e4, mc.cores = 2) {

    snp_index_column_names <- names(df)[grep(snp_index_name, names(df))]

    get_sliding_window_template <- function(df) {
        window_chr <- unique(df[[chr_name]])
        window_start <- seq(min(df[["POS"]]), max(df[["POS"]]) - shift, by = shift)
        window_end <- seq(min(df[["POS"]]), max(df[["POS"]]) - shift, by = shift) + shift
        data.frame(window_chr, window_start, window_end)
    }

    sliding_window_template <- do.call(
        rbind,
        lapply(
            split(df, df[[chr_name]]),
            get_sliding_window_template
        )
    )

    get_window_average <- function(window_chr, window_start, window_end, column_name) {
        window_df <- df[df[[chr_name]] == window_chr & window_start <= df[[pos_name]] & df[[pos_name]] < window_end,]
        if (dim(window_df)[1] == 0) {
            return(NA)
        } else {
            return(mean(window_df[[column_name]], na.rm = T))
        }
    }

    for (column_name in snp_index_column_names) {
        sliding_window_template[[column_name]] <- parallel::mcmapply(
            get_window_average,
            sliding_window_template[["window_chr"]], 
            sliding_window_template[["window_start"]], 
            sliding_window_template[["window_end"]], 
            column_name,
            mc.cores = mc.cores
        )
    }
    return(sliding_window_template)
}

