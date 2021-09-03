### Functions 
subset.haplohh <- function (x, select.hap = NULL, select.mrk = NULL, min_perc_geno.hap = NA, 
    min_perc_geno.mrk = 100, min_maf = NA, max_alleles = NA, 
    verbose = TRUE, ...) 
{
    if (!is.na(min_perc_geno.hap) & (min_perc_geno.hap < 0 | 
        min_perc_geno.hap > 100)) {
        stop("min_perc_geno.hap should lie in the interval [0,100].", 
            call. = FALSE)
    }
    if (is.na(min_perc_geno.mrk) | min_perc_geno.mrk <= 0 | min_perc_geno.mrk > 
        100) {
        stop("min_perc_geno.mrk should lie in the interval (0,100].", 
            call. = FALSE)
    }
    if (!is.na(min_maf)) {
        if (!is.numeric(min_maf) | min_maf < 0 | min_maf > 0.5) {
            stop("min_maf should lie in the interval [0,0.5].", 
                call. = FALSE)
        }
    }
    if (!is.na(max_alleles)) {
        if (!is.numeric(max_alleles) | max_alleles < 2) {
            stop("max_alleles should be at least 2.", call. = FALSE)
        }
    }
    if (!is.haplohh(x)) {
        stop("Data is not a valid object of class haplohh.", 
            call. = FALSE)
    }
    if (min(dim(haplo(x))) == 0) {
        stop("Empty data object.", call. = FALSE)
    }
    if (!is.null(select.hap)) {
        if (verbose) 
            cat("* Subset haplotypes *\n")
        haplo <- x@haplo[select.hap, , drop = FALSE]
        if (nrow(haplo) == nhap(x)) {
            if (verbose) 
                cat("No haplotype discarded.\n")
        }
        else {
            if (verbose) 
                cat(nhap(x) - nrow(haplo), "haplotypes discarded.\n")
            x@haplo <- haplo
            if (verbose) 
                cat(nhap(x), "haplotypes remaining.\n")
        }
        if (nhap(x) == 0) {
            warning("No haplotype left after filtering on 'select.hap'.\n", 
                call. = FALSE, immediate. = TRUE)
            return(x)
        }
    }
    if (!is.null(select.mrk)) {
        if (verbose) 
            cat("* Subset markers *\n")
        haplo <- x@haplo[, select.mrk, drop = FALSE]
        if (ncol(haplo) == nmrk(x)) {
            if (verbose) 
                cat("No marker discarded.\n")
        }
        else {
            if (verbose) 
                cat(nmrk(x) - ncol(haplo), "markers discarded.\n")
            positions <- positions(x)
            names(positions) <- mrk.names(x)
            x@haplo <- haplo
            x@positions <- unname(positions[select.mrk])
            if (verbose) 
                cat(nmrk(x), "markers remaining.\n")
            if (nmrk(x) == 0) {
                warning("No marker left after filtering on 'select.mrk'.\n", 
                  call. = FALSE, immediate. = TRUE)
                return(x)
            }
        }
    }
    if (verbose) 
        cat("* Filtering data *\n")
    if (!is.na(min_perc_geno.hap)) {
        if (verbose) 
            cat("Discard haplotypes with less than", min_perc_geno.hap, 
                "% of genotyped markers.\n")
        hap_sel <- (100 * rowMeans(!is.na(x@haplo))) >= min_perc_geno.hap
        if (sum(hap_sel) == nhap(x)) {
            if (verbose) 
                cat("No haplotype discarded.\n")
        }
        else {
            if (verbose) 
                cat(nhap(x) - sum(hap_sel), "haplotypes discarded.\n")
            x@haplo <- x@haplo[hap_sel, , drop = FALSE]
            if (verbose) 
                cat(nhap(x), "haplotypes remaining.\n")
            if (nhap(x) == 0) {
                warning("No haplotype left after filtering of missing data.\n", 
                  "If applicable, reduce min_perc_geno.hap to allow for more missing data.\n", 
                  call. = FALSE, immediate. = TRUE)
                return(x)
            }
        }
    }
    if (verbose) 
        cat("Discard markers genotyped on less than", min_perc_geno.mrk, 
            "% of haplotypes.\n")
    mrk_sel <- (100 * colMeans(!is.na(x@haplo))) >= min_perc_geno.mrk
    if (sum(mrk_sel) == nmrk(x)) {
        if (verbose) 
            cat("No marker discarded.\n")
    }
    else {
        if (verbose) 
            cat(nmrk(x) - sum(mrk_sel), "markers discarded.\n")
        x@haplo <- x@haplo[, mrk_sel, drop = FALSE]
        x@positions <- x@positions[mrk_sel]
        if (verbose) 
            cat(nmrk(x), "markers remaining.\n")
        if (nmrk(x) == 0) {
            warning("No marker left after filtering of missing data.\n", 
                "If applicable, reduce min_perc_geno.mrk to allow for more missing data.\n", 
                call. = FALSE, immediate. = TRUE)
            return(x)
        }
    }
    if (!is.na(min_maf)) {
        if (verbose) 
            cat("Discard markers with Minor Allele Frequency equal to or below", 
                min_maf, ".\n")
        mrk_sel <- apply(x@haplo, 2, function(x) {
            t <- tabulate(x + 1L)
            if (length(t) == 1) {
                return(FALSE)
            }
            else {
                alleles <- order(t, decreasing = TRUE)[1:2]
                return(min(t[alleles])/sum(t) > min_maf)
            }
        })
        if (sum(mrk_sel) == nmrk(x)) {
            if (verbose) 
                cat("No marker discarded.\n")
        }
        else {
            if (verbose) 
                cat(nmrk(x) - sum(mrk_sel), "markers discarded.\n")
            x@haplo <- x@haplo[, mrk_sel, drop = FALSE]
            x@positions <- x@positions[mrk_sel]
            if (verbose) 
                cat(nmrk(x), "markers remaining.\n")
            if (nmrk(x) == 0) {
                warning("No marker left after filtering on Minor Allele Frequency.\n", 
                  "If applicable, reduce min_maf to allow for less frequent minor alleles.", 
                  call. = FALSE, immediate. = TRUE)
                return(x)
            }
        }
    }
    if (!is.na(max_alleles)) {
        if (verbose) 
            cat("Discard markers with more than", max_alleles, 
                "different alleles.\n")
        mrk_sel <- apply(x@haplo, 2, function(x) {
            length(unique(na.omit(x))) <= max_alleles
        })
        if (sum(mrk_sel) == nmrk(x)) {
            if (verbose) 
                cat("No marker discarded.\n")
        }
        else {
            if (verbose) 
                cat(nmrk(x) - sum(mrk_sel), "markers discarded.\n")
            x@haplo <- x@haplo[, mrk_sel, drop = FALSE]
            x@positions <- x@positions[mrk_sel]
            if (verbose) 
                cat(nmrk(x), "markers remaining.\n")
            if (nmrk(x) == 0) {
                warning("No marker left after filtering on maximal number of different alleles.\n", 
                  call. = FALSE, immediate. = TRUE)
                return(x)
            }
        }
    }
    if (verbose) 
        cat("Data consists of", nhap(x), "haplotypes and", nmrk(x), 
            "markers.\n")
    multicity <- tabulate(apply(x@haplo, 2, function(x) {
        sum(tabulate(x + 1) != 0)
    }))
    names(multicity) <- seq_along(multicity)
    if (verbose) {
        cat("Number of mono-, bi-, multi-allelic markers:\n")
        cat(names(multicity), "\n")
        cat(multicity, "\n")
    }
    return(x)
}
is.haplohh <- function (x) 
{
    check <- (is(x, "haplohh") & validObject(x))
    if (check) {
        check <- is.integer(x@haplo)
        if (!check) {
            warning("Matrix 'haplo' of class 'haplohh' must be of type 'integer'.")
        }
    }
    if (check) {
        check <- nmrk(x) == length(positions(x))
        if (!check) {
            warning("Number of positions must be equal to number of markers.")
        }
    }
    return(check)
}
