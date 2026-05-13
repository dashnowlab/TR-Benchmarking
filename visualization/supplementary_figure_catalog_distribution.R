#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 1000)

suppressPackageStartupMessages({
  library(zoo)
})

# =========================================================
# USER SETTINGS
# =========================================================
bed_file   <- "G:/Shared drives/DashnowLab/_Projects and Writing/TR-benchmarking-paper/data/benchmark_v2_output/benchmark-catalog-v2.longtr.bed"
cyto_file  <- "G:/Shared drives/DashnowLab/_Projects and Writing/TR-benchmarking-paper/data/benchmark_v2_output/cytoBandIdeo.txt"

out_dir    <- "G:/Shared drives/DashnowLab/_Projects and Writing/TR-benchmarking-paper/data/benchmark_v2_output"

# 4-row layout including chrX and chrY
contigs.row1 <- c(1:6)
contigs.row2 <- c(7:12)
contigs.row3 <- c(13:18)
contigs.row4 <- c(19:22, "X", "Y")

# appearance
density_fill <- "gray35"
density_line <- "black"
spacer_bp    <- 16000000

# plot settings
png_width  <- 3000
png_height <- 1260
png_res    <- 300

pdf_width  <- 15
pdf_height <- 6.4

# smoothing
smooth_k_100kb <- 7
smooth_k_1mb   <- 5

# optional fixed y limits
ymax_100kb <- NULL
ymax_1mb   <- NULL

# =========================================================
# HELPERS
# =========================================================
safe_rollmean <- function(x, k = 5) {
  if (length(x) == 0) return(x)
  if (length(x) == 1) return(x)
  
  good <- is.finite(x)
  if (!any(good)) return(rep(NA_real_, length(x)))
  
  x2 <- x
  x2[!good] <- median(x2[good], na.rm = TRUE)
  
  pad_n <- floor(k / 2)
  pad_left  <- rep(x2[1], pad_n)
  pad_right <- rep(x2[length(x2)], pad_n)
  xpad <- c(pad_left, x2, pad_right)
  
  sm <- zoo::rollmean(xpad, k = k, align = "center")
  sm <- sm[seq_len(length(x2))]
  sm <- pmax(sm, 0)
  
  sm
}

make_windows_one_chr <- function(chr, chr_len, bin_size) {
  starts <- seq(1, chr_len, by = bin_size)
  ends   <- pmin(starts + bin_size - 1, chr_len)
  data.frame(chr = chr, start = starts, end = ends)
}

interval_overlap_count <- function(win_start, win_end, feat_start, feat_end) {
  n <- length(win_start)
  if (length(feat_start) == 0) return(rep(0, n))
  
  first_idx <- findInterval(feat_start, win_end)
  first_idx[first_idx < 1] <- 1
  
  last_idx <- findInterval(feat_end, win_start)
  last_idx[last_idx < 1] <- 0
  
  valid <- which(first_idx <= last_idx & last_idx >= 1 & first_idx <= n)
  diff_vec <- rep(0, n + 1)
  
  if (length(valid) > 0) {
    fi <- first_idx[valid]
    li <- pmin(last_idx[valid], n)
    diff_vec[fi] <- diff_vec[fi] + 1
    diff_vec[li + 1] <- diff_vec[li + 1] - 1
  }
  
  cumsum(diff_vec)[seq_len(n)]
}

band_col <- function(stain) {
  if (stain == "gneg")    return("white")
  if (stain == "gpos25")  return("gray85")
  if (stain == "gpos50")  return("gray65")
  if (stain == "gpos75")  return("gray40")
  if (stain == "gpos100") return("black")
  if (stain == "gvar")    return("gray75")
  if (stain == "stalk")   return("gray90")
  if (stain == "acen")    return("#b2182b")
  return("white")
}

build_density_matrix <- function(bed, chr_lengths, bin_size, smooth_k) {
  chr_ids <- names(chr_lengths)
  out <- vector("list", length(chr_ids))
  names(out) <- chr_ids
  
  for (chr in chr_ids) {
    chr_len <- chr_lengths[chr]
    wins <- make_windows_one_chr(chr, chr_len, bin_size)
    
    x <- bed[bed$chr == chr, , drop = FALSE]
    counts <- interval_overlap_count(
      win_start  = wins$start,
      win_end    = wins$end,
      feat_start = x$start,
      feat_end   = x$end
    )
    
    dens_sm <- safe_rollmean(counts, k = smooth_k)
    
    out[[chr]] <- data.frame(
      chr = chr,
      start = wins$start,
      end = wins$end,
      density = dens_sm
    )
  }
  
  out
}

draw_cytoband_ideogram <- function(chr_name, xoff, y0, y1, chr_len, cyto) {
  chr_cyto <- cyto[cyto$chr == chr_name, , drop = FALSE]
  
  if (nrow(chr_cyto) == 0) {
    rect(xoff, y0, xoff + chr_len, y1, col = "white", border = "gray35", lwd = 0.5, xpd = NA)
    return(invisible(NULL))
  }
  
  for (j in seq_len(nrow(chr_cyto))) {
    rect(
      xleft = xoff + chr_cyto$start[j],
      xright = xoff + chr_cyto$end[j],
      ybottom = y0,
      ytop = y1,
      col = band_col(chr_cyto$gieStain[j]),
      border = NA,
      xpd = NA
    )
  }
  
  rect(
    xleft = xoff,
    xright = xoff + chr_len,
    ybottom = y0,
    ytop = y1,
    col = NA,
    border = "gray35",
    lwd = 0.5,
    xpd = NA
  )
}

draw_top_legend <- function(legend_text) {
  par(mar = c(0, 0, 0, 0), bg = "white")
  plot.new()
  
  x0 <- 0.31
  y0 <- 0.58
  box_w <- 0.022
  box_h <- 0.22
  
  rect(
    xleft = x0,
    ybottom = y0 - box_h / 2,
    xright = x0 + box_w,
    ytop = y0 + box_h / 2,
    col = adjustcolor(density_fill, alpha.f = 0.45),
    border = density_line,
    xpd = NA
  )
  
  text(
    x = x0 + box_w + 0.015,
    y = y0,
    labels = legend_text,
    adj = c(0, 0.5),
    cex = 0.95
  )
}

draw_one_row <- function(contigs, mat, chr_lengths, cyto, ymax,
                         fill_col = density_fill, line_col = density_line,
                         show_y = TRUE, spacer = spacer_bp) {
  
  chr_names <- ifelse(grepl("^chr", contigs), contigs, paste0("chr", contigs))
  row_lengths <- as.numeric(chr_lengths[chr_names]) + spacer
  total_len <- sum(row_lengths) - spacer
  
  par(
    mar = c(0.65, ifelse(show_y, 2.7, 0.45), 1.15, 0.25),
    bty = "n",
    bg = "white"
  )
  
  plot(
    x = c(0, total_len),
    y = c(-0.26 * ymax, 1.18 * ymax),
    type = "n",
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    xaxs = "i"
  )
  
  rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
       col = "white", border = NA)
  
  y.at <- pretty(c(0, ymax), n = 5)
  y.at <- y.at[y.at >= 0 & y.at <= ymax]
  
  abline(h = y.at, col = "gray82", lty = 3, lwd = 0.6)
  abline(h = 0, col = "black", lwd = 0.8)
  
  if (show_y) {
    axis(2, at = y.at, labels = round(y.at, 2), las = 2, tck = -0.012, cex.axis = 0.75)
  }
  
  starts_offset <- c(0, cumsum(row_lengths))[-length(chr_names) - 1]
  ends_offset <- starts_offset + as.numeric(chr_lengths[chr_names])
  
  for (i in seq_along(chr_names)) {
    chr_name <- chr_names[i]
    chr_len <- as.numeric(chr_lengths[chr_name])
    xoff <- starts_offset[i]
    chr_dat <- mat[[chr_name]]
    
    xvals <- ((chr_dat$start + chr_dat$end) / 2) + xoff
    yvals <- chr_dat$density
    
    ok <- is.finite(yvals) & xvals >= xoff & xvals <= (xoff + chr_len)
    
    if (sum(ok) > 1) {
      x_use <- xvals[ok]
      y_use <- pmax(yvals[ok], 0)
      
      if (length(y_use) >= 3) {
        y_use[length(y_use) - 1] <- 0.75 * y_use[length(y_use) - 1] + 0.25 * y_use[length(y_use)]
        y_use[length(y_use)]     <- 0.85 * y_use[length(y_use)]
      }
      
      if (x_use[1] > xoff) {
        x_use <- c(xoff, x_use)
        y_use <- c(y_use[1], y_use)
      }
      
      if (tail(x_use, 1) < (xoff + chr_len)) {
        x_use <- c(x_use, xoff + chr_len)
        y_use <- c(y_use, max(0, 0.7 * tail(y_use, 1)))
      }
      
      polygon(
        x = c(x_use, rev(x_use)),
        y = c(y_use, rep(0, length(y_use))),
        col = adjustcolor(fill_col, alpha.f = 0.45),
        border = NA
      )
      
      lines(x_use, y_use, col = line_col, lwd = 0.55, lend = "round", ljoin = "round")
    }
    
    segments(x0 = xoff, y0 = 0, x1 = xoff + chr_len, y1 = 0, col = "black", lwd = 0.7, lend = "round")
    
    text(
      x = mean(c(xoff, ends_offset[i])),
      y = 1.11 * ymax,
      labels = chr_name,
      cex = 0.60,
      xpd = NA
    )
    
    y0 <- -0.18 * ymax
    y1 <- -0.05 * ymax
    
    draw_cytoband_ideogram(
      chr_name = chr_name,
      xoff = xoff,
      y0 = y0,
      y1 = y1,
      chr_len = chr_len,
      cyto = cyto
    )
    
    rect(
      xleft = ends_offset[i],
      xright = xoff + row_lengths[i],
      ybottom = par("usr")[3],
      ytop = par("usr")[4],
      col = "white",
      border = NA
    )
  }
  
  box(col = "gray30", lwd = 0.6)
}

save_plot_set <- function(mat, chr_lengths, cyto, ymax, legend_text, png_file, pdf_file) {
  
  png(png_file, width = png_width, height = png_height, res = png_res, bg = "white")
  layout(matrix(c(1, 2, 3, 4, 5), nrow = 5, byrow = TRUE),
         heights = c(0.16, 1, 1, 1, 1))
  
  draw_top_legend(legend_text)
  draw_one_row(contigs.row1, mat, chr_lengths, cyto, ymax, show_y = TRUE)
  draw_one_row(contigs.row2, mat, chr_lengths, cyto, ymax, show_y = TRUE)
  draw_one_row(contigs.row3, mat, chr_lengths, cyto, ymax, show_y = TRUE)
  draw_one_row(contigs.row4, mat, chr_lengths, cyto, ymax, show_y = TRUE)
  
  dev.off()
  
  pdf(pdf_file, width = pdf_width, height = pdf_height, bg = "white")
  layout(matrix(c(1, 2, 3, 4, 5), nrow = 5, byrow = TRUE),
         heights = c(0.16, 1, 1, 1, 1))
  
  draw_top_legend(legend_text)
  draw_one_row(contigs.row1, mat, chr_lengths, cyto, ymax, show_y = TRUE)
  draw_one_row(contigs.row2, mat, chr_lengths, cyto, ymax, show_y = TRUE)
  draw_one_row(contigs.row3, mat, chr_lengths, cyto, ymax, show_y = TRUE)
  draw_one_row(contigs.row4, mat, chr_lengths, cyto, ymax, show_y = TRUE)
  
  dev.off()
}

# =========================================================
# READ BED
# =========================================================
bed <- read.table(bed_file, header = FALSE, sep = "\t", quote = "", comment.char = "")
if (ncol(bed) < 3) stop("BED file must have at least 3 columns")
bed <- bed[, 1:3]
colnames(bed) <- c("chr", "start", "end")

valid_chr <- c(paste0("chr", 1:22), "chrX", "chrY")
bed <- bed[bed$chr %in% valid_chr, , drop = FALSE]
bed$start <- as.numeric(bed$start) + 1
bed$end   <- as.numeric(bed$end)
bed <- bed[!is.na(bed$start) & !is.na(bed$end) & bed$start <= bed$end, , drop = FALSE]

# =========================================================
# READ CYTOBANDS
# =========================================================
cyto <- read.table(cyto_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
if (ncol(cyto) < 5) stop("cytoBandIdeo file must have 5 columns")
colnames(cyto) <- c("chr", "start", "end", "name", "gieStain")

cyto <- cyto[cyto$chr %in% valid_chr, , drop = FALSE]
cyto$start <- as.numeric(cyto$start) + 1
cyto$end   <- as.numeric(cyto$end)

chr_lengths <- tapply(cyto$end, cyto$chr, max)
chr_order <- c(paste0("chr", 1:22), "chrX", "chrY")
chr_lengths <- chr_lengths[chr_order]
chr_lengths <- chr_lengths[!is.na(chr_lengths)]

# =========================================================
# VERSION 1: 100 kb bins
# =========================================================
mat_100kb <- build_density_matrix(
  bed = bed,
  chr_lengths = chr_lengths,
  bin_size = 100000,
  smooth_k = smooth_k_100kb
)

all_y_100kb <- unlist(lapply(mat_100kb, function(x) x$density))
all_y_100kb <- all_y_100kb[is.finite(all_y_100kb)]
y_100kb <- if (is.null(ymax_100kb)) {
  if (length(all_y_100kb) == 0) 1 else max(pretty(c(0, max(all_y_100kb, na.rm = TRUE)), n = 4))
} else ymax_100kb

save_plot_set(
  mat = mat_100kb,
  chr_lengths = chr_lengths,
  cyto = cyto,
  ymax = y_100kb,
  legend_text = "TR loci per 100 kb bin, smoothed in rolling windows",
  png_file = file.path(out_dir, "TR_catalog_count_per_100kb_4rows.png"),
  pdf_file = file.path(out_dir, "TR_catalog_count_per_100kb_4rows.pdf")
)

# =========================================================
# VERSION 2: 1 Mb bins
# =========================================================
mat_1mb <- build_density_matrix(
  bed = bed,
  chr_lengths = chr_lengths,
  bin_size = 1000000,
  smooth_k = smooth_k_1mb
)

all_y_1mb <- unlist(lapply(mat_1mb, function(x) x$density))
all_y_1mb <- all_y_1mb[is.finite(all_y_1mb)]
y_1mb <- if (is.null(ymax_1mb)) {
  if (length(all_y_1mb) == 0) 1 else max(pretty(c(0, max(all_y_1mb, na.rm = TRUE)), n = 4))
} else ymax_1mb

save_plot_set(
  mat = mat_1mb,
  chr_lengths = chr_lengths,
  cyto = cyto,
  ymax = y_1mb,
  legend_text = "TR loci per 1 Mb bin, smoothed in rolling windows",
  png_file = file.path(out_dir, "TR_catalog_count_per_1Mb_4rows.png"),
  pdf_file = file.path(out_dir, "TR_catalog_count_per_1Mb_4rows.pdf")
)

message("Done.")
message("Created:")
message(file.path(out_dir, "TR_catalog_count_per_100kb_4rows.png"))
message(file.path(out_dir, "TR_catalog_count_per_100kb_4rows.pdf"))
message(file.path(out_dir, "TR_catalog_count_per_1Mb_4rows.png"))
message(file.path(out_dir, "TR_catalog_count_per_1Mb_4rows.pdf"))