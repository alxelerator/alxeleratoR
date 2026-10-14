#' Plot a Sankey diagram of taxonomic composition from a phyloseq object
#'
#' Produces a Sankey diagram flowing from a broad taxonomic rank down to a
#' finer rank, using the \pkg{sankey} package (Csárdi & Weiner). Each rank
#' column is scaled to 100% so that all axes have equal visual height,
#' allowing proportional comparison across ranks. Flow width is proportional
#' to abundance aggregated across all samples with the chosen summary function.
#'
#' @section Aggregation:
#' Each ASV is first summarised across all samples with \code{fun}. Identical
#' taxonomy paths are then summed to obtain one abundance value per unique path
#' through the selected rank range. Flows represent the chosen statistic
#' (sum / mean / median) of read counts per taxonomy group across samples.
#'
#' @section Column scaling (100 percent):
#' Every rank column is independently scaled so that the combined height of
#' all node bars equals a constant reference value. This normalises away
#' rank-to-rank conservation differences and makes the proportional breakdown
#' at each level immediately comparable by eye.
#'
#' @section Top-N filtering and Other bundling:
#' When \code{top > 0}, only the top-N taxa at \code{top_rank} (by total
#' aggregated abundance) are shown as real taxonomy. All excluded taxa are
#' bundled into a single "Other" group at every rank - meaning a non-top path
#' becomes Other-Phylum > Other-Class > ... > Other-Genus. This keeps
#' the number of nodes at every axis manageable. By default, \code{top_rank}
#' is the finest rank (\code{to_rank}) so that the top-N filter is applied at
#' genus (or equivalent) level; the broader ranks then display only the real
#' taxonomy traced back from those top genera.
#'
#' @section Ordering:
#' With \code{sort.by = "abundance"} (default), nodes at every rank are
#' sorted by descending total abundance within their parent group, with the
#' most abundant taxon at the bottom. Set \code{sort.by = "alpha"} for
#' alphabetical order.
#'
#' @section Output and return value:
#' The plot is rendered as a side effect to the current graphics device or
#' saved to a file if \code{file} is set. The format is inferred from the
#' file extension (\code{.svg} or \code{.png}). The \code{sankey} object is
#' returned invisibly and can be re-plotted with \code{sankey::sankey()}.
#'
#' @param ps phyloseq object
#' @param from_rank highest (broadest) taxonomic rank to start from.
#'   Default: first rank in \code{rank_names(ps)}
#' @param to_rank lowest (finest) taxonomic rank to end at.
#'   Default: last rank in \code{rank_names(ps)}
#' @param fun aggregation function applied per ASV across all samples. One of
#'   \code{"sum"}, \code{"mean"}, or \code{"median"}. Default \code{"mean"}
#' @param top limit display to top-N taxa at \code{top_rank} by total
#'   aggregated abundance. \code{0} = show all. Default \code{30}
#' @param top_rank rank at which the top-N filter is applied. Must lie within
#'   the selected range. Default \code{NULL} resolves to \code{to_rank} (the
#'   finest rank), which keeps the number of nodes at all axes manageable.
#' @param sort.by how to order taxa at each axis. \code{"abundance"} (default)
#'   places the largest taxon at the bottom within each parent group.
#'   \code{"alpha"} sorts alphabetically.
#' @param show.other when \code{top > 0}, bundle excluded taxa into a single
#'   "Other" group shown at the top of each axis. Default \code{TRUE}
#' @param other.label prefix for the bundled non-top group label; the rank
#'   name is appended automatically. Default \code{"Other"}
#' @param other.color fill colour for "Other" nodes and flows.
#'   Default \code{"lightgrey"}
#' @param na.label prefix for missing taxonomy entries; rank name is appended.
#'   Default \code{"Unknown"}
#' @param label show taxon name labels on nodes. Default \code{TRUE}
#' @param label.size character expansion factor (\code{cex}) for node labels.
#'   Default \code{0.75}
#' @param label.min.frac suppress labels for nodes whose total abundance is
#'   below this fraction of the overall total. Default \code{0.01} (1 percent)
#' @param legend show a legend for the \code{from_rank} fill colours.
#'   Default \code{TRUE}
#' @param legend.size character expansion factor for the legend. Default \code{0.7}
#' @param title plot title. \code{NULL} generates an automatic title.
#' @param alpha transparency applied to node and flow colours (0-1).
#'   Default \code{0.6}
#' @param file path to save the plot. \code{NULL} = render to current device.
#'   Format inferred from extension: \code{.svg} or \code{.png}.
#' @param width plot width in inches when saving. Default \code{14}
#' @param height plot height in inches when saving. Default \code{9}
#'
#' @returns The \code{sankey} object, invisibly.
#'
#' @examples
#' \dontrun{
#' # Default: Phylum to Genus, mean, top 30 genera
#' ps_plot_sankey(ps, from_rank = "Phylum", to_rank = "Genus")
#'
#' # Top 8 phyla (filter at Phylum level), save as SVG
#' ps_plot_sankey(ps, from_rank = "Phylum", to_rank = "Family",
#'                top = 8, top_rank = "Phylum",
#'                file = "sankey_phylum.svg")
#'
#' # Top 15 families (default filter at to_rank), Phylum to Genus
#' ps_plot_sankey(ps, from_rank = "Phylum", to_rank = "Genus",
#'                top = 15, top_rank = "Family")
#'
#' # Alphabetical sort, no Other bundling
#' ps_plot_sankey(ps, from_rank = "Phylum", to_rank = "Genus",
#'                top = 0, sort.by = "alpha", show.other = FALSE)
#' }
#'
#' @author Alex Bossers \email{a.bossers@@uu.nl}
#'
#' @import phyloseq
#' @import dplyr
#' @import sankey
#' @importFrom stats median
#' @importFrom grDevices svg png dev.off adjustcolor
#' @importFrom graphics title legend
#' @importFrom tools file_ext
#' @export
#'
#' @note
#'   version: 2.1
#'   versiondate: 2026-07-01

ps_plot_sankey <- function(
  ps,
  from_rank      = NULL,
  to_rank        = NULL,
  fun            = "mean",
  top            = 30,
  top_rank       = NULL,
  sort.by        = "abundance",
  show.other     = TRUE,
  other.label    = "Other",
  other.color    = "lightgrey",
  na.label       = "Unknown",
  label          = TRUE,
  label.size     = 0.75,
  label.min.frac = 0.01,
  legend         = TRUE,
  legend.size    = 0.7,
  title          = NULL,
  alpha          = 0.6,
  file           = NULL,
  width          = 14,
  height         = 9
) {

  # -- validate and select ranks -----------------------------------------
  all_ranks <- rank_names(ps)
  if (is.null(from_rank)) from_rank <- all_ranks[1]
  if (is.null(to_rank))   to_rank   <- all_ranks[length(all_ranks)]

  if (!from_rank %in% all_ranks)
    stop(paste0("from_rank '", from_rank, "' not in rank_names(ps): ",
                paste(all_ranks, collapse = ", ")))
  if (!to_rank %in% all_ranks)
    stop(paste0("to_rank '", to_rank, "' not in rank_names(ps): ",
                paste(all_ranks, collapse = ", ")))

  idx_from <- which(all_ranks == from_rank)
  idx_to   <- which(all_ranks == to_rank)
  if (idx_from >= idx_to)
    stop("from_rank must be a higher (broader) rank than to_rank.")

  selected_ranks <- all_ranks[idx_from:idx_to]

  # Default: filter at the finest rank so intermediate bars stay manageable
  if (is.null(top_rank)) top_rank <- to_rank
  if (!top_rank %in% selected_ranks)
    stop(paste0("top_rank '", top_rank, "' must be one of: ",
                paste(selected_ranks, collapse = ", ")))

  # -- extract and clean taxonomy ----------------------------------------
  tax_df <- as.data.frame(tax_table(ps))[, selected_ranks, drop = FALSE]
  for (r in selected_ranks) {
    vals <- as.character(tax_df[[r]])
    vals[is.na(vals) | vals == "" | vals == "NA"] <- paste0(na.label, " [", r, "]")
    tax_df[[r]] <- vals
  }

  # -- per-ASV abundance summary across all samples ----------------------
  otu_mat <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) {
    abund <- switch(fun,
      sum    = rowSums(otu_mat),
      mean   = rowMeans(otu_mat),
      median = apply(otu_mat, 1, median),
      stop(paste0("fun must be 'sum', 'mean', or 'median'. Got: '", fun, "'"))
    )
  } else {
    abund <- switch(fun,
      sum    = colSums(otu_mat),
      mean   = colMeans(otu_mat),
      median = apply(otu_mat, 2, median),
      stop(paste0("fun must be 'sum', 'mean', or 'median'. Got: '", fun, "'"))
    )
  }
  tax_df$Abundance <- abund[rownames(tax_df)]

  # -- aggregate: one row per unique taxonomy path -----------------------
  full_alluv_df <- tax_df %>%
    group_by(across(all_of(selected_ranks))) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    filter(Abundance > 0)

  # -- top-N filter at top_rank ------------------------------------------
  top_vals <- NULL
  if (top > 0) {
    top_vals <- full_alluv_df %>%
      group_by(across(all_of(top_rank))) %>%
      summarise(.total = sum(Abundance), .groups = "drop") %>%
      arrange(desc(.total)) %>%
      slice_head(n = top) %>%
      pull(all_of(top_rank))
    alluv_df <- full_alluv_df %>% filter(!!sym(top_rank) %in% top_vals)
  } else {
    alluv_df <- full_alluv_df
  }

  # -- Other bundling: collapse ALL ranks for non-top paths --------------
  # This is important: preserving real taxonomy above top_rank would expose
  # rare phyla as separate colors and fragment the legend.
  if (top > 0 && show.other) {
    non_top_df <- full_alluv_df %>% filter(!(!!sym(top_rank) %in% top_vals))
    if (nrow(non_top_df) > 0) {
      for (r in selected_ranks)
        non_top_df[[r]] <- paste0(other.label, " [", r, "]")
      other_agg <- non_top_df %>%
        group_by(across(all_of(selected_ranks))) %>%
        summarise(Abundance = sum(Abundance), .groups = "drop")
      alluv_df <- bind_rows(alluv_df, other_agg)
    } else {
      show.other <- FALSE
    }
  }

  # -- sort: abundance (recursive within parent) or alphabetical ---------
  other_pat <- paste0("^", other.label, " \\[")
  real_rows  <- alluv_df %>% filter(!grepl(other_pat, !!sym(from_rank)))
  other_rows <- alluv_df %>% filter( grepl(other_pat, !!sym(from_rank)))

  if (sort.by == "abundance") {
    from_totals <- real_rows %>%
      group_by(across(all_of(from_rank))) %>%
      summarise(.total = sum(Abundance), .groups = "drop") %>%
      arrange(.total)
    real_rows[[from_rank]] <- factor(real_rows[[from_rank]],
                                     levels = from_totals[[from_rank]])

    for (i in seq_along(selected_ranks)[-1]) {
      rank_i       <- selected_ranks[i]
      parent_ranks <- selected_ranks[seq_len(i - 1)]
      rank_totals  <- real_rows %>%
        group_by(across(all_of(c(parent_ranks, rank_i)))) %>%
        summarise(.total = sum(Abundance), .groups = "drop")
      for (pr in parent_ranks)
        if (is.factor(real_rows[[pr]]))
          rank_totals[[pr]] <- factor(rank_totals[[pr]], levels = levels(real_rows[[pr]]))
      rank_totals <- rank_totals %>% arrange(across(all_of(parent_ranks)), .total)
      real_rows[[rank_i]] <- factor(real_rows[[rank_i]],
                                    levels = unique(as.character(rank_totals[[rank_i]])))
    }
    real_rows <- real_rows %>% arrange(across(all_of(selected_ranks)))
  } else {
    real_rows <- real_rows %>% arrange(across(all_of(selected_ranks)))
  }

  alluv_df <- bind_rows(real_rows, other_rows)

  # Per-rank factor levels: Other always first (= top of column), real taxa below
  for (r in selected_ranks) {
    all_lv   <- unique(as.character(alluv_df[[r]]))
    other_lv <- paste0(other.label, " [", r, "]")
    real_lv  <- setdiff(all_lv, other_lv)
    alluv_df[[r]] <- factor(alluv_df[[r]],
                            levels = c(if (other_lv %in% all_lv) other_lv, real_lv))
  }

  # -- colour palette for from_rank groups --------------------------------
  other_fill_lv  <- paste0(other.label, " [", from_rank, "]")
  fill_levels    <- levels(alluv_df[[from_rank]])
  regular_levels <- setdiff(fill_levels, other_fill_lv)
  col_pal        <- setNames(colorDistinct(n = length(regular_levels)), regular_levels)
  if (other_fill_lv %in% fill_levels) col_pal[other_fill_lv] <- other.color

  # -- build node data frame (plain data.frame required by simplegraph) --
  # Node ID = "taxon [rank]" for cross-rank uniqueness.
  # Use factor levels (Other first, then ascending abundance) so y="simple"
  # stacks them correctly: last in data frame = bottom with gravity="bottom".
  node_rows <- lapply(selected_ranks, function(r) {
    vals <- if (is.factor(alluv_df[[r]])) levels(alluv_df[[r]]) else
              unique(as.character(alluv_df[[r]]))
    data.frame(
      id   = paste0(vals, " [", r, "]"),
      rank = r,
      name = vals,
      stringsAsFactors = FALSE
    )
  })
  nodes <- do.call(rbind, node_rows)
  nodes <- nodes[!duplicated(nodes$id), ]

  # Node x position (rank index) used for 100% normalization below
  nodes$x_pos <- match(nodes$rank, selected_ranks)

  # Node colour: inherit the from_rank ancestor colour
  fill_lookup <- lapply(selected_ranks, function(r) {
    alluv_df %>%
      select(all_of(c(from_rank, r))) %>%
      distinct() %>%
      mutate(node_id = paste0(as.character(!!sym(r)), " [", r, "]"),
             fill    = as.character(!!sym(from_rank))) %>%
      select(node_id, fill)
  })
  fill_lookup_df <- do.call(rbind, lapply(fill_lookup, as.data.frame))
  fill_lookup_df <- fill_lookup_df[!duplicated(fill_lookup_df$node_id), ]

  nodes$fill_group <- fill_lookup_df$fill[match(nodes$id, fill_lookup_df$node_id)]
  nodes$col <- ifelse(
    is.na(nodes$fill_group) | !(nodes$fill_group %in% names(col_pal)),
    other.color,
    col_pal[nodes$fill_group]
  )

  # -- node abundance totals and label suppression -----------------------
  total_flow <- sum(alluv_df$Abundance)
  node_abunds <- sapply(seq_len(nrow(nodes)), function(j) {
    r  <- nodes$rank[j]
    nm <- nodes$name[j]
    sum(alluv_df$Abundance[as.character(alluv_df[[r]]) == nm])
  })

  nodes$label <- ifelse(!label | node_abunds < label.min.frac * total_flow,
                        "", nodes$name)
  nodes$cex   <- label.size

  # -- 100% column normalisation -----------------------------------------
  # Scale node sizes within each rank column to a common reference so that
  # all axes have equal visual height.
  col_totals   <- tapply(node_abunds, nodes$x_pos, sum)
  nodes$size   <- node_abunds / col_totals[as.character(nodes$x_pos)]

  # -- build edge data frame (plain data.frame required by simplegraph) --
  chr_df <- alluv_df %>% mutate(across(all_of(selected_ranks), as.character))

  edge_list <- lapply(seq_along(selected_ranks)[-length(selected_ranks)], function(i) {
    r1 <- selected_ranks[i]
    r2 <- selected_ranks[i + 1]
    chr_df %>%
      group_by(across(all_of(c(r1, r2)))) %>%
      summarise(.w = sum(Abundance), .groups = "drop") %>%
      filter(.w > 0) %>%
      transmute(
        from       = paste0(.data[[r1]], " [", r1, "]"),
        to         = paste0(.data[[r2]], " [", r2, "]"),
        weight     = .w / total_flow,   # normalise to match 100% column scale
        colorstyle = "gradient"
      )
  })
  # bind_rows returns tibble; simplegraph needs plain data.frame
  # (tibble[,1] is a tibble, not a vector, so is.character() == FALSE)
  edges <- as.data.frame(bind_rows(edge_list), stringsAsFactors = FALSE)
  edges$from <- as.character(edges$from)
  edges$to   <- as.character(edges$to)

  # Apply alpha to node colours
  nodes$col <- adjustcolor(nodes$col, alpha.f = alpha)

  # -- build sankey object -----------------------------------------------
  # Pass only id, col, size, label, cex. Do NOT pass x_pos (numeric second
  # column fails simplegraph's character-id validation on nodes).
  sk_nodes <- as.data.frame(
    nodes[, c("id", "col", "size", "label", "cex")],
    stringsAsFactors = FALSE
  )
  # y = "simple" stacks nodes in nodes data frame order (our abundance sort).
  # Much faster than "optimal" for large graphs; "optimal" is O(n^2).
  sk <- make_sankey(
    nodes   = sk_nodes,
    edges   = edges,
    y       = "simple",
    gravity = "bottom"
  )

  # -- auto title --------------------------------------------------------
  fun_label <- switch(fun, sum = "Sum", mean = "Mean", median = "Median")
  if (is.null(title)) {
    top_info <- if (top > 0) paste0(", top ", top, " ", top_rank) else ""
    title <- paste0("Taxonomy Sankey  -  ", from_rank, " to ", to_rank,
                    "  (", fun_label, " across samples", top_info, ")")
  }

  # -- render ------------------------------------------------------------
  do_plot <- function() {
    sankey(sk)
    graphics::title(main = title, cex.main = 0.9, line = 0.5)
    if (legend) {
      legend_cols <- adjustcolor(col_pal, alpha.f = alpha)
      graphics::legend(
        "topright",
        legend = names(col_pal),
        fill   = legend_cols,
        title  = from_rank,
        cex    = legend.size,
        bty    = "n"
      )
    }
  }

  if (!is.null(file)) {
    ext <- tolower(tools::file_ext(file))
    if (ext == "svg") {
      grDevices::svg(file, width = width, height = height)
    } else {
      grDevices::png(file, width = width, height = height,
                     units = "in", res = 150)
    }
    tryCatch(do_plot(), finally = grDevices::dev.off())
    message("Saved: ", file)
  } else {
    do_plot()
  }

  invisible(sk)
}
