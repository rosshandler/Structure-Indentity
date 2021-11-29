## Load packages
library(miloR)
library(scran)
library(scater)
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)

library("scProportionTest")

## Load plotting colours

leiden_colours <- c(
  "0"="#cb459b",
  "1"="#5db645",
  "2"="#c267d3",
  "3"="#a6b547",
  "4"="#6d64d5",
  "5"="#d39c3f",
  "6"="#6e83cc",
  "7"="#d25337",
  "8"="#4aadd6",
  "9"="#dc446c",
  "10"="#5bbe7d",
  "11"="#9761aa",
  "12"="#3a7f4b",
  "13"="#e088aa",
  "14"="#50bba7",
  "15"="#a24a5f",
  "16"="#747a32",
  "17"="#b77345",
  "18"="#3c4dbb",
  "19"="#9e4d70",
  "20"="#bbb633",
  "21"="#df8172",
  "22"="#677328",
  "23"="#a6542e",
  "24"="#91702c"
)

sample_colours <- c(
  "B6Kit"="#787ec7",
  "B7Kit"="#64ac48",
  "B8Kit"="#7062d0",
  "B6Fast"="#9a963f",
  "B7Fast"="#be5dc2",
  "B8Fast"="#4aaa86",
  "B6diss"="#d13d7e",
  "B7diss"="#4facd8",
  "B8diss"="#d04f33",
  "B6Unembed"="#be6c9b",
  "B8Unembed_A"="#c88742",
  "B8Unembed_B"="#c2605f"
)

condition_colours <- c(
  "Kit"="#e12f2f",
  "Fast"="#6889ff",
  "Diss"="#9a9100",
  "Unembed"="#00733a"
)

population_colours <- c(
"Ventricular RG" = "#02d7ec",
"Cycling RG" = "#4aadd6",
"Mitotic RG" = "#005dd8",
"Glycolytic RG" = "#0086cb",
"Transcriptionally active RG" = "#4aac8d",
"Differentiating RG" = "#6083d1",
"IPC"   = "#ff8ba6",
"Committed neurons" = "#90b600",
"Migrating excitatory neurons" = "#a6023e",
"DL enriched neurons"   = "#9e5d56",
"UL enriched neurons"   = "#fa2274",
"Mature excitatory neurons"="#d3b000",
"Inhibitory neurons"="#228B22",
"Cajal Retzius cells" = "#00ff00",
"Cortical hem" = "#9e4d70",
"Choroid plexus"  = "#935de6",
"High metabolism/protein translation"   = "#ff6500"
)

population_colours_manuscript <- c(
  "Mitotic RG" = "#005dd8",
  "Cycling RG" = "#4aadd6",
  "Differentiating RG" = "#6083d1",
  "Apical RG" = "#02d7ec",
  "Glycolytic RG" = "#0086cb",
  "High trancription RG" = "#4aac8d",
  "Chp" = "#935de6",
  "Cortical hem" = "#9e4d70",
  "Cycling IPCs" = "#d74897",
  "IPCs"   = "#ff8ba6",
  #"High biosynthesis"   = "#ff6500",
  "IN"="#228B22",
  "Mature IN" = "#90b600",
  "CR cells" = "#00ff00",
  "Migrating excitatory neurons" = "#a6023e",
  "UL neurons"   = "#fa2274",
  "DL neurons"   = "#9e5d56",
  "Mature excitatory neurons"="#d3b000"
)

## Modifying milo plotting colours
library(igraph)
library(ggraph)

#' @importFrom igraph is_igraph
.valid_graph <- function(x){
    # check for a valid graph
    if(isTRUE(is_igraph(x))){
        TRUE
    } else{
        FALSE
    }
}
.valid_graph <- function(x){
    # check for a valid graph
    if(isTRUE(is_igraph(x))){
        TRUE
    } else{
        FALSE
    }
}
plotNhoodGraph_adapted <- function (x, layout = "UMAP", colour_by = NA, subset.nhoods = NULL, 
    size_range = c(0.5, 3), node_stroke = 0.3, ...) 
{
    if (!.valid_graph(nhoodGraph(x))) {
        stop("Not a valid Milo object - neighbourhood graph is missing. Please run buildNhoodGraph() first.")
    }
    if (is.character(layout)) {
        if (!layout %in% names(reducedDims(x))) {
            stop(layout, "isn't in readucedDim(x) - choose a different layout")
        }
    }
    nh_graph <- nhoodGraph(x)
    if (!is.null(subset.nhoods)) {
        nh_graph <- igraph::induced_subgraph(nh_graph, vids = which(as.numeric(V(nh_graph)$name) %in% 
            unlist(nhoodIndex(x)[subset.nhoods])))
    }
    nh_graph <- permute(nh_graph, order(vertex_attr(nh_graph)$size, 
        decreasing = TRUE))
    if (is.character(layout)) {
        redDim <- layout
        layout <- reducedDim(x, redDim)[as.numeric(vertex_attr(nh_graph)$name), 
            ]
        if (!any(class(layout) %in% c("matrix"))) {
            warning("Coercing layout to matrix format")
            layout <- as(layout, "matrix")
        }
    }
    if (!is.na(colour_by)) {
        if (colour_by %in% colnames(colData(x))) {
            col_vals <- colData(x)[as.numeric(vertex_attr(nh_graph)$name), 
                colour_by]
            if (!is.numeric(col_vals)) {
                col_vals <- as.character(col_vals)
            }
            V(nh_graph)$colour_by <- col_vals
        }
        else {
            stop(colour_by, "is not a column in colData(x)")
        }
    }
    else {
        V(nh_graph)$colour_by <- V(nh_graph)$size
        colour_by <- "Nhood size"
    }
    if (colour_by %in% c("logFC")) {
        plot.g <- simplify(nh_graph)
        pl <- ggraph(simplify(nh_graph), layout = layout) + geom_edge_link0(aes(width = weight), 
            edge_colour = "grey66", edge_alpha = 0.2) + geom_node_point(aes(fill = colour_by, 
            size = size), shape = 21, stroke = node_stroke) + 
            scale_size(range = size_range, name = "Nhood size") + 
            scale_edge_width(range = c(0.2, 3), name = "overlap size") + 
            theme_classic(base_size = 14) + theme(axis.line = element_blank(), 
            axis.text = element_blank(), axis.ticks = element_blank(), 
            axis.title = element_blank())
    }
    else {
        pl <- ggraph(simplify(nh_graph), layout = layout) + geom_edge_link0(aes(width = weight), 
            edge_colour = "grey66", edge_alpha = 0.2) + geom_node_point(aes(fill = colour_by, 
            size = size), shape = 21, stroke = node_stroke) + 
            scale_size(range = size_range, name = "Nhood size") + 
            scale_edge_width(range = c(0.2, 3), name = "overlap size") + 
            theme_classic(base_size = 14) + theme(axis.line = element_blank(), 
            axis.text = element_blank(), axis.ticks = element_blank(), 
            axis.title = element_blank())
    }
    if (is.numeric(V(nh_graph)$colour_by)) {
        #pl <- pl + scale_fill_gradientn(colours = wesanderson::wes_palette("Zissou1", 11, type = "continuous"))
        #pl <- pl + scale_fill_gradientn(colours = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(unique(V(nh_graph)$colour_by))))
        #mycolors <- wesanderson::wes_palette("Zissou1", 11, type = "continuous")
        #mycolors[6] <- "white"
        #pl <- pl + scale_fill_gradientn(colours = mycolors)
        pl <- pl + scale_fill_gradientn(colours = c("red","red","white","blue","blue"))
    }
    else {
        mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(V(nh_graph)$colour_by)))
        pl <- pl + scale_fill_manual(values = mycolors, name = colour_by, 
            na.value = "white")
    }
    pl
}

plotNhoodGraphDA_adapted <- function (x, milo_res, alpha = 0.05, res_column = "logFC", ...) 
{
    if (!.valid_graph(nhoodGraph(x))) {
        stop("Not a valid Milo object - neighbourhood graph is missing. Please run buildNhoodGraph() first.")
    }
    if (is.character(layout)) {
        if (!layout %in% names(reducedDims(x))) {
            stop(layout, "is not in readucedDim(x) - choose a different layout")
        }
    }
    signif_res <- milo_res
    signif_res[signif_res$SpatialFDR > alpha, res_column] <- 0
    colData(x)[res_column] <- NA
    colData(x)[unlist(nhoodIndex(x)[signif_res$Nhood]), res_column] <- signif_res[, 
        res_column]
    plotNhoodGraph_adapted(x, colour_by = res_column, ...)
}


