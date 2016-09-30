library(dplyr); library(ggplot2); library(gridExtra); library(vegan); library(magrittr)

# this function combines a list of final output tables form sag-mg-recruit
combine_smr_tables <- function(tablelist){
    all<-data.frame()
    for (t in tablelist){
        df = read.table(t, sep="\t", header=TRUE)
        all <- rbind(all, df) 
        }
    return(all)
    }

# returns a table with a row per metagenome that summarises metagenome overall representation
# reads_per_mbp: number of metagenomic reads recruited per mpb
# prop_mgreads_per_mbp: proportion total meteganome reads recruited per total_sag_mbp
# pct_metagenome: percent metagenomic reads recruited to all sags** this value will be the bar plot 
# underneath the y-axis of the heatmap

summarise_mgs <- function(df){
    mg_summary <- df %>% filter(sag == "concatenated_sags") %>% .[,c(1:8)]
    total_sag_bp <- df %>% .[,c(1,10)] %>% unique %>% filter(sag != "concatenated_sags") %>% .[,2] %>% sum 
    mg_summary$total_sag_bp <- total_sag_bp
    mg_summary$sag_size_mbp <- total_sag_bp/1000000
    mg_summary <- mg_summary %>% mutate(reads_per_mbp=total_reads_recruited/sag_size_mbp, 
                                        prop_mgreads_per_mbp=(total_reads_recruited/sag_size_mbp)/mg_read_count,
                                       pct_metagenome=(total_reads_recruited/mg_read_count)*100)
    underplot <- mg_summary[,c(2,ncol(mg_summary))]
    return(underplot)
}



# summarise the number of total mg reads recruited to each SAG
# for the side plot of the heatmap

summarise_sags <- function(pairs){
    total_mg_reads <- pairs[,8] %>% unique %>% sum
    sideplot <- pairs %>% group_by(sag) %>% summarise(pct_reads_recruited_allmgs = (sum(total_reads_recruited)/total_mg_reads)*100)
    return(sideplot)
}


# process an input file defining the order in which sags or mgs appear 
# in the final heatmap
import_order <- function(file_path){
    order <- read.table(file_path, header=TRUE)
    order$levels <- seq(1, nrow(order))
    return(order)
}

# Plot fuctions 
main_heatmap <- function(pairs, 
                         lowcolor="white",
                         midcolor="lightblue", 
                         highcolor="darkblue",
                        sagorder="None",
                        mgorder="None"){
    if (mgorder != "None"){
        pairs$metagenome <- with(pairs, factor(metagenome, levels=mgorder[,1]))
    }
    if (sagorder != "None"){
        pairs$sag <- with(pairs, factor(sag, levels=sagorder[,1]))
        }
    p1 <- pairs %>% ggplot(aes(metagenome, sag)) + geom_tile(aes(fill=prop_mgreads_per_mbp)) + theme_bw() 
    pmain <- p1 + theme(axis.ticks = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
    pmain <- pmain + labs(x=" ", fill="") + scale_fill_gradientn(colours=c(lowcolor,midcolor,highcolor), values=c(0,0.00001,1))
    pmain <- pmain  + theme(legend.key.height=unit(50, "pt"), legend.position="none") 
    pmain <- pmain + theme(panel.border = element_blank(), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), 
                           panel.background = element_blank()) 
    return(pmain)
}

# plot legend for heatmap separately from the heatmap itself: 
hm_legend <- function(pairs, lowcolor="white",midcolor="lightblue", highcolor="darkblue"){
    
    legend <- pairs %>% ggplot(aes(metagenome, sag)) + geom_blank(aes(fill=prop_mgreads_per_mbp))
    legend <- legend + theme(axis.text = element_blank(),
                    axis.title = element_blank(),
                    line = element_blank(),
                    panel.background = element_blank()) + labs(fill="")
    legend <- legend + theme(legend.key.width=unit(80, "pt"), 
                             legend.position=c(.5, .5), 
                     legend.direction="horizontal") 
    legend <- legend + scale_fill_gradientn(colours=c(lowcolor,midcolor,highcolor), values=c(0,0.00001,1))
    return(legend)

}

# function to create the bar plot that sits underneath the heatmap
mg_plot <- function(df, mgorder="None"){
    underplot <- summarise_mgs(df)
    if (mgorder != "None"){
        underplot$metagenome <- with(underplot, factor(metagenome, levels=mgorder[,1]))
    }
    u <- underplot %>% ggplot(aes(metagenome, pct_metagenome)) + geom_bar(stat="identity", fill="white", color="black") +  scale_y_reverse()
    u <- u + theme(axis.line.y = element_line(color="black"))+ scale_x_discrete() 
    punder <- u + labs(y="", x="") + theme(axis.text.x  = element_text(angle=90, face="bold"))
    punder <- punder + theme(panel.border = element_blank(), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), 
                             panel.background = element_blank(), 
                             axis.ticks = element_blank(),
                             axis.line.y = element_line(color="black"),
                             axis.text.x = element_blank()) # remove this if you want axis labels for the y-axis
    return(punder)
    }

# function to create the bar plot that sits to the left of the heatmap
sag_plot <- function(df, sagorder="None"){
    sideplot <- summarise_sags(df)
    if (sagorder != "None"){
        sideplot$sag <- with(sideplot, factor(sag, levels=sagorder[,1]))
    }
    s <- sideplot %>% ggplot(aes(sag, pct_reads_recruited_allmgs)) + geom_bar(stat="identity", fill="white", color="black") +  theme_classic()+ coord_flip()
    pside <- s + theme(axis.ticks = element_blank(), axis.text.y = element_blank(), axis.title = element_blank())
    pside <- pside + theme(axis.text.x = element_text(angle = -90))
    pside <- pside + theme(panel.border = element_blank(), 
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(), 
                           panel.background = element_blank(),
                           axis.line.x = element_line(color="black")) 
    return(pside)

}

# empty plot required for creating the plot array
empty_plot <- function(df){
    pempty <-ggplot(df, aes(sag, metagenome)) + geom_blank() + theme(axis.text = element_blank(),
                axis.title = element_blank(),
                line = element_blank(),
                panel.background = element_blank())
    return(pempty)
}

# creates an adjacency table from an edgelist; needed for clustering
create_adjtbl <- function(edgelist){
    yaxis = unique(edgelist[,2])
    xaxis = unique(edgelist[,1])
    together <- data.frame(xaxis)
    colnames(together)[1] <- colnames(edgelist)[1]
    for (m in seq(1,length(yaxis))){
        es <- edgelist %>% filter(.[,2]==yaxis[m])
        colnames(es)[3] <- as.character(yaxis[m])
        to_combine <- data.frame(es[,3])
        colnames(to_combine)[1] <- colnames(es)[3]
        #together = merge(together, to_combine, by=colnames(together)[1])
        together = cbind(together, to_combine)
    }
    return(together)
    }

# sets the order based on heirarchical clustering
clust_order_x <- function(edgelist){
    adj <- create_adjtbl(edgelist)
    rownames(adj) <- adj[,1]
    adj[,1] <- NULL
    mat1 <- as.matrix(adj)
    xdist <- mat1 %>% vegdist(method="bray") 
    xdist[is.na(xdist)] <- 1
    xdist %>% hclust -> xclust
    xdend <- as.dendrogram(xclust)
    xorder1 <- order.dendrogram(xdend)
    data <- mat1[xorder1,]
    levels <- seq(1, nrow(data))
    xorder <- cbind.data.frame(rownames(data), levels)
    colnames(xorder)[1] <- colnames(edgelist)[1]
    return(xorder)
    }

plot_array <- function(df,
                       sag_ax_aln = 0,
                       hm_aln = 28,
                       mg_ax_aln = 105,
                       lowcolor="white",
                       midcolor="lightblue", 
                       highcolor="darkblue",
                       hclust_sags = FALSE,
                       hclust_mgs = FALSE,
                       sag_order = "None",
                       mg_order = "None"){
    mgplot <- mg_plot(df)
    pairs <- df %>% filter(sag != "concatenated_sags")
    mgorder="None"
    sagorder="None"
    if (hclust_sags==TRUE){
        sagorder = clust_order_x(pairs[,c(1,2,13)])
        sagplot <- sag_plot(pairs, sagorder=sagorder)
        pairs$sag <- with(pairs, factor(sag, levels=sagorder[,1]))
    } else if (sag_order != "None"){
        sagorder = import_order(sag_order)
        sagplot <- sag_plot(pairs, sagorder = sagorder)
    } else {
        sagplot <- sag_plot(pairs)
    }
    
    if (hclust_mgs==TRUE){
        mgorder = clust_order_x(pairs[,c(2,1,13)])
        mgplot <- mg_plot(df, mgorder=mgorder)
        pairs$metagenome <- with(pairs, factor(metagenome, levels=mgorder[,1]))
    } else if (mg_order != "None"){
        mg_order = import_order(mg_order)
        mgplot <- mg_plot(df, mgorder = mg_order)
    } else {
        mgplot <- mg_plot(df)
    }
    # pairs %>% head %>% print
    heatmap <- main_heatmap(pairs, lowcolor=lowcolor, midcolor=midcolor, highcolor=highcolor, mgorder=mgorder, sagorder=sagorder)
    empty <- empty_plot(df)
    legend <- hm_legend(pairs, lowcolor=lowcolor, midcolor=midcolor, highcolor=highcolor)
    heatmap <- heatmap + theme(plot.margin=unit(c(0,0,hm_aln,0), "pt"))
    mgplot <- mgplot + theme(plot.margin=unit(c(-hm_aln,0,0,mg_ax_aln),"pt"))
    sagplot <- sagplot + theme(plot.margin=unit(c(0,0,sag_ax_aln,0), "pt")) #+ theme(plot.margin=unit(c(0,0,-25,0), "pt"))
    empty <- empty + theme(plot.margin=unit(c(50,0,0,0), "pt"))
    legend <- legend + theme(plot.margin=unit(c(0,0,0,mg_ax_aln), "pt"))
    
    final <- grid.arrange(heatmap, sagplot, mgplot, empty, legend, empty,
             ncol = 2, nrow = 3, widths = c(6, 1), heights = c(8, 2, 1))
    return(final)
    
    }