
wordcloud_heatmap <- function(mnn_list){
  
  
  source("https://gist.githubusercontent.com/jokergoo/bfb115200df256eeacb7af302d4e508e/raw/14f315c7418f3458d932ad749850fd515dec413b/word_cloud_grob.R")
  
  ht <- Heatmap(mat, 
                name = "mat",
                show_row_names = F,
                show_column_names = F,
                top_annotation = column_ha,
                left_annotation = row_ha,
                row_km = mnn_list$metadata$cluster,
                column_km = mnn_list$metadata$cluster)
  
  input_dat <- cbind(mnn_list$metadata, mnn_list$reduction_coords)
  tfidf.celltype <- tf_idf(input_dat, "celltype", terms_per_cluster = 2)
  keywords <- lapply(unique(tfidf.celltype$cluster), function(x){
    dat <- data.frame(tfidf.celltype[tfidf.celltype$cluster==x,])
    row.names(dat) <- dat$word
    return(dat)
  }) %>% `names<-`(unique(tfidf.celltype$cluster))
  
  align_to = split(seq_len(nrow(mat)), mnn_list$metadata$cluster)
  # align_to = align_to[names(align_to) != "0"]
  align_to = align_to[names(align_to) %in% names(keywords)]
  # align_to
  
  
  fontsize_range = c(4, 16)
  gbl = lapply(names(align_to), function(nm) {
    kw = keywords[[nm]][, 1]
    freq = keywords[[nm]][, 2]
    fontsize = keywords[[nm]]$tf_idf
    # scale_fontsize(freq, rg = c(1, max(10, freq)), fs = fontsize_range)
    word_cloud_grob(text = kw, fontsize = fontsize)
  })
  names(gbl) = names(align_to)
  # gbl
  
  margin = unit(8, "pt")
  gbl_h = lapply(gbl, function(x) convertHeight(grobHeight(x), "cm") + margin)
  gbl_h = do.call(unit.c, gbl_h)
  
  gbl_w = lapply(gbl, function(x) convertWidth(grobWidth(x), "cm"))
  gbl_w = do.call(unit.c, gbl_w)
  gbl_w = max(gbl_w) + margin
  
  panel_fun = function(index, nm) {
    # background
    grid.rect(gp = gpar(fill = "#DDDDDD", col = NA))
    # border
    grid.lines(c(0, 1, 1, 0), c(0, 0, 1, 1), gp = gpar(col = "#AAAAAA"), 
               default.units = "npc")
    gb = gbl[[nm]]
    # a viewport within the margins
    pushViewport(viewport(x = margin/2, y = margin/2, 
                          width = grobWidth(gb), height = grobHeight(gb),
                          just = c("left", "bottom")))
    grid.draw(gb)
    popViewport()
  }
  
  ht = ht + rowAnnotation(keywords = anno_link(align_to = align_to, 
                                               which = "row", panel_fun = panel_fun, 
                                               size = gbl_h, gap = unit(2, "mm"), 
                                               width = gbl_w + unit(5, "mm"), # 5mm for the link
                                               link_gp = gpar(fill = "#DDDDDD", col = "#AAAAAA"), 
                                               internal_line = FALSE)) # you can set it to TRUE to see what happens
  draw(ht, ht_gap = unit(2, "pt"))
}


