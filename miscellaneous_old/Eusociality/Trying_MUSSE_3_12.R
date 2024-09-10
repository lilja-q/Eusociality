


library(ape)

library(phytools)
library(stringr)
library(readr)

library(devtools)

#install_github("uyedaj/treeplyr")

library(treeplyr)
library(ggtree)



###############################################
write_MuHiSSE_nex=function (x, file, format = "dna", datablock = TRUE, interleaved = TRUE, 
                            charsperline = NULL, gap = NULL, missing = NULL, nchars=NULL) 
{
  
  
  if(is.null(nchars)){
    symbols= " symbols=\"0123456789\";\n"
  }else if(nchars==2){
    symbols= " symbols=\"1234\";\n"
    
  }else if(nchars==3){
    symbols= " symbols=\"12345678\";\n"
    
    
  }
  
  format <- match.arg(toupper(format), c("DNA", "PROTEIN", 
                                         "STANDARD", "CONTINUOUS"))
  if (inherits(x, "DNAbin") && format != "DNA") {
    format <- "DNA"
    warning("object 'x' is of class DNAbin: format forced to DNA")
  }
  if (inherits(x, "AAbin") && format != "PROTEIN") {
    format <- "PROTEIN"
    warning("object 'x' is of class AAbin: format forced to PROTEIN")
  }
  indent <- "  "
  maxtax <- 5
  defcharsperline <- 80
  defgap <- "-"
  defmissing <- "?"
  if (is.matrix(x)) {
    if (inherits(x, "DNAbin")) 
      x <- as.list(x)
    else {
      xbak <- x
      x <- vector("list", nrow(xbak))
      for (i in seq_along(x)) x[[i]] <- xbak[i, ]
      names(x) <- rownames(xbak)
      rm(xbak)
    }
  }
  ntax <- length(x)
  nchars <- length(x[[1]])
  zz <- file(file, "w")
  if (is.null(names(x))) 
    names(x) <- as.character(1:ntax)
  fcat <- function(..., file = zz) cat(..., file = file, sep = "", 
                                       append = TRUE)
  find.max.length <- function(x) max(nchar(x))
  print.matrix <- function(x, dindent = "    ", collapse = "") {
    Names <- names(x)
    printlength <- find.max.length(Names) + 2
    if (!interleaved) {
      for (i in seq_along(x)) {
        sequence <- paste(x[[i]], collapse = collapse)
        taxon <- Names[i]
        thestring <- sprintf("%-*s%s%s", printlength, 
                             taxon, dindent, sequence)
        fcat(indent, indent, thestring, "\n")
      }
    }
    else {
      ntimes <- ceiling(nchars/charsperline)
      start <- 1
      end <- charsperline
      for (j in seq_len(ntimes)) {
        for (i in seq_along(x)) {
          sequence <- paste(x[[i]][start:end], collapse = collapse)
          taxon <- Names[i]
          thestring <- sprintf("%-*s%s%s", printlength, 
                               taxon, dindent, sequence)
          fcat(indent, indent, thestring, "\n")
        }
        if (j < ntimes) 
          fcat("\n")
        start <- start + charsperline
        end <- end + charsperline
        if (end > nchars) 
          end <- nchars
      }
    }
  }
  if (inherits(x, "DNAbin") || inherits(x, "AAbin")) 
    x <- as.character(x)
  fcat("#NEXUS\n[Data written by write.nexus.data.R, ", date(), 
       "]\n")
  NCHAR <- paste("NCHAR=", nchars, sep = "")
  NTAX <- paste0("NTAX=", ntax)
  DATATYPE <- paste0("DATATYPE=", format)
  if (is.null(charsperline)) {
    if (nchars <= defcharsperline) {
      charsperline <- nchars
      interleaved <- FALSE
    }
    else charsperline <- defcharsperline
  }
  if (is.null(missing)) 
    missing <- defmissing
  MISSING <- paste0("MISSING=", missing)
  if (is.null(gap)) 
    gap <- defgap
  GAP <- paste0("GAP=", gap)
  INTERLEAVE <- if (interleaved) 
    "INTERLEAVE=YES"
  else "INTERLEAVE=NO"
  if (datablock) {
    fcat("BEGIN DATA;\n")
    fcat(indent, "DIMENSIONS ", NTAX, " ", NCHAR, ";\n")
    if (format != "STANDARD") {
      fcat(indent, "FORMAT", " ", DATATYPE, " ", MISSING, 
           " ", GAP, " ", INTERLEAVE, ";\n")
    }
    else {
      fcat(indent, "FORMAT", " ", DATATYPE, " ", MISSING, 
           " ", GAP, " ", INTERLEAVE, symbols)
    }
    fcat(indent, "MATRIX\n")
    if (format != "CONTINUOUS") {
      print.matrix(x)
    }
    else {
      print.matrix(x, collapse = "\t")
    }
    fcat(indent, ";\nEND;\n\n")
  }
  else {
    fcat("BEGIN TAXA;\n")
    fcat(indent, "DIMENSIONS", " ", NTAX, ";\n")
    fcat(indent, "TAXLABELS\n")
    fcat(indent, indent)
    j <- 0
    for (i in seq_len(ntax)) {
      fcat(names(x[i]), " ")
      j <- j + 1
      if (j == maxtax) {
        fcat("\n", indent, indent)
        j <- 0
      }
    }
    fcat("\n", indent, ";\n")
    fcat("END;\n\nBEGIN CHARACTERS;\n")
    fcat(indent, "DIMENSIONS", " ", NCHAR, ";\n")
    if (format != "STANDARD") {
      fcat(indent, "FORMAT", " ", MISSING, " ", GAP, " ", 
           DATATYPE, " ", INTERLEAVE, ";\n")
    }
    else {
      fcat(indent, "FORMAT", " ", MISSING, " ", GAP, " ", 
           DATATYPE, " ", INTERLEAVE, " symbols=\"0123456789\";\n")
    }
    fcat(indent, "MATRIX\n")
    if (format != "CONTINUOUS") {
      print.matrix(x)
    }
    else {
      print.matrix(x, collapse = "\t")
    }
    fcat(indent, ";\nEND;\n\n")
  }
  close(zz)
}

################################################



soc_cleaned_da_silva <- read_csv("soc_names_cleaned_da_silva.csv") # load Da Silva Data
colnames(soc_cleaned_da_silva)


cbind(soc_cleaned_da_silva$...1, soc_cleaned_da_silva$`Full Name`, soc_cleaned_da_silva$Sociality)

tree=ape::read.tree("BEE_mat7_fulltree.nwk") # Load phylogeny data

tree=ape::read.tree("BEE_mat7_fulltree_tplo35_sf20lp.nwk") # Load phylogeny data

is.ultrametric(tree)

pdf("test_tree.pdf", width=20, height=20)

  plot(tree, cex=0.5)

dev.off()



string_split <- function(x) { # Define a function that will collapse names to 'Genus_species'
  split_text = str_split_fixed(x, pattern = "_", 4)[1:2]
  paste0(split_text, collapse = "_")
}


new_labels <- unlist(lapply(tree$tip.label, string_split)) # Apply this function

tree$tip.label <- new_labels # Re-assign new labels in format 'Genus_species'
soc_cleaned_da_silva$Sociality <- gsub("\\?","/ambiguous", soc_cleaned_da_silva$Sociality) # change "?" to "/ambiguous" for easier processing
species <- as.vector(soc_cleaned_da_silva[["Full Name"]]) # get species list from Da Silva




#Prune tree to species in Da Silva data
intersects <- intersect(species, tree$tip.label)
tips_to_drop <- as.character(union(setdiff(intersects, tree$tip.label), setdiff(tree$tip.label, intersects)))
pruned.tree<-drop.tip(tree,tips_to_drop) # Prune the tree



length(pruned.tree$node.label)



pruned.tree = force.ultrametric(pruned.tree)

#Also make a separate, subsetted Da Silva table containing just the relevant species
subset_da_silv <- soc_cleaned_da_silva[soc_cleaned_da_silva$`Full Name` %in% unlist(intersects),]


subset_da_silv$`Full Name`
pruned.tree$tip.label


physorted_subset_da_silv= matrix(NA, nrow = 0, ncol =ncol(subset_da_silv ) )

for (i in 1:length(pruned.tree$tip.label)){
  
  physorted_subset_da_silv=rbind(physorted_subset_da_silv, subset_da_silv[pruned.tree$tip.label[[i]]==subset_da_silv$`unlist(soc_treedat$phy$tip.label)`,])
 
}



#physorted_subset_da_silv$`unlist(soc_treedat$phy$tip.label)`==pruned.tree$tip.label

soc_treedat=make.treedata(pruned.tree, subset_da_silv)

soc_treedat$phy



pdf("test_pruned_tree.pdf", width=20, height=20)

plot(soc_treedat$phy, cex=0.5)
nodelabels()

dev.off()


is.ultrametric(soc_treedat$phy)
soc_treedat$phy= force.ultrametric(soc_treedat$phy)

#write.tree(soc_treedat$phy, "pruned_dadsilva_tree.tre")

subset_da_silv=cbind(unlist(soc_treedat$phy$tip.label), soc_treedat$dat)

#Class species as either solitary, parasocial, advanced eusocial, or parasitic
soc_col <- paste(subset_da_silv[['Sociality']]) # get "raw" sociality values from df
split_soc_col <- strsplit(paste(soc_col), "/")
split_soc_col <- unlist(split_soc_col)
unique_socs <- unique(split_soc_col)

unique_socs=unique_socs[unique_socs!="solitary (uncertain, likelysubsocial)"]

unique_socs=c( "solitary", "subsocial", "semisocial", "communal", "primitively eusocial", "advanced eusocial","kleptoparasite", "social parasite", "ambiguous" )

soc_mat=matrix(0, nrow=length(soc_col), ncol=length(unique_socs))

colnames(soc_mat)=unique_socs
rownames(soc_mat)=subset_da_silv$`unlist(soc_treedat$phy$tip.label)`

for (soc in 1:length(unique_socs)){

  soc_mat[,soc] = as.integer(grepl(unique_socs[soc], soc_col))
  
}


#
###exploratory
#for (i in seq_along(soc_col)) {
#  if (any(grepl(klepto_cats, soc_col[i]))) {
#    print(subset_da_silv[['Full Name']][i])
#  }
#}



{
  pdf(file="new_Social_traits_tree.pdf", width = 75, height = 100)
  
  
  library(diversitree)
  col_pal=list(
    c("white", "black"   ),
    c("white", "light blue" ),
    c("white", "blue"   ),
    c("white", "dark blue"   ),
    c("white", "light green"    ),
    c("white", "dark green" ),
    c("white", "red"  ),
    c("white", "orange red"  ),
    c("white", "gray"))
  
  names(col_pal)=colnames(soc_mat)
  
  soc_df = as.data.frame(soc_mat[,1:9])
  
  class((soc_df$subsocial))
  names(soc_df)
  
  tree=soc_treedat$phy
  
  diversitree::trait.plot(tree = soc_treedat$phy, dat=as.data.frame(as.matrix(soc_df)), col_pal, cex.legend = 11,type = "p",w = 0.25, cex=5)
  
  dev.off()
}











class(soc_df[,1:3])

`[.data.frame`((soc_df[,1:3]), soc_treedat$phy$tip.label,  drop = FALSE) 

seq_along(soc_df[,1:3])










##Parasite categories
klepto_cats  <- c("kleptoparasite", "social parasite")
##Presocial categories
presoc_cats  <- c("solitary")
##Eusocial categories
eusoc_cats   <- c("primitively eusocial", "advanced eusocial")
##Parasocial categories
parasoc_cats <- c("communal", "semisocial", "subsocial")

simple_socs  = c("Presocial" , "Parasocial", "Eusocial", "Kleptoparasite")

simple_soc_mat=matrix(0, nrow=length(soc_col), ncol=length(simple_socs))

colnames(simple_soc_mat)=simple_socs
rownames(simple_soc_mat)=subset_da_silv$`unlist(soc_treedat$phy$tip.label)`


##exploratory

for (i in seq_along(soc_col)) {
  
  if(any(soc_mat[i,klepto_cats]==1)){
    
    simple_soc_mat[i,"Kleptoparasite"] = 1
    
  }
  
  if(any(soc_mat[i,presoc_cats]==1)){
    
    simple_soc_mat[i,"Presocial" ] = 1
    
  }
  
  
  if(any(soc_mat[i,eusoc_cats]==1)){
    
    simple_soc_mat[i,"Eusocial" ] = 1
    
  }
  
  
  if(any(soc_mat[i,parasoc_cats]==1)){
    
    simple_soc_mat[i,"Parasocial" ] = 1
    
  }
  
  
  
}



{
  pdf(file="Simple_Social_traits_tree.pdf", width = 75, height = 100)
  
  
  library(diversitree)
  col_pal=list(
    c("white", "black"   ),
    c("white", "blue"   ),
    c("white", "light green"    ),
    c("white", "red"  )
    )
  
  names(col_pal)=colnames(simple_soc_mat)
  
  simple_soc_df = as.data.frame(simple_soc_mat[,1:4])
  
  #class((simple_soc_df$subsocial))
  names(simple_soc_df)
  
  tree=soc_treedat$phy
  
  diversitree::trait.plot(tree = soc_treedat$phy, dat=as.data.frame(as.matrix(simple_soc_df)), col_pal, cex.legend = 11,type = "p",w = 0.25, cex=5)
  
  dev.off()
}


#aveRDS(soc_treedat, file = "soctreedat")


#just get the set of soc items that correspond to each cat and then do an intersect and see whether the intersects are zero or not


#########convert from simplie_soc matrix to nexus file


#full_PFA_df = read.csv(paste(path, "Data/", subclade_labels[[clade]],"_", tree_threshold_vec[[thresh]],"_thresh_", arbor_vec[[a]] ,"Arb", "_full.csv", sep=""))
#phy_PFA_df  = read.csv(paste(path, "Data/", subclade_labels[[clade]],"_", tree_threshold_vec[[thresh]],"_thresh_", arbor_vec[[a]] ,"Arb", "_phy.csv", sep=""))
#FullTree=read.tree(paste(path,"Data/preadapt_tree", sep=""))


#write.csv(phy_PFA_df,paste("Data/",tree_threshold_vec[[thresh]],"_thresh_preadapt_phy.csv", sep=""), row.names = F)
full_df = full_PFA_df[,trait_sets[[traits]]]
phy_df  =  phy_PFA_df[,trait_sets[[traits]]]

#modify arboreal to be truly discrete
#full_PFA_df$multi.arboreal[full_PFA_df$multi.arboreal %in% c("semi.arboreal", "arboreal")] =  1
#full_PFA_df$multi.arboreal[full_PFA_df$multi.arboreal %in% c("not.arboreal")]              =  0
#phy_PFA_df$multi.arboreal[phy_PFA_df$multi.arboreal %in% c("semi.arboreal", "arboreal")] =  1

#full_df_string=do.call(paste, lapply(1:ncol(full_df), function(char) full_df[,char] ))
#phy_df_string= do.call(paste, lapply(1:ncol(phy_df), function(char) phy_df[,char] ))

#samp_frac=table(phy_df_string)/table(full_df_string)

phy_PFA_df$Binomial=gsub(" ","_",  phy_PFA_df$Binomial)



tip_states=convert_traitdf2traitsIntvec(phy_df, state_space = state_space)

names(tip_states)=phy_PFA_df$Binomial

write_MuHiSSE_nex(tip_states, file = paste(path, "Data/", trait_set_labels[[traits]],"_",subclade_labels[[clade]],"_", tree_threshold_vec[[thresh]],"_thresh_", arbor_vec[[a]] ,"Arb", "_phy.nex", sep="")
                  , format = "standard", interleaved = F, nchars=length(trait_sets[[traits]]) )










