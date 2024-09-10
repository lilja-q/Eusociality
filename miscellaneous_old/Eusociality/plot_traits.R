library(ape)
library(phytools)
library(stringr)
library(treeplyr)


#tree=read.tree("BEE_mat7_fulltree.nwk")
traits=read.csv("cleaned_da_silva_SWM_mod.csv", header = F)

tree=read.tree("BEE_mat7_fulltree_tplo35_sf20lp.nwk")
tree_speces_mat=str_split_fixed(tree$tip.label, pattern = "_", 4)[,1:2]
tree_species=apply(tree_speces_mat, 1, function(x){ paste0(x, collapse="_")})

tree$tip.label=tree_species


tree_dat=make.treedata(tree, traits)



social_traits=unique(as.vector(str_split_fixed(tree_dat$dat$V2, "/", n=10)))
social_traits=social_traits[social_traits!=""]

#create trait matrix, each column is a sociality state, if a species has it, then we socre 1, if not its 0, we create the empty matrix here and fill it in with the loop below
trait_mat=matrix(0, nrow = length(tree_dat$phy$tip.label), ncol=length(social_traits))
colnames(trait_mat)=social_traits

#split of the sociality score to a bunch of strings (some of these are janky + spelling errors)
trait_stringmat=str_split_fixed(tree_dat$dat$V2, "/", n=15)

rownames(trait_stringmat) = tree_dat$phy$tip.label
rownames(trait_mat)       = tree_dat$phy$tip.label

# check to make sure its a kelptoparasite
trait_stringmat["Epeolus_pusillus",]


i=1



tree_dat$phy$tip.label


for ( i in 1:ncol(trait_mat)){
  colnames(trait_mat)[[i]]
  
  trait_mat[(1:nrow(trait_stringmat))[grep(social_traits[[i]] ,x = tree_dat$dat$V2, pattern = T)],i]=1
}


#trait_mat[trait_mat[,"primitively eusocial"]==1,"primatively eusocial" ]=1
#trait_mat[trait_mat[,"primitively eusocial?"]==1,"primatively eusocial?" ]=1
#trait_mat[trait_mat[,"communal"]==1," communal" ]=1

social_traits_ordered=c(
"solitary"   ,                                               
#"solitary (uncertain, likely subsocial)",
"subsocial"   ,
#"subsocial?"  , 
"semisocial"  ,
#"semisocial?" ,
"communal"    ,
"primitively_eusocial" ,
#"primitively eusocial?",
"advanced_eusociality" ,   
"social_parasite"  ,
"kleptoparasite"  ,
"?"
)                          

final_trait_mat=as.data.frame(trait_mat[,social_traits_ordered])

colnames(final_trait_mat)=social_traits_ordered

rownames(final_trait_mat)=tree_dat$phy$tip.label


final_trait_mat$`?`

Solitary=c("solitary")
Solitary_sp=Reduce(`+`, lapply(Solitary, function(i) final_trait_mat[i]))

Social=c("solitary (uncertain, likely subsocial)", "subsocial", "subsocial?", "semisocial", "semisocial?", "communal", "primitively eusocial", "primitively eusocial?")

Social_sp=Reduce(`+`, lapply(Social, function(i) final_trait_mat[i]))

Advanced_Eusocial=c("advanced eusociality")
A_Eusocial_sp=Reduce(`+`, lapply(Advanced_Eusocial, function(i) final_trait_mat[i]))

Parasite=c("social parasite", "kleptoparasite")
Parasite_sp=Reduce(`+`, lapply(Parasite, function(i) final_trait_mat[i]))


cbind(Solitary_sp, Social_sp, A_Eusocial_sp, Parasite_sp)

#Parasocial=c("semisocial", "solitary (uncertain, likely subsocial)", "subsocial", "subsocial?" )


do.call(sum, lapply(Solitary, function(i) final_trait_mat[i]), )>=1

sum(lapply(Scoail, function(i) final_trait_mat[i]))>=1

Reduce(`+`, lapply(Solitary, function(i) final_trait_mat[i]))

as.numeric((final_trait_mat$solitary +final_trait_mat$`solitary (uncertain, likely subsocial)` +final_trait_mat$subsocial+final_trait_mat$`subsocial?`)>=1)

as.numeric((final_trait_mat$semisocial +final_trait_mat$`solitary (uncertain, likely subsocial)` +final_trait_mat$subsocial+final_trait_mat$`subsocial?`)>=1)


final_trait_mat["Epeolus_pusillus",]

{
  pdf(file="Social_traits_tree.pdf", width = 450, height = 300)
  
  
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
  names(col_pal)=colnames(final_trait_mat)
  
  diversitree::trait.plot(tree_dat$phy, dat=final_trait_mat,col_pal,legend = T,cex.lab = 10, cex.legend = 11,type = "p",w = 0.25 )
  
  dev.off()
}



#######bad code/scrap######


library(ape)
library(phytools)
library(stringr)


tree=read.tree("BEE_mat7_fulltree.nwk")
traits=read.csv("cleaned_da_silva (1).csv", header = F)


tree=read.tree("BEE_mat7_fulltree_tplo35_sf20lp.nwk")

tree_speces_mat=str_split_fixed(tree$tip.label, pattern = "_", 4)[,1:2]

tree_species=apply(tree_speces_mat, 1, function(x){ paste0(x, collapse="_")})

trait_species=traits$V1

shared_sp=intersect(tree_species, trait_species)

shared_traits=traits[traits$V1 %in% shared_sp,]

unique(shared_traits$V2)

shared_sp

pruned_tree=drop.tip( tree, tree$tip.label[!(tree_species %in% shared_sp)])

#plot(pruned_tree)
social_traits=unique(as.vector(str_split_fixed(shared_traits$V2, "/", n=10)))[-6]

trait_mat=matrix(0, nrow = length(pruned_tree$tip.label), ncol=length(social_traits))
colnames(trait_mat)=social_traits

trait_stringmat=str_split_fixed(shared_traits$V2, "/", n=15)

for ( i in 1:ncol(trait_mat)){
  trait_mat[(1:nrow(trait_stringmat))[grep(social_traits[[i]] ,x = shared_traits$V2,)],i]=1
}

trait_mat[trait_mat[,7]==1,5]=1
trait_mat[trait_mat[,3]==1,12]=1

final_trait_mat=as.data.frame(trait_mat[,c(1,10, 3,14, 15, 4, 9, 5, 12, 6, 2,8)])

colnames(final_trait_mat)=social_traits[c(1,10, 3,14, 15, 4, 9, 5, 12, 6, 2,8)]


rownames(final_trait_mat)=pruned_tree$tip.label

{
  pdf(file="Social_traits_tree.pdf", width = 450, height = 300)
  
  
  library(diversitree)
  col_pal=list(
    c("white", "black"   ),
    c("gray", "black"   ),
    c("white", "light blue" ),
    c("gray", "light blue"   ),
    c("white", "blue"   ),
    c("white", "dark blue"   ),
    c("gray", "dark blue"   ),
    c("white", "light green"    ),
    c("gray", "light green"    ),
    c("white", "dark green" ),
    c("white", "red"  ),
    c("white", "orange red"  ))
  names(col_pal)=colnames(final_trait_mat)
  
  diversitree::trait.plot(pruned_tree, dat=final_trait_mat,col_pal,legend = T,cex.lab = 10, cex.legend = 11,type = "p",w = 0.25 )
  
  dev.off()
}









