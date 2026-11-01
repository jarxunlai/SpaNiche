#' prepare_for_LRintegratedmatrix
#'
#' @description
#' This function prepares an integrated expression matrix for ligand–receptor analysis in spatial transcriptomics data.
#' It filters CellChat ligand–receptor interactions based on genes present in the Spatial assay of a Seurat object and returns both the integrated expression matrix and the filtered interaction table.
#' @param seu_obj A Seurat object containing spatial transcriptomics data, with 'Spatial' as the assay name. The object should be normalized using the NormalizeData function.
#' @param ccc_db This parameter represents the reference database of cell-cell interactions, e.g., CellChatDB.human, containing interaction, complex, and geneInfo tables.
#' @param assay Assay to read normalized expression from. Defaults to
#'   \code{"Spatial"} if present, otherwise the default assay.
#'
#' @return A list with two elements: (1) a gene/complex-by-spot expression matrix, and (2) a filtered CellChat interaction data frame.
#' @export
#'
prepare_for_LRintegratedmatrix=function(
    seu_obj,
    ccc_db,
    assay = NULL
){
  expr = .get_layer_data(seu_obj, layer = "data", assay = assay)
  wehave_genes=rownames(expr)


  ### Extract reference information from the CellChat database
  cc_interaction=ccc_db$interaction
  cc_complex=ccc_db$complex
  cc_geneInfo=ccc_db$geneInfo
  # Identify how many ligands/receptors are complexes
  complex1=setdiff(union(cc_interaction$ligand, cc_interaction$receptor) ,cc_geneInfo$Symbol)


  ### Filter the reference information based on genes present in ST data
  cc_geneInfo=cc_geneInfo%>%dplyr::filter(Symbol %in% wehave_genes)
  tmp_complex=cc_complex%>%apply(2,function(x){x %in% c(wehave_genes,"")})
  cc_complex=cc_complex[!rowSums(tmp_complex) < 4,]
  # Number of candidate complexes after filtering and the gene composition matrix of complexes
  complex1=intersect(complex1,rownames(cc_complex))
  cc_complex=cc_complex[complex1,]


  ### Extract genes
  # Genes involved in complexes
  gene1=unlist(cc_complex)
  gene1=gene1[gene1 != ""]
  gene1=as.character(gene1)
  # Genes involved in non-complex ligands and receptors
  gene2=intersect(union(cc_interaction$ligand, cc_interaction$receptor) ,cc_geneInfo$Symbol)
  # All genes involved in CellChat interactions, constrained by the current ST dataset
  gene.used=union(gene1,gene2)


  ### Extract expression data
  gene.used=intersect(rownames(expr),gene.used) # Genes present in the ST dataset; adjust the order accordingly
  spatial.data=as.matrix(expr[gene.used,]) %>% as.data.frame() # Subset the expression matrix to relevant genes


  ### Add expression data for complexes
  new.data=matrix(nrow = length(complex1),ncol = dim(seu_obj)[2]) # Expression matrix for complexes

  for (i in 1:length(complex1)) {
    one_complex=cc_complex[complex1[i],] %>% as.character()
    one_complex=one_complex[one_complex!=""]

    one_complex_mean=colMeans(spatial.data[one_complex,])
    new.data[i,] = one_complex_mean
  }
  new.data=as.data.frame(new.data)
  colnames(new.data)=colnames(spatial.data)
  rownames(new.data)=complex1

  new.data=rbind(new.data,spatial.data[gene2,])


  ### Filter the original interaction data frame
  cc_interaction=cc_interaction[
    cc_interaction$ligand %in% rownames(new.data) &
      cc_interaction$receptor %in% rownames(new.data),
    ]


  ### Export results
  list(new.data,cc_interaction)
}
