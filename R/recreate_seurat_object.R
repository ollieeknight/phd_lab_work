alldata <- readRDS('RA_atlas/objects/alldata.rds')
alldata <- JoinLayers(alldata, assay = 'RNA', layer = c('counts', 'data', 'scale.data'))
alldata <- JoinLayers(alldata, assay = 'ADT', layer = c('counts', 'data', 'scale.data'))

metadata <- alldata@meta.data

rna_counts <- LayerData(alldata, assay = 'RNA', layer = 'counts')
rna_data <- LayerData(alldata, assay = 'RNA', layer = 'data')
adt_counts <- LayerData(alldata, assay = 'ADT', layer = 'counts')
adt_data <- LayerData(alldata, assay = 'ADT', layer = 'data')

alldata <- CreateSeuratObject(counts = rna_counts, data = rna_data, meta.data = metadata)
alldata[['ADT']] <- CreateAssay5Object(counts = adt_counts, data = adt_data)

alldata <- SetIdent(alldata, value = 'celltype_l1')
DefaultAssay(alldata) <- 'RNA'

rm(rna_counts, rna_data, adt_counts, adt_data)

gc()