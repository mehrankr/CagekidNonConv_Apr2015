options(stringsAsFactors=FALSE)
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0){
	if(args[1]=="-h" | args[1]=="--help"){
	print("1. Enter a cutoff for AbHAC FDR. By default: 0.05")
	print("2. Enter the minimum number of patients a gene must be mutated in for consideration")
	}
}
if(length(args)==0){
	print("Using 0.05 as FDR cutoff of pathway analysis")
	cutoff = 0.05
	min_pats = 1
}else{
	cutoff = as.numeric(args[1])
	min_pats = as.numeric(args[2])
	if(is.na(cutoff) | is.na(min_pats)){
		stop("Non-numeric values not accepted for cutoff")
	}
}
require(KeggPA)
require(pheatmap)
require(ggplot2)
require(AbHAC)


load_files = function(){
	print(paste("loading files from ../data/tsv"))
	setwd("../data/tsv")
	files=dir()[grep("tsv.gz",dir())]
	list.files=list()
	for(i in 1:length(files)){
		list.files[[i]] = read.csv(files[i], header=T, sep="\t", row.names=1)
		names(list.files)[i] = unlist(strsplit(files[i], ".tsv", T))[1]
		print(paste(files[i], "loaded"))
	}
        list.files$Pap.Unclassified = rbind(list.files$unclassified, list.files$Pap)
	setwd("../../src")
	return(list.files)
}

get_genes = function(df, min_pats){
	genes = colnames(df)[which(apply(df,2,sum)>=min_pats)]
	return(genes)
}

abhac_analysis = function(genes, cutoff, file_name){
	proteins = ids.to.uniprot(genes)
	abhac_result = abhac.brief(snv=proteins, enrichment.categories="snv",
		fisher.fdr="Permutation", fisher.fdr.cutoff=cutoff)
	gz1 = gzfile(file_name, "w")
	write.table(abhac_result, file=gz1, quote=F, row.names=F, sep="\t")
	close(gz1)
	return(abhac_result)
}

get_abhac_genes = function(abhac_file){
	return(abhac_file[which(abhac_file[,2]<=0.5),4])
}

pathway_analysis = function(genes, project_name){
	if(any(genes%in%kegg.list[,3])){
		kpa = keggPA(hgnc=genes, project_name = project_name)
	}else{
		kpa = NA
	}
	return(kpa)
}

visualize_pathways = function(list.pathways, cutoff, file_name){
	ind = which(unlist(lapply(list.pathways, is.na)))
	ind = setdiff(1:length(list.pathways),ind)
	ind = ind[1]
	mat.sig.pathways = matrix(0,
		nrow = nrow(list.pathways[[ind]]),
		ncol = length(list.pathways),
		dimnames = c(list(rownames(list.pathways[[ind]])),
			list(names(list.pathways))))
	for(i in 1:length(list.pathways)){
		temp_kpa = list.pathways[[i]]
		if(!is.na(temp_kpa)){
			name_kpa = names(list.pathways)[i]
			pathways = rownames(temp_kpa)[which(temp_kpa$BH.FDR<=cutoff)]
			if(length(pathways)>0){
				mat.sig.pathways[pathways,i] = 10
			}
		}
	}
	mat.sig.pathways = mat.sig.pathways[apply(mat.sig.pathways>0,1,any),]
	if(nrow(mat.sig.pathways)>0){
		cutoff_name = paste(cutoff*100,"%", sep="")
		pdf(file_name, width = 5, height = nrow(mat.sig.pathways) * 0.5)
		pheatmap(mat.sig.pathways,
			cluster_rows=FALSE,
			cluster_cols=FALSE,
			legend=FALSE,
			color = c("gray90", "skyblue"),
			main = paste("FDR cutoff:",cutoff_name))
		dev.off()
	}
}

visualize_abhac_targets = function(list.abhac, cutoff, file_name){
	genes = unique(unlist(lapply(list.abhac,function(x){
		return(x[which(x[,2]<=0.5),4])
	})))
	mat.sig.targets = matrix(0,
                nrow = length(genes),
                ncol = length(list.abhac),
                dimnames = c(list(genes),
                        list(names(list.abhac))))
	for(i in 1:length(list.abhac)){
                temp_abhac = list.abhac[[i]]
                name_abhac = names(list.abhac)[i]
                targets = temp_abhac[which(temp_abhac[,2]<=0.05),4]
                if(length(targets)>0){
                        mat.sig.targets[targets,i] = 10
                }
        }
	cutoff_name = paste(cutoff*100,"%", sep="")
        pdf(file_name, width = 5, height = nrow(mat.sig.targets) * 0.5)
        pheatmap(mat.sig.targets,
                cluster_rows=FALSE,
                cluster_cols=FALSE,
                legend=FALSE,
                color = c("gray90", "skyblue"),
                main = paste("FDR cutoff:",cutoff_name))
        dev.off()
}

list.files = load_files()
setwd("../results/AbHACNetworkAnalysis")
dir_name = paste("FDR",paste(cutoff*100,"p", sep=""),"NumberPatients", min_pats, sep="_")
dir.create(dir_name)
setwd(dir_name)
list.abhac = vector("list", length(list.files))
for(i in 1:length(list.files)){
	if(min_pats==1){
		project_name = paste("abhac_analysis",names(list.files)[i], min_pats, "minimum_patient_affected.tsv.gz", sep="_")
	}else{
		project_name = paste("abhac_analysis",names(list.files)[i], min_pats, "minimum_patients_affected.tsv.gz", sep="_")
	}
	list.abhac[[i]] = abhac_analysis(get_genes(list.files[[i]], min_pats),cutoff, project_name)
	names(list.abhac)[i] = names(list.files)[i]
}

list.pathways = vector('list', length(list.abhac))
for(i in 1:length(list.abhac)){
	abhac_file = list.abhac[[i]]
	if(min_pats==1){
                project_name = paste(names(list.files)[i], min_pats, "minimum_patient_affected", sep="_")
        }else{
                project_name = paste(names(list.files)[i], min_pats, "minimum_patients_affected", sep="_")
        }

	list.pathways[[i]] = pathway_analysis(get_abhac_genes(abhac_file), project_name)
	names(list.pathways)[i] = names(list.abhac)[i]
}


file_name = paste("AbHACNetworkAnalysis_pathway_analysis_Heatmap_",min_pats,"patients", paste(cutoff*100,"p", sep=""), "cutoff.pdf", sep="_")
visualize_pathways(list.pathways, cutoff, file_name)
file_name = paste("AbHACNetworkAnalysis_targets_Heatmap_",min_pats,"patients", paste(cutoff*100,"p", sep=""), "cutoff.pdf", sep="_")
visualize_abhac_targets(list.abhac, cutoff, file_name)
print(paste("successfully created the results in:", getwd()))

