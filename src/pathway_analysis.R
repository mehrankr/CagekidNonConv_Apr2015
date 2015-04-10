options(stringsAsFactors=FALSE)
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0){
	if(args[1]=="-h" | args[1]=="--help"){
	print("1. Enter a cutoff for kegg pathway analysis FDR. By default: 0.05")
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
	setwd("../../src")
	return(list.files)
}

get_genes = function(df, min_pats){
	genes = colnames(df)[which(apply(df,2,sum)>=min_pats)]
	return(genes)
}


pathway_analysis = function(genes, project_name){
	kpa = keggPA(hgnc=genes, project_name = project_name)
	return(kpa)
}

visualize_pathways = function(list.pathways, cutoff, file_name){
	mat.sig.pathways = matrix(0,
		nrow = nrow(list.pathways[[1]]),
		ncol = length(list.pathways),
		dimnames = c(list(rownames(list.pathways[[1]])),
			list(names(list.pathways))))
	for(i in 1:length(list.pathways)){
		temp_kpa = list.pathways[[i]]
		name_kpa = names(list.pathways)[i]
		pathways = rownames(temp_kpa)[which(temp_kpa$BH.FDR<=cutoff)]
		if(length(pathways)>0){
			mat.sig.pathways[pathways,i] = 10
		}
	}
	mat.sig.pathways = mat.sig.pathways[apply(mat.sig.pathways>0,1,any),]
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

	

list.files = load_files()
setwd("../results/KeggPathwayAnalysis")
dir_name = paste("FDR",paste(cutoff*100,"p", sep=""),"NumberPatients", min_pats, sep="_")
dir.create(dir_name)
setwd(dir_name)
list.pathways = vector("list", length(list.files))
for(i in 1:length(list.files)){
	if(min_pats==1){
		project_name = paste(names(list.files)[i], min_pats, "minimum_patient_affected", sep="_")
	}else{
		project_name = paste(names(list.files)[i], min_pats, "minimum_patients_affected", sep="_")
	}
	list.pathways[[i]] = pathway_analysis(get_genes(list.files[[i]], min_pats), project_name)
	names(list.pathways)[i] = names(list.files)[i]
}


file_name = paste("KeggPathwayAnalysis_Heatmap_",min_pats,"patients", paste(cutoff*100,"p", sep=""), "cutoff.pdf", sep="_")
visualize_pathways(list.pathways, cutoff, file_name)
print(paste("successfully created the results in:", getwd()))
