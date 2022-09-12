## Complete Pipeline for Network Construction from Transcriptome Dataset ##
setwd('e:/writing/NetBID/')
############### Step 0: 准备 ###############
# 安装包，如果不顺利，在保持原有R版本的基础上，下载最新版本的R，调用后一般可以安装好
# 装包过程中避免出些一切warning错
# 我才用到github下载后本地安装
#devtools::install_local('d:/R/NetBID-master.zip')
# 加载包
library(NetBID2)

# 设定主要目录和项目名称，在目录下创建好
project_main_dir <- './test' 

# 工作日期，可以避免重复文件名
current_date <- format(Sys.time(), "%Y-%m-%d") 
project_name <- sprintf('project_%s',current_date) 

# 如果进入以前的，可以用以下，和上面那句选择其一运行
#project_name <- sprintf('project_%s',"2022-03-15") 


# 创建一个文件夹
network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)

############### Step 1: 数据集下载和载入 ###############

# GEO下载，同时下载GPL，注意会有点慢，因为我们一般是不下载GPL的
net_eset <- load.exp.GEO(out.dir=network.par$out.dir.DATA,GSE='GSE116028',GPL='GPL6480',getGPL=TRUE,update=FALSE)

# ID转换，依然是必须的，NetBID提供了快捷的转ID的方式
net_eset <- update_eset.feature(use_eset=net_eset,use_feature_info=fData(net_eset),from_feature='ID',to_feature='GENE_SYMBOL',merge_method='median')
# NetBID还提供了把pData表型信息附加到对象中的方法
net_eset <- update_eset.phenotype(use_eset=net_eset,use_phenotype_info=pData(net_eset),use_sample_col='geo_accession',use_col='GEO-auto')

# net_eset加入到network.par对象中，关键步骤
network.par$net.eset <- net_eset

# 标准化前的质控
# 出去文件夹看QC文件，可以看到这个时候表达量是不平均的
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=NULL,
             do.logtransform=FALSE,prefix='beforeQC_',
             pre_define=c('WNT'='blue','SHH'='red','G4'='green'),emb_plot_type='2D.interactive')

# 保存文件（载入步骤）
NetBID.saveRData(network.par = network.par,step='exp-load')

############### Step 2: 数据标准化 ###############

# 读入 Step 1
NetBID.loadRData(network.par = network.par,step='exp-load')

# 表达矩阵
mat <- exprs(network.par$net.eset)

# boxplot看一下，并不齐
boxplot(mat)
# 处理缺失值，本数据集没有缺失值
sample_na_count <- apply(mat,1,function(x){length(which(is.na(x)==TRUE))})
print(table(sample_na_count))
gene_na_count <- apply(mat,2,function(x){length(which(is.na(x)==TRUE))})
print(table(gene_na_count))
# impute.knn填补缺失值
if(sum(sample_na_count)+sum(gene_na_count)>0) mat <- impute.knn(mat)$data

## log2转换，当然是选择性的，下面代码自动判断
med_val <- median(apply(mat,2,median)); print(med_val)
if(med_val>16){mat <- log2(mat)}

## 四分位数标准化，我们很熟悉了，相当于normalizebetweenarray
mat <- normalizeQuantiles(mat)
# 齐了
boxplot(mat)
## 过滤到低表达基因
choose1 <- apply(mat<= quantile(mat, probs = 0.05), 1, sum)<= ncol(mat) * 0.90
print(table(choose1))
mat <- mat[choose1,]

# 用标化后的矩阵换掉原来的
net_eset <- generate.eset(exp_mat=mat, phenotype_info=pData(network.par$net.eset)[colnames(mat),],
                          feature_info=fData(network.par$net.eset)[rownames(mat),],
                          annotation_info=annotation(network.par$net.eset))
# 更新network.par
network.par$net.eset <- net_eset

# 再次画QC图，标化后的
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='afterQC_',
             pre_define=c('WNT'='blue','SHH'='red','G4'='green'),emb_plot_type='2D.interactive')

# 保存
NetBID.saveRData(network.par = network.par,step='exp-QC')

############### Step 2: 数据标准化 ###############

# 读入 Step 1
NetBID.loadRData(network.par = network.par,step='exp-load')

# 表达矩阵
mat <- exprs(network.par$net.eset)

# boxplot看一下，并不齐
boxplot(mat)
# 处理缺失值，本数据集没有缺失值
sample_na_count <- apply(mat,1,function(x){length(which(is.na(x)==TRUE))})
print(table(sample_na_count))
gene_na_count <- apply(mat,2,function(x){length(which(is.na(x)==TRUE))})
print(table(gene_na_count))
# impute.knn填补缺失值
if(sum(sample_na_count)+sum(gene_na_count)>0) mat <- impute.knn(mat)$data

## log2转换，当然是选择性的，下面代码自动判断
med_val <- median(apply(mat,2,median)); print(med_val)
if(med_val>16){mat <- log2(mat)}

## 四分位数标准化，我们很熟悉了，相当于normalizebetweenarray
mat <- normalizeQuantiles(mat)
# 齐了
boxplot(mat)
## 过滤到低表达基因
choose1 <- apply(mat<= quantile(mat, probs = 0.05), 1, sum)<= ncol(mat) * 0.90
print(table(choose1))
mat <- mat[choose1,]

# 用标化后的矩阵换掉原来的
net_eset <- generate.eset(exp_mat=mat, phenotype_info=pData(network.par$net.eset)[colnames(mat),],
                          feature_info=fData(network.par$net.eset)[rownames(mat),],
                          annotation_info=annotation(network.par$net.eset))
# 更新network.par
network.par$net.eset <- net_eset

# 再次画QC图，标化后的
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='afterQC_',
             pre_define=c('WNT'='blue','SHH'='red','G4'='green'),emb_plot_type='2D.interactive')

# 保存
NetBID.saveRData(network.par = network.par,step='exp-QC')

############### Step 3: 样本聚类 ###############

# 载入 Step 2
NetBID.loadRData(network.par = network.par,step='exp-QC')

# 选择高变基因
mat <- exprs(network.par$net.eset)
choose1 <- IQR.filter(exp_mat=mat,use_genes=rownames(mat),thre = 0.5)
print(table(choose1))
mat <- mat[choose1,]

# 临时的eset
tmp_net_eset <- generate.eset(exp_mat=mat, phenotype_info=pData(network.par$net.eset)[colnames(mat),],
                              feature_info=fData(network.par$net.eset)[rownames(mat),], annotation_info=annotation(network.par$net.eset))
# QC plot
draw.eset.QC(tmp_net_eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='Cluster_',
             pre_define=c('WNT'='blue','SHH'='red','G4'='green'),emb_plot_type='2D.interactive')

# 看样本好不好
# 抽提表型信息
phe <- pData(network.par$net.eset)
intgroup <- get_int_group(network.par$net.eset)
# 2维聚类
for(i in 1:length(intgroup)){
  print(intgroup[i])
  pred_label <- draw.emb.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,intgroup[i]),
                                pre_define=c('WNT'='blue','SHH'='red','G4'='green'))
}
# 3维聚类
for(i in 1:length(intgroup)){
  print(intgroup[i])
  pred_label <- draw.emb.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,intgroup[i]),plot_type='3D')
  print(table(list(pred_label=pred_label,obs_label=get_obs_label(phe,intgroup[i]))))
}
# Pick one phenotype column "subgroup" from the demo eset and show various plots NetBID2 can create
use_int <- 'subgroup'
pred_label <- draw.emb.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='2D',pre_define=c('WNT'='blue','SHH'='red','G4'='green'))
pred_label <- draw.emb.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='2D.ellipse',pre_define=c('WNT'='blue','SHH'='red','G4'='green'))
pred_label <- draw.emb.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='2D.text',pre_define=c('WNT'='blue','SHH'='red','G4'='green'))
pred_label <- draw.emb.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='3D',pre_define=c('WNT'='blue','SHH'='red','G4'='green'))
print(table(list(pred_label=pred_label,obs_label=get_obs_label(phe, use_int))))
draw.clustComp(pred_label,obs_label=get_obs_label(phe,use_int),outlier_cex=1,low_K=10)

############### Step 4: 准备 SJARACNe 所需的文件 ###############

# 载入 Step 2
NetBID.loadRData(network.par = network.par, step='exp-QC') ## 不要载入第三步cluster的文件哦！
load('./test/project_2022-03-15/DATA/network.par.Step.exp-QC.RData')

# 载入数据库，都是内置的
db.preload(use_level='gene',use_spe='human',update=FALSE)

# gene ID得到转录因子/信号通路列表中
use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)

# 挑选样本，这边13个全上
phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) # here is using all samples, users can modify
prj.name <- network.par$project.name # if use different samples, need to change the project name
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0.5,IQR.loose_thre = 0.1,
                 SJAR.project_name=prj.name,SJAR.main_dir=network.par$out.dir.SJAR)


tf=read.table('test/project_2022-03-15/SJAR/project_2022-03-15/tf.txt',header = F)
tf$V1=paste0('MSC_',tf$V1)
write.table(tf,file='test/project_2022-03-15/SJAR/project_2022-03-15/demo1_tf.txt',col.names = F,row.names = F,quote = F)
sig=read.table('test/project_2022-03-15/SJAR/project_2022-03-15/sig.txt',header = F)
sig$V1=paste0('MSC_',sig$V1)
write.table(sig,file='test/project_2022-03-15/SJAR/project_2022-03-15/MSC_sig.txt',col.names = F,row.names = F,quote = F)

exp=read.table('test/project_2022-03-15/SJAR/project_2022-03-15/MSC_input.exp',header = T,check.names = F)
exp$isoformId=paste0('MSC_',exp$isoformId)
write.table(exp,file='test/project_2022-03-15/SJAR/project_2022-03-15/MSC_input.exp',col.names = T,row.names = F,quote = F,sep = '\t')

exp=read.table('./BRCA100.exp',header = T,check.names = F,sep = '\t')
