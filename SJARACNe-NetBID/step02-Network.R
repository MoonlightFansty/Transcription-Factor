## Complete Pipeline for Driver Estimation and Master Table Creation ##
setwd('e:/writing/NetBID/')
############### Step 0: 准备###############

#加载R包
library(NetBID2)

# 使用示例network
network.dir <- sprintf('%s/demo1/network/',system.file(package = "NetBID2")) # use demo network in the package
network.project.name <- 'project_2019-02-14'# 若用示例数据，则不要改
# 定义主要工作目录
project_main_dir <- 'test/' # 依然在之前的test目录下
current_date <- format(Sys.time(), "%Y-%m-%d") # 防止重复
project_name <- sprintf('driver_%s',current_date) 

# 创立文件夹
# 关键步骤
# 这一步尤其难点的是，要注意D:/R/R_Library/NetBID2/demo1/network/下一定放着
# /SJAR/project_2019-02-14/output_tf_sjaracne_project_2019-02-14_out_.final/consensus_network_ncol_.txt
# 要对应日期并且要格式完全一样，才能保证读取
# 等于说network.dir和network.project.name是连锁到SJAR的两个project的
# 不理解清楚这一点，你很难将他们换成自己的数据
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

############### Step 1: Load in gene expression dataset for analysis (exp-load, exp-cluster, exp-QC) ###############
#sprintf的用法
#name <- 'max'
#sprintf('my name is %s',name)
#[1] "my name is max"

# 虽然我们已经制作了自己的表达矩阵（在第一节），但是此处用他的先
# 所以load进来的是内置包里面的exp-QC
# 我们自己的exp-QC在NetBID\test\project_2022-03-15\DATA下
load(sprintf('%s/DATA/network.par.Step.exp-QC.RData',network.dir)) # RData saved after QC in the network construction step
# 把exp-QC给到analysis.par中
analysis.par$cal.eset <- network.par$net.eset

# Save Step 1 network.par as RData，在\test\driver_2022-03-23\DATA中
NetBID.saveRData(analysis.par=analysis.par,step='exp-QC')

# 现在的analysis-par里面，有QC后的表达矩阵，也有network，他们都来源于内置包
# 重申一遍，我们将在下节换成自己运行的表达矩阵和network

############### Step 2: 读入网络文件并计算驱动活动（act-get） ###############

# 读取exp-QC
#analysis.par <- list()
#analysis.par$out.dir.DATA <- 'test//driver_2019-05-06//DATA/'
NetBID.loadRData(analysis.par=analysis.par,step='exp-QC')

# 获取网络
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)

# network质控
draw.network.QC(analysis.par$tf.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='TF_net_',html_info_limit=FALSE)
draw.network.QC(analysis.par$sig.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='SIG_net_',html_info_limit=TRUE)


# 把TF和SIG合并
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)

# 根据合并后的，获取活性
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')

# 制作活性矩阵
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

# 活性矩阵质控
# 这几部质控你需要理解的是，QC自动生成带有临床数据的质控图
# WNT,SHH,G4为分子亚型，需要去看下原文
?draw.eset.QC
draw.eset.QC(analysis.par$merge.ac.eset,outdir=analysis.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='AC_',
             pre_define=c('WNT'='blue','SHH'='red','G4'='green'),emb_plot_type ='2D.interactive')


# 保存活性矩阵
NetBID.saveRData(analysis.par=analysis.par,step='act-get')

############### Step 3: 获取driver的差异表达 (DE) / 差异活性 (DA) (act-DA) ###############

# 把活性矩阵读进来
#analysis.par <- list()
#analysis.par$out.dir.DATA <- 'test//driver_2019-05-06//DATA/'
NetBID.loadRData(analysis.par=analysis.par,step='act-get')

# 来个空列表
analysis.par$DE <- list()
analysis.par$DA <- list()

# 先来比较G4和WNT
comp_name <- 'G4.Vs.WNT' # Each comparison must has a name
# Get sample names from each compared group
phe_info <- pData(analysis.par$cal.eset)
# G1,G0为分组
G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`subgroup`=='WNT')] # Control group
# 差异分析，包括差异表达（表达矩阵，针对所有基因）和差异活性（活性矩阵，针对TF和SIG）
# 当然差异表达就是我们熟知的
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='G4',G0_name='WNT')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='G4',G0_name='WNT')
# 把差异表达和差异活性的结果保存入analysis.par
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid

# 比较G4和SHH
comp_name <- 'G4.Vs.SHH' # Each comparison must has a name
# Get sample names from each compared group
phe_info <- pData(analysis.par$cal.eset)
G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`subgroup`=='SHH')] # Control group
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='G4',G0_name='SHH')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='G4',G0_name='SHH')
# 把差异表达和差异活性的结果保存入analysis.par
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid

## G4 和其他类别(三种不同的比较方式)
# Combine the comparison results from `G4.Vs.WNT` and `G4.Vs.SHH`
comp_name <- 'G4.Vs.otherTwo' # Each comparison must has a name
DE_gene_comb <- combineDE(DE_list=list('G4.Vs.WNT'=analysis.par$DE$`G4.Vs.WNT`,'G4.Vs.SHH'=analysis.par$DE$`G4.Vs.SHH`))
DA_driver_comb <- combineDE(DE_list=list('G4.Vs.WNT'=analysis.par$DA$`G4.Vs.WNT`,'G4.Vs.SHH'=analysis.par$DA$`G4.Vs.SHH`))
analysis.par$DE[[comp_name]] <- DE_gene_comb$combine
analysis.par$DA[[comp_name]] <- DA_driver_comb$combine


# 画差异最大基因的图表
draw.combineDE(DE_gene_comb)
# 输出
draw.combineDE(DE_gene_comb,pdf_file=sprintf('%s/combineDE.pdf',analysis.par$out.dir.PLOT)) # Save it as PDF

# 画最大差异活性driver图表
draw.combineDE(DA_driver_comb)

# 另一种方式比较G4和其他，就是不合并，重新命名比较
comp_name <- 'G4.Vs.others'
phe_info <- pData(analysis.par$cal.eset)

G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')] # Experiment group
G0  <- rownames(phe_info)[which(phe_info$`subgroup`!='G4')] # Combine other groups as the Control group
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='G4',G0_name='others')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='G4',G0_name='others')
# 存入
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid

# 保存所有差异
NetBID.saveRData(analysis.par=analysis.par,step='act-DA')

# 可视化 top drivers，高级
draw.NetBID(DA_list=analysis.par$DA,DE_list=analysis.par$DE,main_id='G4.Vs.others')
draw.NetBID(DA_list=analysis.par$DA,DE_list=analysis.par$DE,main_id='G4.Vs.others',pdf_file=sprintf('%s/NetBID_TOP.pdf',analysis.par$out.dir.PLOT),text_cex=0.8) # Save as PDF

############### Step 4:为driver生成 master table ###############

# 载入
#analysis.par <- list()
#analysis.par$out.dir.DATA <- 'test//driver_自己的时间//DATA/'
NetBID.loadRData(analysis.par=analysis.par,step='act-DA')

# Reload data into R workspace, and saves it locally under db/ directory with specified species name and analysis level.
db.preload(use_level='gene',use_spe='human',update=FALSE)
# 把所有比较组列出
all_comp <- names(analysis.par$DE) # Users can use index or name to get target ones
# 准备换算表
use_genes <- unique(c(analysis.par$merge.network$network_dat$source.symbol,analysis.par$merge.network$network_dat$target.symbol))
transfer_tab <- get_IDtransfer2symbol2type(from_type = 'external_gene_name',use_genes=use_genes)
analysis.par$transfer_tab <- transfer_tab
# 创建主表
analysis.par$final_ms_tab <- generate.masterTable(use_comp=all_comp,DE=analysis.par$DE,DA=analysis.par$DA,
                                                  target_list=analysis.par$merge.network$target_list,
                                                  tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                                                  main_id_type='external_gene_name')

# 输出文件名
out_file <- sprintf('%s/%s_ms_tab.xlsx',analysis.par$out.dir.DATA,analysis.par$project.name)
# Highlight marker genes（要看原文，你感兴趣哪些基因）
mark_gene <- list(WNT=c('WIF1','TNC','GAD1','DKK2','EMX2'),
                  SHH=c('PDLIM3','EYA1','HHIP','ATOH1','SFRP1'),
                  G4=c('KCNA1','EOMES','KHDRBS2','RBM24','UNC5D'))
# 配色
#mark_col <- get.class.color(names(mark_gene)) # this randomly assign color codes
mark_col <- list(G4='green','WNT'='blue','SHH'='red')
# Save the final master table as EXCEL file
out2excel(analysis.par$final_ms_tab,out.xlsx = out_file,mark_gene,mark_col)

# 保存
NetBID.saveRData(analysis.par=analysis.par,step='ms-tab')

## Complete Pipeline for Driver Estimation and Master Table Creation ##
setwd('e:/writing/NetBID/')
############### Step 0: 准备###############

#加载R包
library(NetBID2)

# 不再用示例network，用自己输出的两个200多M的文件

## 文件放置架构（最重要）
## 模仿以下结构，在自己的test文件下，分别作network，SJAR,project，output_sig_sjaracne_project文件
## 这里是非常死板的！
## D:\R\R_Library\NetBID2\demo1\network\SJAR\project_2019-02-14\output_sig_sjaracne_project_2019-02-14_out_.final
## 一定要再三理解这里的逻辑！
## 最后存放的是服务器SJAR输出的TF和sig调控网络
## 2022-03-21是告诉别人/自己是什么时间得到这两个文件的
## 不理解的一定要去你装NetBID2这个包的文件看demo1数据里的结构
## D:\R\R_Library\NetBID2\demo1\network\SJAR\project_2019-02-14\output_sig_sjaracne_project_2019-02-14_out_.final
## network存放SJAR，SJAR不可省，SJAR存放project_时间，再里面存放输出文件
## 输出文件TF和sig只要是服务器里的都为consensus_network_ncol_.txt，不需要重命名
# 如果建好了，你将会非常容易理解下面的内容
network.dir <- './test/network/'
network.project.name <- 'project_2022-03-21'# 若用示例数据，则不要改
# 定义主要工作目录，这个不变，是我们做工作下的test文件
project_main_dir <- 'test/' # 依然在之前的test目录下
# 防止重复，例如3月23号我们用了示例数据，今天用自己跑的数据，就避免重复了
# 这个project时间不是network创立的时间，是analysis的时间
current_date <- format(Sys.time(), "%Y-%m-%d") # 防止重复
project_name <- sprintf('driver_%s',current_date) 

# 修改一下ID
tf=read.table('./test/network/SJAR/project_2022-03-21/output_tf_sjaracne_project_2022-03-21_out_.final/consensus_network_ncol_.txt',header = T,check.names = F)
tf$source=stringr::str_remove(tf$source,pattern = 'MSC_')
tf$target=stringr::str_remove(tf$target,pattern = 'MSC_')
write.table(tf, file ='./test/network/SJAR/project_2022-03-21/output_tf_sjaracne_project_2022-03-21_out_.final/consensus_network_ncol_.txt',col.names = T,row.names = F,quote = F,sep = '\t')

sig=read.table('./test/network/SJAR/project_2022-03-21/output_sig_sjaracne_project_2022-03-21_out_.final/consensus_network_ncol_.txt',header = T,check.names = F)
sig$source=stringr::str_remove(sig$source,pattern = 'MSC_')
sig$target=stringr::str_remove(sig$target,pattern = 'MSC_')
write.table(sig, file ='./test/network/SJAR/project_2022-03-21/output_sig_sjaracne_project_2022-03-21_out_.final/consensus_network_ncol_.txt',col.names = T,row.names = F,quote = F,sep = '\t')

# 创立文件夹
# 关键步骤
analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_name,
                                            network_dir=network.dir, network_project_name=network.project.name)

############### Step 1: 载入表达矩阵 (exp-load, exp-cluster, exp-QC) ###############
# 全部用自己的生成的文件做，不依赖任何内置数据

# 我们已经制作了自己的表达矩阵（在第一节project_2022-03-15里面），此处就用自己的
# 我们自己的exp-QC在\test\project_2022-03-15\DATA下



# 当然限于时间你跟我不一样，所以你下面要换你的时间




load('./test/project_2022-03-15/DATA/network.par.Step.exp-QC.RData') # RData saved after QC in the network construction step
# 把exp-QC给到analysis.par中
analysis.par$cal.eset <- network.par$net.eset

# 保存勿忘！在\test\driver_时间\DATA中
NetBID.saveRData(analysis.par=analysis.par,step='exp-QC')

# 注意留心观察，我们没有依赖任何内置数据
draw.combineDE(DA_driver_comb,pdf_file=sprintf('%s/combineDA.pdf',analysis.par$out.dir.PLOT)) # Save it as PDF
