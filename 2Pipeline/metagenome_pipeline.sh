# 一、宏基因组有参分析流程 metagenome reference-based pipeline (humann2)

# 系统要求 System: Linux Ubuntu 18.04 / CentOS 7.5
# 依赖软件 Sofware: KneadData、metaphlan2、humann2
# 运行前准备
# 1. 按7software, 8database目录中软件和数据库按课件说明安装，并添加环境变量
# 2. 学员U盘复制测序数据3metagenome目录到C:或服务器~目录
# 3. Rstudio打开pipeline_ref_humann2.sh文件，Terminal中切换至工作目录

# 中文教程：https://mp.weixin.qq.com/s/XkfT5MAo96KgyyVaN_Fl7g
# 英文教程：https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-Tutorial-(Humann2)

# 上课演示Linux服务器：数据己经拷备到服务器/db/3metagenome目录，
# 用户登陆Rstudio网页版：192.168.1.107:8787 即可


# 1. 了解工作目录和文件

# 建立并进入工作目录
mkdir -p 3metagenome
cd 3metagenome

# 目录
mkdir -p doc # 实验设计 design file
ln -s /db/3metagenome/doc/design.txt doc/
head -n3 doc/design.txt
mkdir -p seq # 原始数据 raw data
ln -s /db/3metagenome/seq/*.fq seq/
mkdir -p temp # 临时文件 temp directory for intermediate files

# 文件说明

# pipeline_ref_humann2.sh 有参分析主流程

# seq/*.fq 原始测序数据，公司返回的测序结果，通常为一个样品一对fastq格式文件
# 如果测序数据是.gz结尾的压缩文件，使用gunzip解压，注意测序文件命名，结果可以是fastq/fq，左端可以_1/R1.fq
# gunzip seq/* # 如果压缩文件还需要解压
# 测序数据有12个样本，共24个文件。只取1%抽样用于测试。
head -n4 seq/p136C_1.fq


# 2. 序列质控和去宿主 Qaulity control and remove host contamination

# kneaddata是一套工作流，依赖trimmomatic进行质控和去接头；依赖bowtie2比对宿主，并筛选非宿主序列
# kneaddata -h # 显示帮助
# kneaddata_database # 查看可用数据库
# kneaddata_database --download human_genome bowtie2 ./ # 如下载人类基因组bowtie2索引至当前目录，并脚本

# 以p136C单样品合并为例
kneaddata -i seq/p144C_1.fq -i seq/p144C_2.fq \
  -o temp/qc -v -t 4 --remove-intermediate-output \
  --trimmomatic /conda2/share/trimmomatic-0.36-3/ --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" \
  --bowtie2-options "--very-sensitive --dovetail" -db /db/bowtie2/Homo_sapiens   
    
# 现实中是有一大堆样品，你可以逐个修改样品名运行，如果服务器性能允许可以并行加速分析
# parallel --citation # 打will cite以后不再提醒
parallel -j 3 --xapply \
  'kneaddata -i {1} -i {2} \
  -o temp/qc -v -t 3 --remove-intermediate-output \
  --trimmomatic /conda2/share/trimmomatic-0.36-3/ --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" \
  --bowtie2-options "--very-sensitive --dovetail" -db /db/bowtie2/Homo_sapiens' \
 ::: seq/*_1.fq ::: seq/*_2.fq

# 质控结果统计
kneaddata_read_count_table --input temp/qc --output seq/kneaddata_read_counts.txt
cat seq/kneaddata_read_counts.txt

# 合并质控后样品文件：有参宏基因组不考虑双端，将单双端和双端不完全序列共4个文件合并
mkdir -p temp/concat
for i in `tail -n+2 doc/design.txt | cut -f 1`;do \
  cat temp/qc/${i}*data_paired* > temp/concat/${i}.fq; done


# 3. 计算功能和代谢通路

humann2_config # 查看参数和数据库位置是否正确
# metaphlan2数据库默认位于程序所在目录的db_v20和databases下各一份
humann2 --input temp/concat/p136C.fq  \
  --output temp/ \
  --threads 8


# 并行处理所有样品
time parallel -j 12 \
  'humann2 --input {}  \
  --output temp/ ' \
  ::: temp/concat/*.fq
  

# 4. 物种组成分析

# 样品结果合并
merge_metaphlan_tables.py temp/*_humann2_temp/*_metaphlan_bugs_list.tsv | sed 's/_metaphlan_bugs_list//g' > result/metaphlan2.tsv # ;/^\#/d

# 测序数据量有限，结果不够丰度，请用原始物种代表简化表下游分析更精彩
cp full/metaphlan2.tsv result/metaphlan2.tsv

# 转换为spf格式方便stamp分析
metaphlan_to_stamp.pl result/metaphlan2.tsv > result/metaphlan2.spf

# 下载design.txt和metaphlan2.spf使用stamp分析


# 5. 物种组成分析和可视化进阶

# 绘制热图
metaphlan_hclust_heatmap.py --in result/metaphlan2.tsv --out result/metaphlan2.pdf #  -c bbcry --top 25 --minv 0.1 -s log 
# c设置颜色方案，top设置物种数量，minv最大相对丰度，s标准化方法，文件名结尾可先pdf/png/svg三种图片格式。更多说明详见 metaphlan_hclust_heatmap.py -h

# GraPhlAn图
# metaphlan2 to graphlan
export2graphlan.py --skip_rows 1,2 -i result/metaphlan2.tsv \
  --tree temp/merged_abundance.tree.txt \
  --annotation temp/merged_abundance.annot.txt \
  --most_abundant 100 --abundance_threshold 1 --least_biomarkers 10 \
  --annotations 5,6 --external_annotations 7 --min_clade_size 1
# graphlan annotation
graphlan_annotate.py --annot temp/merged_abundance.annot.txt temp/merged_abundance.tree.txt \
  temp/merged_abundance.xml
# output PDF figure, annoat and legend
graphlan.py temp/merged_abundance.xml result/metaphaln2_graphlan.pdf --external_legends --dpi 300 

# LEfSe差异分析和Cladogram
# 修改样本品为组名
sed '1 s/p[0-9]*//g' result/metaphlan2.tsv | grep -v '#' > metaphlan2/lefse.txt
# 格式转换为lefse内部格式
lefse-format_input.py metaphlan2/lefse.txt temp/input.in -c 1 -o 1000000
# 运行lefse
run_lefse.py temp/input.in temp/input.res
# 绘制物种树注释差异
lefse-plot_cladogram.py temp/input.res metaphlan2/lefse_cladogram.pdf --format pdf --dpi 600 
# 绘制所有差异features柱状图
lefse-plot_res.py temp/input.res metaphlan2/lefse_res.pdf --format pdf --dpi 600
# 绘制单个features柱状图(同STAMP中barplot)
sort -k3,3n temp/input.res |less -S # 查看差异features列表
lefse-plot_features.py -f one --feature_name "k__Bacteria.p__Firmicutes" --format pdf \
  temp/input.in temp/input.res metaphlan2/Firmicutes.pdf 
# 批量绘制所有差异features柱状图
lefse-plot_features.py -f diff --archive none --format pdf \
  temp/input.in temp/input.res metaphlan2/features


# 5. 功能组成分析

# 合并所有样品，通路包括各功能和具体的物种组成，还有基因家族(太多)，通路覆盖度层面可以进分析
humann2_join_tables --input temp/ --file_name pathabundance --output result/pathabundance.tsv
sed -i 's/_Abundance//g' result/pathabundance.tsv

# 标准化为相对丰度relab或百万分数cpm，原始为什么值？用R统计
humann2_renorm_table --input result/pathabundance.tsv --units relab --output result/pathabundance_relab.tsv

# 分层结果
humann2_split_stratified_table --input result/pathabundance_relab.tsv --output result/
# 结果stratified(每个菌的功能组成)和unstratified(功能组成)两个


# 筛选某个通路结果用STAMP展示
head -n 1 result/pathabundance_relab_stratified.tsv | \
  sed 's/# Pathway/MetaCyc_pathway/' \
  > result/pathabundance_relab_stratified_LACTOSECAT-PWY.spf
grep "LACTOSECAT-PWY" result/pathabundance_relab_unstratified.tsv \
  >> result/pathabundance_relab_stratified_LACTOSECAT-PWY.spf









# 合并所有样品至同一文件
cat temp/*.merged.fq > temp/all.fq
ls -l temp/all.fq # 639 Mb
# 删除中间文件
rm temp/*.merged.fq
# 压缩原始文件节省空间
# gzip seq/*

# 3. 切除引物与质控 Cut primers and quality filter
# Cut barcode 10bp + V5 19bp in left and V7 18bp in right
vsearch --fastx_filter temp/all.fq \
  --fastq_stripleft 29 --fastq_stripright 18 \
  --fastqout temp/stripped.fq
# 10S, 39Mb

# 质量控制fastq filter, keep reads error rates less than 1%
vsearch --fastx_filter temp/stripped.fq \
  --fastq_maxee_rate 0.01 \
  --fastaout temp/filtered.fa
# 16S, 6Mb


# 4. 有参比对，如Greengenes，可用于picurst, bugbase分析

# 基于greengenes的结果全部保存于gg目录
mkdir -p gg/

# 生成OTU表 Create OTUs table
vsearch --usearch_global temp/filtered.fa \
  --db db/97_otus.fasta \
  --id 0.97 \
  --otutabout gg/otutab.txt --threads 4 

# 与GG所有97% OTUs比对，i5笔记本4线程用时17:42，内存729Mb
usearch -otutab_stats gg/otutab.txt \
	-output gg/otutab.stat 
cat gg/otutab.stat 

# 接下来的主成分、物种、组间比较，可用STAMP，或基于无参的R脚本方法


# 5. 生成STAMP输入物种文件

mkdir -p gg/stamp

# OTU表

# 筛选高丰度OTU，关注重点，减少计算量，降低FDR校正效应，推荐HiSeq选万分之五，可选万一，千一，千五
usearch -otutab_trim gg/otutab.txt -min_otu_freq 0.00005 -output temp/otutab_trim.txt
usearch -otutab_stats temp/otutab_trim.txt \
	-output temp/otutab_trim.stat 
cat temp/otutab_trim.stat

# OTU表：删除#和空格，序列前添加字母
sed '1 s/#OTU //;s/^/S/' temp/otutab_trim.txt > gg/otutab.spf

# 获取OTUs
cut -f 1 gg/otutab.txt | tail -n+2 > temp/otu.id
usearch -fastx_getseqs db/97_otus.fasta -labels temp/otu.id -fastaout gg/otu.fa

# 物种注释
# usearch -sintax gg/otu.fa -db db/gg_16s_13.5.fa \
#   -strand both -tabbedout gg/sintax.txt -sintax_cutoff 0.6
# Memory limit of 32-bit process exceeded, 64-bit build required
usearch -sintax gg/otu.fa -db db/rdp_16s_v16_sp.fa \
  -strand both -tabbedout gg/sintax.txt -sintax_cutoff 0.6
# 各分类级
for i in p c o f g;do
  usearch -sintax_summary gg/sintax.txt \
  -otutabin gg/otutab.txt \
  -rank ${i} \
  -output gg/tax_sum_${i}.txt
done  
# 可直接stamp统计分析，或用R绘制堆叠柱状图


# 6. OTU表转换为KO表

biom convert -i gg/otutab.txt -o gg/otutab.biom --table-type="OTU table" --to-json
biom summarize-table -i gg/otutab.biom
# 校正拷贝数
normalize_by_copy_number.py -i gg/otutab.biom -o gg/otutab_norm.biom -c /db/picrust/16S_13_5_precalculated.tab.gz
# 预测宏基因组KO表
predict_metagenomes.py -i gg/otutab_norm.biom -o gg/ko.biom -c /db/picrust/ko_13_5_precalculated.tab.gz
predict_metagenomes.py -f -i gg/otutab_norm.biom -o gg/ko.txt  -c /db/picrust/ko_13_5_precalculated.tab.gz
# 按功能级别分类汇总, -c指输出类型，有KEGG_Pathways, COG_Category, RFAM三种，-l是级别，分4级，初始KO为4级，可全并为1-3级
categorize_by_function.py -f -i gg/ko.biom -c KEGG_Pathways -l 3 -o gg/ko3.txt
categorize_by_function.py -f -i gg/ko.biom -c KEGG_Pathways -l 2 -o gg/ko2.txt
categorize_by_function.py -f -i gg/ko.biom -c KEGG_Pathways -l 1 -o gg/ko1.txt

# 生成stamp可分析表
sed  -i '/# Constru/d;s/#OTU //' gg/ko*.txt # 删除表头多作注释
# 将最后一列注释调整为第一列即为stamp使用的spf格式，大家也可用excel手动调整
num=`head -n1 gg/ko3.txt|wc -w`
paste <(cut -f $num gg/ko.txt) <(cut -f 1-$[num-1] gg/ko.txt) > gg/ko.spf #  | cut -f 1 -d ';' 删除低级冗余信息stamp发现有unclassified无法打开
paste <(cut -f $num gg/ko3.txt) <(cut -f 1-$[num-1] gg/ko3.txt) > gg/ko3.spf #  | cut -f 1 -d ';' 删除低级冗余信息stamp发现有unclassified无法打开
paste <(cut -f $num gg/ko3.txt) <(cut -f 1-$[num-1] gg/ko3.txt) > gg/ko3.spf #  | cut -f 1 -d ';' 删除低级冗余信息stamp发现有unclassified无法打开
paste <(cut -f $num gg/ko1.txt) <(cut -f 1-$[num-1] gg/ko3.txt) > gg/ko1.spf #  | cut -f 1 -d ';' 删除低级冗余信息stamp发现有unclassified无法打开








# 二、扩增子无参分析流程 16S Amplicon de novo pipeline


# 4. (可选Denovo方法)去冗余与生成OTUs Dereplication and cluster otus 

# 4.1 序列去冗余，推荐使用vsearch，并添加miniuniqusize为8，去除低丰度，增加计算速度

# 去冗余Find unique read sequences and abundances
vsearch --derep_fulllength temp/filtered.fa \
  --sizeout --minuniquesize 8 \
  --output temp/uniques.fa
# 5S, 368Mb

# 4.2 生成OTU

# 有两种方法选择，如最新的unoise3推荐，也可以使用传统的97%聚类，供备选。

# 可选97%聚类，不推荐，除非reviewer要求
# usearch -cluster_otus temp/uniques.fa \
#   -otus temp/otus.fa \
#   -relabel OTU_
# 4S, 16Mb, 984 OTUs, 428 chimeras

# 预测生物学序列OTU并去除嵌合 Denoise: predict biological sequences and filter chimeras
usearch -unoise3 temp/uniques.fa \
  -zotus temp/zotus.fa
# 43S, 56Mb, 3326 OTUs, 348 chimeras

# 修改序列名：格式调整 format OTU prefix
sed 's/Zotu/OTU_/g' temp/zotus.fa > temp/otus.fa


# 4.3 基于参考去嵌合
# 准备数据库 http://www.drive5.com/sintax
# cd db
# wget http://www.drive5.com/sintax/silva_16s_v123.fa.gz
# gunzip silva_16s_v123.fa.gz
# cd ..
vsearch --uchime_ref temp/otus.fa \
  --db db/silva_16s_v123.fa \
  --nonchimeras result/otus.fa


# 4.4 生成OTU表 Creat OTUs table
vsearch --usearch_global temp/filtered.fa \
  --db result/otus.fa \
  --id 0.97 --threads 4 \
  --otutabout result/otutab.txt
# Matching query sequences: 639169 of 740493 (86.32%)

# 5. OTU表统计和标准化

# OTU表简单统计 Summary OTUs table
usearch -otutab_stats result/otutab.txt \
	-output result/otutab_report.txt 
cat result/otutab_report.txt 

# 等量抽样标准化 normlize by subsample to 10000
# 我们看到最小样品数据量为3.1万，可以抽样至3万
usearch -otutab_norm result/otutab.txt \
	-sample_size 31541 \
	-output result/otutab_norm.txt 
usearch -otutab_stats result/otutab_norm.txt \
	-output result/otutab_norm_report.txt 
cat result/otutab_norm_report.txt 



# 6. 物种注释 Assign taxonomy
# wget http://www.drive5.com/sintax/rdp_16s_v16_sp.fa.gz
# gunzip rdp_16s_v16_sp.fa.gz
usearch -sintax result/otus.fa -db db/rdp_16s_v16.fa \
  -strand both -tabbedout result/sintax.txt -sintax_cutoff 0.8
# 21s, 103Mb
# usearch -sintax result/otus.fa -db db/rdp_16s_v16.fa \
#   -strand plus -tabbedout result/sintax.txt -sintax_cutoff 0.8
# 单链下18s， 103Mb

head -n2 result/sintax.txt
# 统计门纲目科属，使用 rank参数 p c o f g，为phylum, class, order, family, genus缩写
# 我们先计算门phylum, 目order和属genus水平汇总表方便观察
mkdir -p tax
for i in p c o f g;do
  usearch -sintax_summary result/sintax.txt \
  -otutabin result/otutab_norm.txt \
  -rank ${i} \
  -output tax/sum_${i}.txt
done
# usearch -sintax_summary result/sintax.txt -otutabin result/otutab_norm.txt \
#   -rank p -output tax/tax_phylum.txt
# usearch -sintax_summary result/sintax.txt -otutabin result/otutab_norm.txt \
#   -rank o -output tax/tax_order.txt
# usearch -sintax_summary result/sintax.txt -otutabin result/otutab_norm.txt \
#   -rank g -output tax/tax_genus.txt
  
# Taxonomy中异常字符
sed -i 's/(//g;s/)//g;s/\"//g;s/\/Chloroplast//g' tax/sum_*.txt

# 格式化物种注释：去除sintax中置信值，只保留物种注释，替换:为_，删除引号
cut -f 1,4 result/sintax.txt|sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' > result/otus.tax
head -n2 result/otus.tax
# 注意注释是非整齐的，由于新物种只是相近而不同
# 生成物种表格：注意OTU中会有末知为空白，补齐分类未知新物种为Unassigned
awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";\
  split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
  print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"];}' \
  result/otus.tax >temp/otus.tax
sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-7| sed '1 s/^/#OTUID\tKindom\tPhylum\tClass\tOrder\tFamily\tGenus\n/' > tax/taxtab.txt
head -n3 tax/taxtab.txt



# 7. Alpha多样性 Alpha diversity
# Calculate all alpha diversity, details in http://www.drive5.com/usearch/manual/alpha_metrics.html
mkdir -p alpha
usearch -alpha_div result/otutab_norm.txt \
  -output alpha/alpha.txt 
# 稀释曲线：取1%-100%的序列中OTUs数量 Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method fast / with_replacement / without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
usearch -alpha_div_rare result/otutab_norm.txt \
  -output alpha/alpha_rare.txt -method without_replacement
# 9s, 6.8Mb



# 8. Beta多样性 Beta diversity
# 基于OTU构建进化树 Make OTU tree
usearch -cluster_agg result/otus.fa -treeout result/otus.tree
# 25s, 81Mb
# 结果有多个文件，需要目录
mkdir -p beta/ 
# 生成5种距离矩阵：bray_curtis, euclidean, jaccard, manhatten, unifrac
usearch -beta_div result/otutab_norm.txt -tree result/otus.tree \
  -filename_prefix beta/





# 二、统计绘图

# 绘图所需R脚本位于script目录中

# Windows中需要Rstudio打开脚本，设置工作目录(data)，全选运行即可，
# 参数修改位于"解析命令行"部分，可修改每个参数的default值改变输入、输入文件等
# Linux用需目前命令行操作设置参数非常方便


## 1. Alpha多样性指数箱线图 Alpha index in boxplot
# 脚本放在环境变量可直接使用文件名，否则指定目录
# 显示脚本帮助
Rscript ./script/alpha_boxplot.r -h
# 默认画richness
Rscript ./script/alpha_boxplot.r
# 完整参数，输出文件名默认为alpha指数类型
Rscript ./script/alpha_boxplot.r -i alpha/alpha.txt -t richness \
  -d doc/design.txt -n group \
  -o alpha/richness \
  -w 4 -e 2.5
# 绘制chao1
Rscript ./script/alpha_boxplot.r -t chao1 



## 2. Alpha稀释曲线 Rarefraction curve
# 显示帮助
Rscript ./script/alpha_rare.r -h
# 默认绘制4x2.5英寸的样品和组均值图
Rscript ./script/alpha_rare.r
# 指定输入文件和实验设计，实验组列，图片长宽和输出文件前缀
Rscript ./script/alpha_rare.r -i alpha/alpha_rare.txt \
  -d doc/design.txt -n group \
  -o alpha/rare_ \
  -w 4 -e 2.5 
# 只想用两组展示PCoA，删除B组，使用新的design2.txt(手动删除了B组行)
Rscript ./script/alpha_rare.r -d doc/design2.txt -n group \
  -o alpha/rare_2
  


## 3. Beta主坐标轴分析 PCoA

# 展示样品间距离分布，统计组间是否显著，也用于异常样品筛选

# 显示脚本帮助
Rscript ./script/beta_pcoa.r -h
# 默认基于bray_curtis距离
Rscript ./script/beta_pcoa.r
# 完整默认参数
Rscript ./script/beta_pcoa.r -i beta/bray_curtis.txt -t bray_curtis \
  -d doc/design.txt -n group \
  -o beta/pcoa_bray_curtis \
  -w 4 -e 2.5 
# 基于unifrac距离
Rscript ./script/beta_pcoa.r -t unifrac
# 只想用两组展示PCoA，删除B组，使用新的design2.txt(手动删除了B组行)
Rscript ./script/beta_pcoa.r -i beta/bray_curtis.txt \
  -d doc/design2.txt -n group \
  -o beta/pcoa_bray_curtis_AC



## 4. 限制性主坐标轴分析 Constrained PCoA

# 展示组间最大差异

# 显示脚本帮助
Rscript ./script/beta_cpcoa.r -h

# 基于bray距离计算CCA，默认bray方法
Rscript ./script/beta_cpcoa.r

# 基于jaccard距离计算CCA
Rscript ./script/beta_cpcoa.r -t jaccard

# 附完整参数
Rscript ./script/beta_cpcoa.r -i result/otutab_norm.txt -t bray \
  -d doc/design.txt -n group \
  -o beta/cpcoa_bray



## 5. 物种丰度柱状图 Taxonomy barplot

# 显示帮助，主要是参数说明
Rscript ./script/tax_stackplot.r -h

# 默认按phylum和前8类展示, 4X2.5
Rscript ./script/tax_stackplot.r 
# Legend too long, main text overlap. Increase figure size.

# 按目前10，图片宽6 x 4
Rscript ./script/tax_stackplot.r -t order \
  -b 10 -w 6 -e 4




# 第二天，第1讲

## 6. 组间差异比较

# 需要otu表、实验设计和物种注释

# 结果目录
mkdir -p compare

# 显示帮助
Rscript ./script/compare_edgeR.r -h # 显示帮助

# 默认参数：计算group分类下A-B比较
Rscript ./script/compare_edgeR.r

# 计算A-C
Rscript ./script/compare_edgeR.r -c A-C

# 按genotype分组下KO-WT
Rscript ./script/compare_edgeR.r -n genotype -c KO-WT

## 7. 绘制火山图、热图和曼哈顿图，同上

# 数据矩阵在edgeR_KO-WT_sig.txt文件中
# 样品注释在design.txt中，用于列分组注释
# OTU物种注释来自taxtab.txt文件(可选)# Vsearch 16S Amplicon pipeline

# Usearch32位版分析>4GB文件受限，64位收费。Vsearch免费且不受任何限制

# 命令计算时间基于18个5万条序列样品，于win10 2.3G i5双核4线程笔记本

# 1. 了解工作目录和文件

# 建立临时和结果目录
mkdir -p temp # 临时文件 temp directory for intermediate files
mkdir -p result # 最终结果 important results

# 文件
# vsearch.sh 分析主流程
# db/*.fa # 参考数据库 database files
# seq/*.fq.gz 压缩的原始测序数据
# doc/design.txt 实验设计文件



# 2. 合并双端序列与样品拆分 Merge paired reads and label samples

# 测序数据解压
# gunzip seq/*

# 依照实验设计批处理并合并
for i in `tail -n+2 doc/design.txt | cut -f 1`;do
  vsearch --fastq_mergepairs seq/${i}_1.fq --reverse seq/${i}_2.fq \
  --fastqout temp/${i}.merged.fq --relabel ${i}.
done 

# 合并所有样品至同一文件
cat temp/*.merged.fq > temp/all.fq
ls -l temp/all.fq



# 3. 切除引物与质控 Cut primers and quality filter
# 请按实际修改，如Cut barcode 10bp + V5 19bp in left and V7 18bp in right
time vsearch --fastx_filter temp/all.fq \
  --fastq_stripleft 29 --fastq_stripright 18 \
  --fastqout temp/stripped.fq # 1m5s
# 质量控制fastq filter, keep reads error rates less than 1%
time vsearch --fastx_filter temp/stripped.fq \
  --fastq_maxee_rate 0.01 \
  --fastaout temp/filtered.fa # 55s
#761431 sequences kept (of which 0 truncated), 5627 sequences discarded.



# 4. 去冗余与生成OTUs Dereplication and cluster otus
# 4.1 序列去冗余，推荐使用vsearch，并添加miniuniqusize为8，去除低丰度，增加计算速度
time vsearch --derep_fulllength temp/filtered.fa \
  --sizeout --minuniquesize 8 \
  --output temp/uniques.fa # 5s


## 此处我们用基于reference的去嵌合，下载rdp_gold.fa作
#为reference数据库
#wget http://drive5.com/uchime/rdp_gold.fa

# 聚类方式生成OTU
time vsearch --cluster_fast temp/uniques.fa \
  --id 0.97 --centroids temp/otus.fa \
  --relabel OTU_ # 3s Clusters: 1244 --uc temp/clusters.uc
dos2unix temp/otus.fa # 去除windows换行符


# 细菌可用Usearch作者整理的RDP Gold数据库去除嵌合体
# wget http://drive5.com/uchime/rdp_gold.fa
time vsearch --uchime_ref temp/otus.fa \
  --db db/rdp_gold.fa \
  --nonchimeras result/otus.fa # 9s, 1041 non-chimeras,
dos2unix result/otus.fa # 去除windows换行符
	
# Create OTUs table
time vsearch --usearch_global temp/filtered.fa \
  --db result/otus.fa \
  --id 0.97 \
  --otutabout result/otutab.txt --threads 4 # 8m54s
dos2unix result/otutab.txt
# 检查是否还有windows换行符
cat -A result/otutab.txt | head

# 物种注释
# vsearch --usearch_global result/otus.fa --db db/rdp_16s_v16.fa --biomout out_tax.txt --id 0.97
