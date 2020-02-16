# 一、扩增子有参分析流程 16S Amplicon Pipeline -- Reference-based 

# 系统要求 System: Windows 10 / Linux Ubuntu 18.04
# 依赖软件 Sofware: usearch, vsearch2.8, R3.4 and Rstudio1.1, gitforwidnows(仅win10需要)
# 预测宏基因组部分需要 picurst + greengene13.5数据库

# 运行前准备
# 1. 按7software目录中课件说明安装依赖软件并添加环境变量
# 2. 学员U盘复制测序数据2amplicon目录到C:或服务器~目录
# 3. Rstudio打开2amplicon/pipeline_reference.sh文件，Terminal中切换至工作目录

## Windows用户：切换至工作目录
# cd /c/2amplicon

# Linux服务器
# 用户登陆Rstudio网页版：内网 192.168.1.107:8787，外网 210.75.224.32:8787 进行练习
# mkdir 2amplicon # 建立项目目录
# cd 2amplicon/ # 进入工作目录
# ln -s /db/2amplicon/pipeline_reference.sh ./ # 链接流程至工作目录，保持更新
# 右侧打开流程文件，链接无权限修改，保存时会弹出另存即可

# 1. 了解工作目录和文件

# 目录：文件归类思想
# 进入扩增子分析的工作目录 2amplicon
mkdir -p seq # 原始数据 raw data
ln -s /db/2amplicon/seq/* ./seq
mkdir -p doc # 实验设计及相关文档 design / metadata
ln -s /db/2amplicon/doc/* ./doc
mkdir -p temp # 临时文件 temp directory for intermediate files
mkdir -p gg # 以greengene13.5为参考数据库的结果目录

# 文件说明

# pipeline*.sh 分析主流程

# seq/*.fq 原始测序数据，公司返回的测序结果，通常为一个样品一对fastq格式文件
# 如果测序数据是.gz结尾的压缩文件，使用gunzip解压
# gunzip seq/*
head -n4 seq/KO1_1.fq

# doc/design.txt 实验设计文件
head -n3 doc/design.txt

# /db/gg/*.fa  Greengene16S数据库
# greengene 13_8: ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
head -n2 /db/gg/97_otus.fasta
head -n2 /db/gg/97_otu_taxonomy.txt

# 2. 合并双端序列与样品标签 Merge paired reads and label samples

# 以WT1单样品合并为例
vsearch -fastq_mergepairs seq/WT1_1.fq -reverse seq/WT1_2.fq \
	-fastqout temp/WT1.merged.fq \
	-relabel WT1. # 1S, 63Mb

# 依照实验设计批处理并合并
# rstudio中运行会异常中断，中途现Ctrl+C，是其它软件如词典引起
for i in `tail -n+2 doc/design.txt | cut -f 1`;do
  vsearch --fastq_mergepairs seq/${i}_1.fq --reverse seq/${i}_2.fq \
  --fastqout temp/${i}.merged.fq --relabel ${i}. &
done
# 合并所有样品至同一文件
cat temp/*.merged.fq > temp/all.fq
ls -l temp/all.fq # 660 Mb
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


# 4. 有参比对，如Greengenes，可用于picurst, stamp

# 生成OTU表 Create OTUs table
vsearch --usearch_global temp/filtered.fa \
  --db /db/gg/97_otus.fasta \
  --id 0.97 \
  --otutabout gg/otutab.txt --threads 9

# 与GG所有97% OTUs比对，i5笔记本4线程用时17:42，内存729Mb
usearch -otutab_stats gg/otutab.txt \
	-output gg/otutab.stat 
cat gg/otutab.stat 

# 接下来的主成分、物种、组间比较，可用STAMP，或基于无参的R脚本方法


# 5. 整理结果：OTU表、代表序列和物种注释，生成STAMP、lefse输入文件

# 简化分组信息，同目录方便使用
cut -f 1,5-7 doc/design.txt > gg/metadata.txt 

# 筛选高丰度OTU，关注重点，减少计算量，降低FDR校正效应，获取人类可读量的结果，推荐HiSeq选万1，可选千1，十万1
usearch -otutab_trim gg/otutab.txt -min_otu_freq 0.0001 -output temp/otutab_trim.txt
usearch -otutab_stats temp/otutab_trim.txt -output temp/otutab_trim.stat 
cat temp/otutab_trim.stat
# 566487 / 615803 = 91.2% reads, 903 / 6740 = 13.4% 的高丰度OTUs

# OTU表：删除#和空格，序列前添加字母，纯数字行/列名称软件和R语言容易报错
sed '1 s/#OTU //;s/^/S/' temp/otutab_trim.txt > gg/otutab.spf

# 获取OTUs和物种注释
cut -f 1 gg/otutab.txt | tail -n+2 > temp/otu.id
usearch -fastx_getseqs /db/gg/97_otus.fasta -labels temp/otu.id -fastaout gg/otu.fa
awk 'BEGIN{OFS=FS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $1,a[$1]}' \
  /db/gg/97_otu_taxonomy.txt temp/otu.id > gg/otu.tax
sed 's/; /\t/g' gg/otu.tax | sed '1 i ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' > gg/otu.tax8

# Linux下计算各分类级汇总表，用于stamp和lefse统计
Rscript ../6script/taxonomy_summary.R
# Windows下需要使用Rstudio打开脚本，设置工作目录，全选运行
wc -l gg/tax_sum* # 各级数量
# 可直接stamp统计分析，或用R绘制堆叠柱状图
# stamp拉开otu表和metadata，默认多组比较查看pca, heatmap, barplot, boxplot，导出差异; 两组比较，Extended error bar, bar/box


# 6. R脚本差异分析及可视化

## 6.1 Alpha多样性 Alpha diversity
# Calculate all alpha diversity, details in http://www.drive5.com/usearch/manual/alpha_metrics.html
mkdir -p gg/alpha
# 样品抽平至最小值，才能评估alpha多样性
usearch10 -otutab_norm gg/otutab.txt \
	-sample_size 30000 \
	-output gg/otutab_norm.txt 
# 计算15种多样性指数
usearch -alpha_div gg/otutab_norm.txt \
  -output gg/alpha/index.txt 
# 稀释曲线：取1%-100%的序列中OTUs数量 Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method fast / with_replacement / without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
usearch -alpha_div_rare gg/otutab_norm.txt \
  -output gg/alpha/rare.txt -method without_replacement # 9s, 6.8Mb

# 显示脚本帮助
Rscript ../6script/alpha_boxplot.R -h

# 默认画richness
Rscript ../6script/alpha_boxplot.R

# 绘制箱线图+差异统计，方法可选richness, chao1，shannon_e等，详见head -n1 gg/alpha/index.txt
Rscript ../6script/alpha_boxplot.R -i gg/alpha/index.txt \
  -d doc/design.txt -n group -t richness \
  -o gg/alpha/richness \
  -w 4 -e 2.5 

# 绘制chao1
Rscript ../6script/alpha_boxplot.R -t chao1 
  
# 绘制稀释曲线
Rscript ../6script/alpha_rarefaction_curve.R -i gg/alpha/rare.txt \
  -d doc/design.txt -n group \
  -o gg/alpha/rare_ \
  -w 4 -e 2.5 
# Windows中使用Rstudio打开，设置项目工作目录，并修改输入输出文件名，运行即可


# 6.2 Beta多样性 Beta diversity
# 结果有多个文件，需要目录
mkdir -p gg/beta/ 
# 基于OTU构建进化树 Make OTU tree
usearch -cluster_agg gg/otu.fa -treeout gg/beta/otu.tree
# 全长建树比较慢，10线程 06:29 613Mb
# 生成5种距离矩阵：bray_curtis, euclidean, jaccard, manhatten, unifrac
usearch -beta_div gg/otutab.txt -tree gg/beta/otu.tree \
  -filename_prefix gg/beta/

# 完整默认参数
Rscript ../6script/beta_pcoa.R -i gg/beta/bray_curtis.txt -t bray_curtis \
  -d doc/design.txt -n group \
  -o gg/beta/pcoa_bray_curtis \
  -w 4 -e 2.5 
# 基于unifrac距离
Rscript ../6script/beta_pcoa.R -t unifrac

# 限制性PCoA / CCA，找分组间最大差异，距离t类型,可选manhattan, euclidean, canberra, bray, kulczynski, jaccard, gower, altGower, morisita, horn, mountford, raup , binomial, chao, cao or mahalanobis
Rscript ../6script/beta_cpcoa.R -i gg/otutab_norm.txt -t bray \
  -d doc/design.txt -n group \
  -o gg/beta/cpcoa_bray
# 来一个欧式距离试试
Rscript ../6script/beta_cpcoa.R -t euclidean

## 6.3 物种丰度柱状图 Taxonomy barplot

mkdir -p gg/tax/ 

# 显示帮助，主要是参数说明
Rscript ../6script/tax_stackplot.R -h

# 默认按phylum和前8类展示, 4X2.5
Rscript ../6script/tax_stackplot.R
# Legend too long, main text overlap. Increase figure size.

# 按目前10，图片宽6 x 4，按目水平绘图显著高丰度前10类
Rscript ../6script/tax_stackplot.R -t o \
  -b 10 -w 10 -e 5


## 6.4 差异比较

mkdir -p gg/compare

# 显示帮助
Rscript ../6script/compare.R -h # 显示帮助

# 计算A-C
Rscript ../6script/compare.R -c A-C

# 按genotype分组下KO-WT
Rscript ../6script/compare.R -n genotype -c KO-WT

# 目前只有OTU表的count值，需生成界、门、纲、目、科、属、种级别
Rscript ../6script/taxonomy_summary_count.R



# 7. LEfSe差异分析和Cladogram

mkdir -p gg/lefse

# 生成lefse输入，output_lefse.txt
Rscript ../6script/taxonomy_summary.R -i gg/otutab.txt \
  -t gg/otu.tax8 -T 0.5 -o gg/lefse/sum
wc -l gg/lefse/sum*

# 格式转换为lefse内部格式
lefse-format_input.py gg/tax_sum_lefse.txt temp/input.in -c 1 -o 1000000
# 运行lefse
run_lefse.py temp/input.in temp/input.res
# 绘制物种树注释差异
lefse-plot_cladogram.py temp/input.res gg/lefse/cladogram.pdf --format pdf
# 绘制所有差异features柱状图
lefse-plot_res.py temp/input.res gg/lefse/res.pdf --format pdf
# 绘制单个features柱状图(同STAMP中barplot)
head temp/input.res # 查看差异features列表
lefse-plot_features.py -f one --feature_name "k__Bacteria.p__Proteobacteria.c__Alphaproteobacteria.o__Rhizobiales" --format pdf \
  temp/input.in temp/input.res gg/lefse/Rhizobiales.pdf 
# 批量绘制所有差异features柱状图，慎用(几百张差异结果柱状图阅读也很困难)
lefse-plot_features.py -f diff --archive none --format pdf \
  temp/input.in temp/input.res gg/lefse/



# 8. PICRUSt分析：OTU表转换为KO表

# download_picrust_files.py 下载数据数据至 /usr/local/lib/python2.7/dist-packages/picrust/data
# 转换为OTU表通用格式，方便下游分析和统计
biom convert -i gg/otutab.txt -o gg/otutab.biom --table-type="OTU table" --to-json
biom summarize-table -i gg/otutab.biom
# 校正拷贝数
normalize_by_copy_number.py -i gg/otutab.biom -o temp/otutab_norm.biom -c /db/picrust/16S_13_5_precalculated.tab.gz
# 预测宏基因组KO表，biom方便下游归类，txt方便查看分析
predict_metagenomes.py -i temp/otutab_norm.biom -o gg/ko.biom -c /db/picrust/ko_13_5_precalculated.tab.gz
predict_metagenomes.py -f -i temp/otutab_norm.biom -o gg/ko.txt  -c /db/picrust/ko_13_5_precalculated.tab.gz

# KO转换为spf格式用于stamp分析
sed  -i '/# Constru/d;s/#OTU //' gg/ko.txt # 删除表头多余行及修正表头
# 将最后一列注释调整为第一列即为stamp使用的spf格式，大家也可用excel手动调整
num=`head -n1 gg/ko.txt|wc -w`
paste <(cut -f $num gg/ko.txt) <(cut -f 1-$[num-1] gg/ko.txt) > gg/ko.spf
# 现在可以下载ko.spf，配合metadata.txt，使用stamp统计绘图

# 按功能级别分类汇总, -c指输出类型KEGG_Pathways，可合并为1-3级
for i in 1 2 3;do
  categorize_by_function.py -f -i gg/ko.biom -c KEGG_Pathways -l ${i} -o gg/ko${i}.txt
  sed  -i '/# Constru/d;s/#OTU //' gg/ko${i}.txt
  paste <(cut -f $num gg/ko${i}.txt) <(cut -f 1-$[num-1] gg/ko${i}.txt) > gg/ko${i}.spf
done
wc -l gg/ko*.spf
# ko1只有8类, ko2有42个类，推荐ko3级别统计分析：即不多，又不少








# 二、扩增子无参分析流程 16S Amplicon Pipeline -- De novo (选学内容，推荐看扩增子视频教程)


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

# # 获得RDP物种注释
# usearch -sintax gg/otu.fa -db /db/usearch/rdp_16s_v16_sp.fa \
#   -strand both -tabbedout gg/sintax.txt -sintax_cutoff 0.6
# # 各分类级汇总表
# for i in p c o f g;do
#   usearch -sintax_summary gg/sintax.txt \
#   -otutabin gg/otutab.txt \
#   -rank ${i} \
#   -output gg/tax_sum_${i}.txt
# done


