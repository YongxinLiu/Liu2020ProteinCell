[TOC]

# 一、分析准备工作

    # 宏基因组分析流程：版本1.06，更新日期2019/11/1
    # 系统要求 System: Linux Ubuntu 18.04 / CentOS 7
    # 依赖软件 Sofware: Rstudio server 1.2、KneadData v0.6.1、
    #                   MetaPhlAn2 v2.7.5、HUMAnN2 v0.11.2 ......
    # Windows/Linux访问Rstudio服务器，推荐使用Chrome浏览器，可选Edge，IE有兼容问题
    
    # 运行前准备
    # 1. 第一次使用，请参见`db/script/soft_db.sh`脚本，
    #    安装分析所需软件和数据库(大约3-5天)
    
    # 2. 学员U盘和服务器家目录(~)有项目文件夹meta, 
    #    包含测序数据(seq/*.fq.gz)
    #    和实验设计(result/metadata.tsv);若无可手动从U盘上传。
    
    # 3. 使用谷歌浏览器访问服务器，具体IP地址课上通知，
    #    服务器可使用1个月
    
    # 4. Terminal中新建工作目录 mkdir -p meta ，右侧File中进入
    #    meta目录并上传流程pipeline.sh流程文件，完成后点击打开

## 1.0 每次重新进入必须运行的初始化步骤


    # **环境变量设置(每次开始分析前必须运行)**
    # 如果使用的变量未定义，则报错
    # set -u 
    # 设置数据库、软件和工作目录
    # 公共数据库database目录位置，如db公用可能为/db，而自己下载可能为~/db
    # ********
    db=/db
    # Conda软件software安装目录，如soft公用可能为/conda2，而自己下载可能为~/miniconda2
    soft=/conda2
    #export PATH=${soft}/bin:$PATH
    export PATH=${soft}/envs/metagenome_env/bin/:$PATH
    
    # wd为项目工作目录work directory，如meta
    wd=~/meta
    
    # 这里没有使用这个Rscript，是因为metagonome_env里面的R被降级了
    # 默认使用Rstudio的
    # Rscript="${soft}/bin/Rscript --vanilla"
    alias Rscript="/anaconda2/bin/Rscript --vanilla"
    # 创建并进入
    mkdir -p $wd 
    cd $wd
    # 加载宏基因组分析工作环境，新版conda命令为 conda activate
    source ${soft}/bin/activate metagenome_env

## 1.0.1 上传测序数据和实验设计文件

    # -p 表示如果文件夹已存在，则不作任何操作
    mkdir -p seq result tmp
    
    # # 用户filezilla上传测序文件至seq目录，本次从其它位置复制
    
    # # 6个样本数据，只取了10万条PE100数据(20MB)作为演示，通常>2千万条PE150数据(6GB)
    # /bin/cp -rf /db/meta/seq ./
    # -l 列出详细信息
    # -sh 显示文件大小 (s: size; h: human readable)
    ls -lsh seq
    # # 创建临时和结果目录
    # /bin/mkdir -p temp result
    # # 上传实验设计 metadata.tsv 于结果目录
    # /bin/cp -f /db/meta/result/design.txt result/metadata.tsv
    head result/metadata.tsv    


## 1.1 了解工作目录和文件

    # 显示文件结构， -L 2只显示两层
    tree -L 2
    # .
    # ├── pipeline.sh
    # ├── result
    # │   └── metadata.tsv
    # ├── seq
    # │   ├── p136C_1.fq.gz
    # │   ├── ....._2.fq.gz
    # │   └── p153N_2.fq.gz
    # ├── soft_db.sh
    # └── temp
    # pipeline.sh是分析流程代码；
    # seq目录中有6个样本双端测序，共12个序列文件；
    # temp是临时文件夹，存储分析中间文件，结束可全部删除节约空间
    # result是重要节点文件和整理化的分析结果图表，
    # 实验设计metadata.tsv也在此
    
    # 查看测序数据压缩文件的前8行
    
    # | 为管道符，上一个命令的输出，传递给下一个命令做输入
    # gzip: stdout: Broken pipe：管道断开。这里是人为断开，不是错误
    # 运行过程中需要仔细阅读屏幕输出的信息
    # ********
    zcat seq/p153C_2.fq.gz | head -n 8


## 1.2 FastQC质量评估(可选)

    # time统计运行时间，fastqc质量评估，*.gz为原始数据，
    # -t指定多线程，30S
    
    time fastqc seq/*.gz -t 1
    # 结果见seq目录，解读见
    
    # [数据的质量控制软件——fastQC](https://mp.weixin.qq.com/s/tDMih7ISLJcL4F4sWBq3Vw)
    
    # 生成多样品报告比较
    multiqc -f -d seq/ -o result/qc
    
    # 查看右侧Files - Home>meta>result/qc目录中multiqc_report.html，
    # 鼠标左键单击文件，弹出菜单中选择"View in web browser" 可查看交互式报告
    # 正常N和癌症C组GC都不同，组间差异太明显了


## 1.3 序列质控和去宿主 Qaulity control and remove host contamination

    # kneaddata是工作流程，依赖trimmomatic质控和去接头，
    # bowtie2比对宿主并筛选非宿主序列
    # 可只选一行中部分代码点击Run，如选中下行中#号后面命令查看程序帮助
    # kneaddata -h # 显示帮助
    
    # 以p144C单样品质控为例
    # 请务必根据自己软件和数据库安装位置，
    # 修改trimmomatic安装位置和数据库下载位置
    # 多行注释命令运行，可全选，按Ctrl+Shift+C进行注释的取消和添加
    
    # 输入文件：双端FASTQ测序数据，提供给参数-i
    # 输入文件：宿主基因组索引 -db ${db}/kneaddata/human_genome/Homo_sapiens
    # 输出文件：质控后的FASTQ测序数据，在目录temp/qc下面
    #           p153C_1_kneaddata_paired_1.fastq
    #           p153C_1_kneaddata_paired_1.fastq 用于后续分析
    
    time kneaddata -i seq/p144C_1.fq.gz -i seq/p144C_2.fq.gz \
      -o temp/qc -v -t 3 --remove-intermediate-output \
      --trimmomatic ${soft}/envs/metagenome_env/share/trimmomatic/ \
      --trimmomatic-options 'ILLUMINACLIP:${soft}/envs/metagenome_env/share/trimmomatic/adapters/TruSeq3-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:50' \
      --bowtie2-options '--very-sensitive --dovetail' -db ${db}/kneaddata/human_genome/Homo_sapiens
    # 10万条序列质控时间，10-30S
    
    # 查看质控后的结果文件
    ls -sh temp/qc/p144C*
    
    # 现实中是有一大堆样品，你可以逐个修改样品名运行，
    # 如果服务器性能允许可以并行加速分析
    # parallel --citation # 打will cite以后不再提醒
    
    # 运行下面这行，体会下parallel的工作原理
    # ::: 表示传递参数；第一个::: 后面为第一组参数，对应于{1};
    # 第一个::: 后面为第二组参数，对应于{2}，依次替换
    # parallel -j 3 --xapply "echo {1} {2}" ::: seq/*_1.fq.gz ::: seq/*_2.fq.gz
    
    mkdir -p temp/qc
    # 手动设置trimmomatic程序位置和接头位置 /conda2/envs/metagenome_env/share/trimmomatic/
    # 查看multiqc结果中接头污染最严重的样本，再到fastqc报告中查看接头序列，复制前20个碱基检索确定接头文件
    # grep 'GATCGGAAGAGCACACGTCT' /conda2/envs/metagenome_env/share/trimmomatic/adapters/*
    # 根据实际情况选择单端SE或双端PE，一般为 TruSeq3-PE-2.fa，更准确的是问测序公司要接头文件
    time parallel -j 3 --xapply \
      "kneaddata -i {1} -i {2} \
      -o temp/qc -v -t 1 --remove-intermediate-output \
      --trimmomatic ${soft}/envs/metagenome_env/share/trimmomatic/ \
      --trimmomatic-options 'ILLUMINACLIP:${soft}/envs/metagenome_env/share/trimmomatic/adapters/TruSeq3-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:50' \
      --bowtie2-options '--very-sensitive --dovetail' \
      -db ${db}/kneaddata/human_genome/Homo_sapiens" \
      ::: seq/*_1.fq.gz ::: seq/*_2.fq.gz
    # real 1m54s, user 7m14s

    # 质控结果统计
    kneaddata_read_count_table --input temp/qc --output result/01kneaddata_sum.txt
    cat result/01kneaddata_sum.txt


## 1.4 质控后质量再评估 (可选)

    fastqc temp/qc/*_1_kneaddata_paired_* -t 2
    multiqc -f -d temp/qc/ -o result/qc/
    # 整理bowtie2, trimmomatic, fastqc报告，接头和PCR污染率一般小于1%



# 二、宏基因组有参分析流程 Reference-based Metagenomic Pipeline (HUMAnN2)

    # 中文教程：https://mp.weixin.qq.com/s/XkfT5MAo96KgyyVaN_Fl7g
    # 英文教程：https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-Tutorial-(Humann2)


## 2.1 合并质控文件为HUMAnN2输入

    # 有参宏基因组不考虑双端，HUMAnN2要求将双端序列合并的单个文件作为输入文件
    mkdir -p temp/concat
    
    # for循环根据实验设计样本名批量双端序列合并
    # 注意星号和问号，当然大家更不能溜号~~~
    
    # tail -n+2 result/metadata.tsv
    # tail -n+2 result/metadata.tsv | cut -f 1
    # for i in `tail -n+2 result/metadata.tsv | cut -f 1`;do echo $i; done
    # ` 反引号为键盘右上角esc键下面的按键，一般在数字1的左边
    mkdir -p temp/concat
    for i in `tail -n+2 result/metadata.tsv | cut -f 1`; do 
      cat temp/qc/${i}*_1_kneaddata_paired_?.fastq > temp/concat/${i}.fq; 
    done
    
    # 查看样品数量和大小
    ls -sh temp/concat/*.fq 


## 2.2 HUMAnN2计算物种和功能组成


    # HUMAnN is a pipeline for efficiently and accurately profiling the presence/absence 
    # and abundance of microbial pathways in a community from metagenomic or metatranscriptomic 
    # sequencing data (typically millions of short DNA/RNA reads). This process, referred to as 
    # functional profiling, aims to describe the metabolic potential of a microbial community 
    # and its members. More generally, functional profiling answers the question 
    # "What are the microbes in my community-of-interest doing (or capable of doing)?"
    # 物种组成调用MetaPhlAn2, bowtie2比对至核酸序列；
    # 功能组成为humann2调用diamond比对至蛋白库11Gb
    
    # 单个文件运行示例，大约59 min
    # time humann2 --input temp/concat/p136C.fq  \
    #   --output temp/ \
    #   --threads 1
    
    # 并行HUMAnN2处理所有样品
    # parallel --citation # 输入提示文字
    # 运行需要75 min
    
    # 输入文件：temp/concat/*.fq 每个样品质控后双端合并后的fastq序列
    # 输出文件：temp/humann2/ 目录下
    #           p136C_pathabundance.tsv
    #           p136C_pathcoverage.tsv
    #           p136C_genefamilies.tsv
    # 整合后的输出：
    #           result/metaphlan2/taxonomy.tsv 物种丰度表
    #           result/metaphlan2/taxonomy.spf 物种丰度表（用于stamp分析）
    #           result/humann2/pathabundance_relab_stratified.tsv 通路丰度表
    #           result/humann2/pathabundance_relab_unstratified.tsv 通路丰度表
    #           stratified(每个菌的功能组成)和unstratified(功能组成)
    mkdir -p temp/humann2
    time parallel -j 2 \
      'humann2 --input {}  \
      --output temp/humann2/ --threads 8 ' \
      ::: temp/concat/*.fq > temp/log
    cat temp/log
    # 3 jobs X 8p; 时间1-2小时
    # 核心步骤，测序数据3X8=24线程，用时120min，真实数据可能要几小时至几天
    ls -shl temp/humann2/p136C*

## 2.3 物种组成表

	  # MetaPhlAn2整理了超过17000个参考基因组，包括13500个细菌和古菌，3500个病毒和110种真核生物，汇编整理了100万+类群特异的标记基因，可以实现：
	  # * 精确的分类群分配
	  # * 准确估计物种的相对丰度
	  # * 种水平精度
	  # * 株鉴定与追踪
    # * 超快的分析速度

    # MetaPhlAn is a computational tool for profiling the composition of microbial communities 
    # from metagenomic shotgun sequencing data. MetaPhlAn relies on unique clade-specific marker 
    # genes identified from ~17,000 reference genomes (~13,500 bacterial and archaeal, 
    # ~3,500 viral, and ~110 eukaryotic), allowing:

    # * up to 25,000 reads-per-second (on one CPU) analysis speed 
    #    (orders of magnitude faster compared to existing methods);
    # * unambiguous taxonomic assignments as the MetaPhlAn markers are clade-specific;
    # * accurate estimation of organismal relative abundance 
    #   (in terms of number of cells rather than fraction of reads);
    # * species-level resolution for bacteria, archaea, eukaryotes and viruses;
    # * extensive validation of the profiling accuracy on several synthetic datasets 
    #   and on thousands of real metagenomes.


    mkdir -p result/metaphlan2

### 2.3.1 样品结果合并

    merge_metaphlan_tables.py temp/humann2/*_humann2_temp/*_metaphlan_bugs_list.tsv | \
      sed 's/_metaphlan_bugs_list//g' > result/metaphlan2/taxonomy.tsv
    # 查看结果
    head -n6 result/metaphlan2/taxonomy.tsv

### 2.3.2 转换为spf格式方便stamp分析

    metaphlan_to_stamp.pl result/metaphlan2/taxonomy.tsv > result/metaphlan2/taxonomy.spf
    head -n6 result/metaphlan2/taxonomy.spf
    # 下载metadata.tsv和taxonomy.spf使用stamp分析

## 2.3.3 绘制热图

    # 无法运行时，推荐在Excel中筛选数据并使用ImageGP绘图
    metaphlan_hclust_heatmap.py --in result/metaphlan2/taxonomy.tsv \
        --out result/metaphlan2/heatmap.pdf \
        -c jet --top 30 --minv 0.1 -s log
    # c设置颜色方案，top设置物种数量，minv最小相对丰度，
    # s标准化方法，log为取10为底对数，文件名结尾可先pdf/png/svg三种图片格式。
    # 更多说明详见 metaphlan_hclust_heatmap.py -h

## 2.4 功能组成分析

    mkdir -p result/humann2
    # 合并所有样品通路丰度(pathabundance)，包括各功能和具体的物种组成，
    # 还有基因家族(genefamilies 太多)，通路覆盖度(pathcoverage)层面可以分析
    
    # gene family abundance is reported in RPK (reads per kilobase). 
    # This is computed as the sum of the scores for all alignments for a gene family. 
    # An alignment score is based on the number of matches to the reference gene 
    # for a specific sequence. It is divided by the length of the reference gene in kilobases 
    # to normalize for gene length. Each alignment score is also normalized 
    # to account for alignments for a single sequence to multiple reference genes. 
    # Alignments are not considered if they do not pass the e-value, identity, and coverage thresholds.
    
    # Gene families are groups of evolutionarily-related protein-coding sequences 
    # that often perform similar functions.
    # Gene family abundance at the community level is stratified to show 
    # the contributions from known and unknown species. Individual species' abundance contributions 
    # sum to the community total abundance.
    # HUMAnN2 uses the MetaPhlAn2 software along with the ChocoPhlAn database and 
    # translated search database for this computation.
    # Gene family abundance is reported in RPK (reads per kilobase) units to normalize for gene length; 
    # RPK units reflect relative gene (or transcript) copy number in the community. 
    # RPK values can be further sum-normalized to adjust for differences in sequencing depth across samples.

    # The "UNMAPPED" value is the total number of reads which remain unmapped after 
    # both alignment steps (nucleotide and translated search). 
    # Since other gene features in the table are quantified in RPK units, "UNMAPPED" can be interpreted as
    # a single unknown gene of length 1 kilobase recruiting all reads that failed to map to known sequences.
    
    # The UniRef50_unknown values represent the total abundance of reads which map to 
    # ChocoPhlAn nucleotide sequences which do not have a UniRef50 annotation.
    
    # result/humann2/pathabundance.tsv
    # This file details the abundance of each pathway in the community as a function of 
    # the abundances of the pathway's component reactions, with each reaction's abundance 
    # computed as the sum over abundances of genes catalyzing the reaction.

    # Pathway abundance is computed once at the community level and again for each species 
    # (plus the "unclassified" stratum) using community- and species-level gene abundances 
    # along with the structure of the pathway.

    # The pathways are ordered by decreasing abundance with pathways for each species also sorted by
    # decreasing abundance. Pathways with zero abundance are not included in the file.

    # Pathway abundance is proportional to the number of complete "copies" of the pathway in the community. 
    # Thus, for a simple linear pathway RXN1â†’RXN2â†’RXN3â†’RXN4, 
    # if RXN1 is 10 times as abundant as RXNs 2-4, the pathway abundance will be 
    # driven by the abundances of RXNs 2-4.

    # Unlike gene abundance, a pathway's community-level abundance is not necessarily 
    # the sum of its stratified abundance values. For example, continuing with the simple linear pathway example
    # introduced above, if the abundances of RXNs 1-4 are [5, 5, 10, 10] in Species_A and 
    # [10, 10, 5, 5] in Species_B, HUMAnN2 would report that Species_A and Species_B each 
    # contribute 5 complete copies of the pathway. 
    # However, at the community level, the reaction totals are [15, 15, 15, 15], and 
    # thus HUMAnN2 would report 15 complete copies.

    # In greater detail, the abundance for each pathway is a recursive computation of 
    # abundances of sub-pathways with paths resolved to abundances based on the relationships
    # and abundances of the reactions contained in each. Each path, the smallest portion of 
    # a pathway or sub-pathway which can't be broken down into sub-pathways, 
    # has an abundance that is the max or harmonic mean of the reaction abundances 
    # depending on the relationships of these reactions. Optional reactions are only added to 
    # the overall abundance if their abundance is greater than the harmonic mean of the required reactions.

    # Gap filling allows for a single required reaction to have a zero abundance. 
    # For all pathways, the required reaction with the lowest abundance is replaced with 
    # the abundance of the required reaction with the second lowest abundance.
    # By default, HUMAnN2 uses MetaCyc pathway definitions and MinPath to identify 
    # a parsimonious set of pathways which explain observed reactions in the community.

    # The user has the option to provide a custom pathways database to HUMAnN2 and 
    # to use all pathways instead of the minimal pathways computed by MinPath.
    # To account for non-linearity in the conversion of gene copy number to pathway copy number, 
    # we define a "compression constant" (k) equal to the total pathway abundance divided by 
    # the total abundance of genes that contributed to pathways. 
    
    # The "UNMAPPED" value reported in the pathway abundance table is equal to 
    # the total number of unmapped reads scaled by k 
    # (making it more comparable with pathway abundance values). 
    
    # Similarly, we define an "UNINTEGRATED" abundance for 
    # 1) the community, 2) each identified species, and 3) the "unclassified" stratum 
    # equal to the total abundance of genes in that level that did not contribute 
    # to pathways (scaled by k).
    # "UNINTEGRATED" does not appear for stratifications with no detected pathways.
    
    humann2_join_tables --input temp/humann2 --file_name pathabundance \
      --output result/humann2/pathabundance.tsv
    # 预文件头尾格式
    head -n3 result/humann2/pathabundance.tsv
    tail -n3 result/humann2/pathabundance.tsv
    # 删除列名多余信息
    sed -i 's/_Abundance//g' result/humann2/pathabundance.tsv

    # 标准化为相对丰度relab或百万分数cpm
    humann2_renorm_table --input result/humann2/pathabundance.tsv \
      --units relab \
      --output result/humann2/pathabundance_relab.tsv
    head -n5 result/humann2/pathabundance_relab.tsv
    
    # 分层结果：功能对应物种表(stratified)和功能组成表(unstratified)
    humann2_split_stratified_table --output result/humann2/ \
      --input result/humann2/pathabundance_relab.tsv 
    # Split stratified tables created:
    # result/humann2/pathabundance_relab_stratified.tsv
    # result/humann2/pathabundance_relab_unstratified.tsv
    head -n3 result/humann2/pathabundance_relab_stratified.tsv
    head -n3 result/humann2/pathabundance_relab_unstratified.tsv
    # 结果stratified(每个菌的功能组成)和unstratified(功能组成)两个
    # 可以使用stamp对分层，或末分层的结果进行统计分析

    # 具体解释 https://bitbucket.org/biobakery/humann2/wiki/Home#markdown-header-standard-workflow

### 2.4.1 功能组成分析脚本——差异比较
    
    # 输入数据：通路丰度表格 result/humann2/pathabundance_relab_stratified.tsv
    # 输入数据：实验设计信息 result/metadata.tsv
    # 中间数据：PCL文件 如 result/humann2/pathabundance_relab_stratified.pcl 
    #           是包含分组信息的通路丰度表格文件
    # 生成下游分析的pcl文件 ********
    # 如果是自己的数据，则把下面8行的注释去掉
    # 从humann2_associate到humann2_barplot中文件名hmp_pathabund.pcl替换为pathabundance_relab_stratified.pcl
     # awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$1]=$2;}ARGIND==2\
     # {if(FNR==1) {print $0; group="STSite"; for(i=2;i<=NF;i++) group=group"\t"a[$i]; print group;} \
     # else {print $0}}' result/metadata.tsv result/humann2/pathabundance_relab_stratified.tsv \
     # >result/humann2/pathabundance_relab_stratified.pcl
     # 
     # humann2_associate --input result/humann2/pathabundance_relab_stratified.pcl  \
     #    --focal-metadatum STSite --focal-type categorical --last-metadatum STSite --fdr 0.2 \
     #    --output result/humann2/pathabundance_relab_stratified.stats.txt
    
    # https://bitbucket.org/biobakery/biobakery/wiki/humann2#rst-header-downstream-analyses
    # 使用示例数据信息更完整，样本名下面有分组信息，需要手动制作
    # 差异比较：输入input、分组focal、分组类型type、分组行结果、FDR和输出output
    
    # When associating metadata with HUMAnN2 features, it is often beneficial to associate with 
    # community totals and avoid testing each individual feature stratification 
    # (to improve statistical power). 
    # This is the approach used by the HUMAnN2 utility script humann2_associate, 
    # which can compare feature totals across samples with 
    # 1) a single continuous metadatum (via the Spearman Correlation) or 
    # 2) a single categorical metadatum (via the Kruskal-Wallis H-test). 
    # Notably, this is a naive approach to association, but it is useful for tutorial purposes.
    # https://bitbucket.org/biobakery/biobakery/wiki/humann2
    
    humann2_associate --input result/humann2/hmp_pathabund.pcl  \
        --focal-metadatum STSite --focal-type categorical --last-metadatum STSite --fdr 0.2 \
        --output result/humann2/stats.txt
    
    # first pathway in the table (METSYN-PWY: L-homoserine and L-methionine biosynthesis) 
    # tends to be more abundant on average at the three oral sites (buccal mucosa, tongue, and plaque) 
    # than at the three non-oral sites (stool, posterior fornix, and nares).
    # 结果为4列的文件：通路名字，通路在各个分组的丰度，差异P-value，校正后的Q-value
    head -n 3 result/humann2/stats.txt
    
### 2.4.2 功能组成分析脚本——绘制通路物种组成

    #输入文件：通路丰度信息表 result/humann2/hmp_pathabund.pcl
    #输出文件： result/humann2/plot*.pdf 图片
    
    # barplot展示各样品高丝氨酸合成通路的物种组成
    # focal-feature指定会绘制的通路，需要根据需求修改 ********
    humann2_barplot --input result/humann2/hmp_pathabund.pcl --focal-feature METSYN-PWY \
        --focal-metadatum STSite --last-metadatum STSite --output result/humann2/plot1.pdf
    # --sort sum 按丰度排序
    humann2_barplot --sort sum --input result/humann2/hmp_pathabund.pcl --focal-feature METSYN-PWY \
        --focal-metadatum STSite --last-metadatum STSite --output result/humann2/plot_order_by_abundance.pdf
    # --sort sum metadata 按分组 丰度排序
    humann2_barplot --sort sum metadata --input result/humann2/hmp_pathabund.pcl --focal-feature METSYN-PWY \
        --focal-metadatum STSite --last-metadatum STSite \
        --output result/humann2/plot_order_by_group_then_abundance.pdf


## 2.5 GraPhlAn图

    # 输入：result/metaphlan2/taxonomy.tsv
    # 输出：result/metaphlan2/graphlan.pdf
    # metaphlan2 to graphlan
    export2graphlan.py --skip_rows 1,2 -i result/metaphlan2/taxonomy.tsv \
      --tree temp/merged_abundance.tree.txt \
      --annotation temp/merged_abundance.annot.txt \
      --most_abundant 1000 --abundance_threshold 20 --least_biomarkers 10 \
      --annotations 3,4 --external_annotations 7
    # 参数说明见PPT，或运行 export2graphlan.py --help
    # graphlan annotation
    graphlan_annotate.py --annot temp/merged_abundance.annot.txt \
      temp/merged_abundance.tree.txt  temp/merged_abundance.xml
    # output PDF figure, annoat and legend
    graphlan.py temp/merged_abundance.xml result/metaphlan2/graphlan.pdf \
      --external_legends 


## 2.6 LEfSe差异分析和Cladogram

    # 输入文件：物种丰度表result/metaphlan2/taxonomy.tsv
    # 输入文件：样品分组信息 result/metadata.tsv
    # 整合后用于LefSe分析的文件：result/metaphlan2/lefse.txt 
    #            (这个文件可以提供给www.ehbio.com/ImageGP 用于在线LefSE分析)
    # LefSe结果输出：result/metaphlan2/目录下lefse开头和feature开头的文件
    
    # 修改样本品为组名，**非通用代码**，可手动修改
    # 只保留组名C, N
    head -n3 result/metaphlan2/taxonomy.tsv
    
    # 注意$2，代表metadata.tsv的第二列，也就是分组信息
    # 如果你的分组信息不在第二列则相应修改 ********
    awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$1]=$2;}ARGIND==2{if(FNR==1) \
      {for(i=2;i<=NF;i++) $i=a[$i];} print $0;}' result/metadata.tsv result/metaphlan2/taxonomy.tsv \
      | grep -v '#' > result/metaphlan2/lefse.txt
    
    # sed '1 s/p[0-9]*//g' result/metaphlan2/taxonomy.tsv \
    #   | grep -v '#' > result/metaphlan2/lefse.txt
    head -n3 result/metaphlan2/lefse.txt
    
    # LEfSe本地代码供参考，如Rstudio运行报错，可选Xshell下运行或ImageGP在线分析
    # 格式转换为lefse内部格式
    lefse-format_input.py result/metaphlan2/lefse.txt temp/input.in -c 1 -o 1000000
    # 运行lefse
    run_lefse.py temp/input.in temp/input.res
    # 绘制物种树注释差异
    lefse-plot_cladogram.py temp/input.res result/metaphlan2/lefse_cladogram.pdf --format pdf --dpi 600
    # 绘制所有差异features柱状图
    lefse-plot_res.py temp/input.res result/metaphlan2/lefse_res.pdf --format pdf --dpi 600
    # 绘制单个features柱状图(同STAMP中barplot)
    grep -v '-' temp/input.res | sort -k3,3n  # 查看显著差异features，按丰度排序
    lefse-plot_features.py -f one --feature_name "k__Bacteria.p__Firmicutes" \
      --format pdf \
      temp/input.in temp/input.res result/metaphlan2/Firmicutes.pdf
    # 批量绘制所有差异features柱状图
    lefse-plot_features.py -f diff --archive none --format pdf \
      temp/input.in temp/input.res result/metaphlan2/features


## 2.7 kraken2基于NCBI数据库注释reads层面

    # 还可以进行contig、gene、bin层面的序列物种注释

### 2.7.1 物种注释

    cd ${wd}
    mkdir -p temp/kraken2
    
    # 单样本注释，大约1-5m
    # time kraken2 --db ${db}/kraken2 --paired temp/qc/p136C_1_kneaddata_paired*.fastq \
    #   --threads 3 --use-names --use-mpa-style --report-zero-counts \
    #   --report temp/kraken2/p136C_report \
    #   --output temp/kraken2/p136C_output
    
    # 多样本并行
    # 固定输入：-db ${db}/kraken2
    # 样品输入：temp/qc/{1}_1_kneaddata_paired*.fastq 质控后的FASTQ数据
    # 输出结果：每个样本单独输出，temp/kraken2/{1}_report和temp/kraken2/{1}_output
    #           {1}代表样本名字
    # 整合后的输出结果： result/kraken2/taxonomy_count.txt 物种丰度表
    time parallel -j 2 \
      "kraken2 --db ${db}/kraken2 --paired temp/qc/{1}_1_kneaddata_paired*.fastq \
      --threads 3 --use-names --use-mpa-style --report-zero-counts \
      --report temp/kraken2/{1}_report \
      --output temp/kraken2/{1}_output" \
      ::: `tail -n+2 result/metadata.tsv | cut -f 1`
    # 屏幕会输出各样品注释比例，和运行时间 10m

### 2.7.2 汇总样品物种组成表

    mkdir -p result/kraken2
    # 输出结果行数相同，但不一定顺序一致，要重新排序
    parallel -j 1 \
      'sort temp/kraken2/{1}_report | cut -f 2 | sed "1 s/^/{1}\n/" > temp/kraken2/{1}_count ' \
      ::: `tail -n+2 result/metadata.tsv | cut -f 1`
    # 提取第一样本品行名为表行名
    header=`tail -n 1 result/metadata.tsv | cut -f 1`
    echo $header
    sort temp/kraken2/${header}_report | cut -f 1 | \
      sed "1 s/^/Taxonomy\n/" > temp/kraken2/0header_count
    head -n3 temp/kraken2/0header_count
    # paste合并样本为表格
    ls temp/kraken2/*count
    paste temp/kraken2/*count > result/kraken2/taxonomy_count.txt
    
    # 查看物种丰度表
    head result/kraken2/taxonomy_count.txt

### 2.7.3 物种多样性分析

    # 输入数据：kraken2鉴定出的物种丰度表 result/kraken2/taxonomy_count.txt
    # 输出结果：result/kraken2目录下的txt文件和pdf文件，每个代表一种多样性指数
    
    # 提取种级别注释并抽平至最小测序量，计算6种alpha多样性指数
    # Rscript ${db}/script/kraken2alpha.R -h # 查看帮助和参数详细
    # 输入count文件，按最小值抽平为norm文件，计算alpha多样性
    Rscript ${db}/script/kraken2alpha.R --input result/kraken2/taxonomy_count.txt \
      --depth 0 --normalize result/kraken2/taxonomy_count.norm \
      --output result/kraken2/taxonomy_count.alpha
    
    # 绘制Alpha多样性指数，结果为输入文件+类型richness/chao1/ACE/shannon/simpson
    Rscript ${db}/script/alpha_boxplot.R -i result/kraken2/taxonomy_count.alpha.txt \
      -t shannon -d result/metadata.tsv -n group -w 4 -e 2.5
    
    # 批量计算6种指数的箱线图+统计
    for i in richness chao1 ACE shannon simpson invsimpson;do
    Rscript ${db}/script/alpha_boxplot.R -i result/kraken2/taxonomy_count.alpha.txt \
      -t ${i} -d result/metadata.tsv -n group -w 4 -e 2.5
    done
    # 解读：richness和shannon无显著差异，chao1有显著差异，
    # 推测是测序深度不足(测试数据仅有0.25%数据)



# 三、宏基因组无参分析流程 De novo Metagenomic Pipeline (MEGAHIT)


## 3.1 拼接 Assembly

### 3.1.1 Megahit拼接

    # 输入数据：质控后的fastq文件
    # 输出数据：拼装好的序列文件 temp/megahit/final.contigs.fa

    # 删除旧文件夹，否则不能重新计算
    /bin/rm -rf temp/megahit
    # 组装核心程序，计算和内存需求密集，几百G或上T数据需几天至几周
    
    # 方法1. 混合组装。2线程2分钟
    # 优点是简单快速获得得一套参考序列，基因基本不用冗余，增加低丰度菌测序深度提高拼接长度；
    # 缺点需要更大内存，混样提高错误拼接风险
    time megahit -t 2 \
      -1 `ls temp/qc/*_1_kneaddata_paired_1.fastq|tr '\n' ','|sed 's/,$//'` \
      -2 `ls temp/qc/*_1_kneaddata_paired_2.fastq|tr '\n' ','|sed 's/,$//'` \
      -o temp/megahit 
    
    # 高级用户可指定kmer优化时间，最小值越大、步长越大可缩短计算时间 # --k-min 27 --k-max 141 --k-step 12
    # 默认k list: 21,29,39,59,79,99,119
    # `ls...`用于自动获得文件列表，无需手动添写，但是要记住模式
    # 默认使用所有线程，如本机96线程，时间2m, 系统机时192m
    # 设定9线程，时间3m, 系统机时19m；才快了1.5倍，但浪费了10倍资源或机时费用
    
    # 方法2. 单样品并行组装。2线程，4分钟
    # 优点是资源消耗少，防旧混样的相似菌错误组装；
    # 缺点是低丰度菌难组装较完整，样品间基因大量冗余，需要更多的计算时间；
    # time parallel -j 1 \
    # 'time megahit -t 2 \
    #  -1 temp/qc/{1}_1_kneaddata_paired_1.fastq \
    #  -2 temp/qc/{1}_1_kneaddata_paired_2.fastq \
    #  -o temp/megahit/{1}' \
    #  ::: `tail -n+2 result/metadata.tsv | cut -f 1`
    # 有条件可以两种都做，在基因层面去冗余获得最完整的基因集

    # 以下按混合组装方式演示查看拼接结果
    ls -sh temp/megahit/final.contigs.fa
    # 最前2行，前60列字符预览
    head -n2 temp/megahit/final.contigs.fa | cut -c1-60

### 3.1.2 metaSPAdes拼接(备选1)

    # 方法1. 混合组装。9线程7分钟
    time metaspades.py -t 9 -m 500 \
      `ls temp/qc/*_1_kneaddata_paired_1.fastq|sed 's/^/-1 /'| tr '\n' ' '` \
      `ls temp/qc/*_1_kneaddata_paired_2.fastq|sed 's/^/-2 /'| tr '\n' ' '` \
      -o temp/metaspades
    ll metaspades/contigs.fasta
    
    # 方法2. 单样品并行组装。9线程，4分钟
    time parallel -j 1 \
    'metaspades.py -t 9 -m 500 \
      -1 temp/qc/{1}_1_kneaddata_paired_1.fastq \
      -2 temp/qc/{1}_1_kneaddata_paired_2.fastq \
      -o temp/metaspades/{1}' \
      ::: `tail -n+2 result/metadata.tsv | cut -f 1`

### 3.1.3 quast评估

    # 输入文件：拼装好的序列文件 result/megahit/final.contigs.fa
    # 输出文件：result/megahit/quast/report.html 
    
    mkdir -p result/megahit/
    ln -f temp/megahit/final.contigs.fa result/megahit/
    quast.py result/megahit/final.contigs.fa -o result/megahit/quast -t 2
    # 生成report文本tsv/txt、网页html、PDF等格式报告
    
    # 可选metaquast评估，更全面，但需下载相关数据库，受网速影响可能时间很长
    # metaquast based on silva, and top 50 species genome to access
    time metaquast.py result/megahit/final.contigs.fa -o result/megahit/metaquast


## 3.2 基因注释和定量 Gene annotation & quantitfy
    
    # 输入文件：拼装好的序列文件 result/megahit/final.contigs.fa
    # 总结输出文件：Proka预测的基因序列 temp/prokka/mg.ffn
    # 去冗余后的基因和蛋白序列：result/NR/nucleotide.fa
    #                           result/NR/protein.fa
    # Salmon定量后的结果：result/salmon/gene.count
    #                     result/salmon/gene.TPM
    
### 3.2.1 prokka基因注释

    # 查看文件大小，预估时间
    ls -sh temp/megahit/final.contigs.fa # 3.0 Mb
    # 启动虚拟环境，新版本 source 替换为 conda
    # prokka依赖关系复杂，我们使用metawrap中的prokka
    # 直接使用时，当前环境Java版本不符，只能使用另一个虚拟环境
    source ${soft}/bin/activate metawrap
    # 找不到Can't locate XML/Simple.pm，手动指定Perl library
    # locate XML/Simple.pm 找模块位置
    # find . -name Module name
    export PERL5LIB=$PERL5LIB:${soft}/envs/metawrap/lib/perl5/site_perl/5.22.0
    time prokka temp/megahit/final.contigs.fa --outdir temp/prokka \
      --prefix mg --metagenome --kingdom Archaea,Bacteria,Mitochondria,Viruses \
      --force --cpus 3
    # 退出虚拟环境; 退回metagenome_env
    conda deactivate
    # 3线程，耗时2m

### 3.2.2 cd-hit构建非冗余基因集

    mkdir -p temp/NR
    # 统计注释基因数量1948或1950
    grep -c '>' temp/prokka/mg.ffn 
    # 输入文件可由多个样本、组、批次序列合并文件，方便整合分析
    # aS覆盖度，c相似度，G局部比对，M内存0不限制，T多线程，g最优解
    time cd-hit-est -i temp/prokka/mg.ffn -o temp/NR/mg.ffn.nr \
        -aS 0.9 -c 0.95 -G 0 -M 0 -T 9 -g 1
    
    # 统计非冗余基因数量1938
    grep -c '>' temp/NR/mg.ffn.nr
    # 数量下降不大，因为拼接内部去除冗余，多批拼接效果明显
    # 挑选重要结果保存到结果目录
    mkdir -p result/NR
    ln temp/NR/mg.ffn.nr result/NR/nucleotide.fa
    # 翻译核酸为对应蛋白序列， emboss
    transeq -sequence result/NR/nucleotide.fa -outseq result/NR/protein.fa -trim Y 
    # 序列名自动添加了_1，为与核酸对应要去除
    sed -i 's/_1 / /' result/NR/protein.fa

### 3.2.3 基因定量salmon

    mkdir -p temp/salmon
    # 切换到虚拟环境中使用salmon
    # 建索引,-t转录本，--type类型fmd/quasi，-k kmer长度默认31, -i 索引
    # 直接运行找不到库，可使用程序完整路径解决问题(error while loading shared libraries: liblzma.so.0)
    ${soft}/envs/metagenome_env/share/salmon/bin/salmon index -t result/NR/nucleotide.fa -p 9 -k 31 \
        -i temp/salmon/index # --type quasi 
    
    # 定量，l文库类型自动选择，p线程 ，--meta宏基因组模式, 3个任务并行，共20秒
    # 注意parallel中待并行的命令必须是双引号 ********
    time parallel -j 3 \
        "${soft}/envs/metagenome_env/share/salmon/bin/salmon quant -i temp/salmon/index -l A -p 3 --meta \
        -1 temp/qc/{1}_1_kneaddata_paired_1.fastq \
        -2 temp/qc/{1}_1_kneaddata_paired_2.fastq \
        -o temp/salmon/{1}.quant" \
        ::: `tail -n+2 result/metadata.tsv | cut -f 1`
    
    # 合并
    mkdir -p result/salmon
    ${soft}/envs/metagenome_env/share/salmon/bin/salmon quantmerge \
        --quants temp/salmon/*.quant \
        -o result/salmon/gene.TPM
    ${soft}/envs/metagenome_env/share/salmon/bin/salmon quantmerge \
        --quants temp/salmon/*.quant \
        --column NumReads -o result/salmon/gene.count
    sed -i '1 s/.quant//g' result/salmon/gene.*
    # 查看结果表格
    head -n3 result/salmon/gene.*


## 3.3 功能基因注释

    # 输入数据：上一步预测的蛋白序列 result/NR/protein.fa
    # 中间结果：temp/eggnog/protein.emapper.seed_orthologs
    #           temp/eggnog/output.emapper.annotations
    #           temp/eggnog/output
    
    # COG定量表：result/eggnog/cogtab.count
    #            result/eggnog/cogtab.count.spf (用于STAMP)
    
    # KO定量表：result/eggnog/kotab.count
    #           result/eggnog/kotab.count.spf  (用于STAMP)
    
    # CAZy碳水化合物注释和定量：result/dbcan2/cazytab.count
    #                           result/dbcan2/cazytab.count.spf (用于STAMP)
    
    # 抗生素抗性：result/resfam/resfam.count
    #             result/resfam/resfam.count.spf (用于STAMP)
    
    # 这部分可以拓展到其它数据库

### 3.3.1 基因注释eggNOG/COG/KEGG/GO

    # diamond比对基因至eggNOG数据库, 用时97 m (minutes)
    mkdir -p temp/eggnog
    time emapper.py -m diamond --no_annot --no_file_comments \
      --data_dir ${db}/eggnog --cpu 2 -i result/NR/protein.fa \
      -o temp/eggnog/protein --override
    
    # 比对结果功能注释 6 minutes
    time emapper.py --annotate_hits_table \
      temp/eggnog/protein.emapper.seed_orthologs --no_file_comments \
    	-o temp/eggnog/output --cpu 2 --data_dir ${db}/eggnog --override
    
    # 结果注释表头, 重点1序列名，7KO，12COG分类，13注释
    mkdir -p result/eggnog
    sed '1 i Name\teggNOG\tEvalue\tScore\tGeneName\tGO\tKO\tBiGG\tTax\tOG\tBestOG\tCOG\tAnnotation' \
      temp/eggnog/output.emapper.annotations > temp/eggnog/output

    # 1.整理COG表
    # 提取12列COG分类
    cut -f 1,12 temp/eggnog/output|cut -f 1 -d ','|grep -v -P '\t$' \
      >temp/eggnog/1cog.list
    # 基因丰度矩阵末尾添加对应cog编号，没注释的直接删除
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' \
      temp/eggnog/1cog.list result/salmon/gene.count | \
    	sed '/\t$/d' | sed '1 s/COG/KO/' > temp/eggnog/gene_cog.count
    # 按COG类型合并表格，输出count和RPM值，n设置标准化单位，默认1M，可选100/1
    ${Rscript} ${db}/script/mat_gene2ko.R -i temp/eggnog/gene_cog.count \
      -o result/eggnog/cogtab -n 1000000
    # STAMP的spf格式，结果metadata.tsv进行KO或Description差异比较
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2"\t"$3} NR>FNR{print a[$1],$0}' \
      ${db}/eggnog/COG.anno result/eggnog/cogtab.count > \
      result/eggnog/cogtab.count.spf
    
    # 2. 整理KO表
    # 提取基因KO表，基因1对多个KO时只提取第一个KO
    cut -f 1,7 temp/eggnog/output|cut -f 1 -d ','|grep -v -P '\t$' \
      > temp/eggnog/2ko.list
    # 基因丰度矩阵末尾添加对应KO编号，没注释的直接删除
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' \
      temp/eggnog/2ko.list \
      result/salmon/gene.count | sed '/\t$/d' > temp/eggnog/gene_ko.count
    # 合并基因表为KO表，输出count值和tpm值
    ${Rscript} ${db}/script/mat_gene2ko.R -i temp/eggnog/gene_ko.count \
      -o result/eggnog/kotab -n 1000000
    # STAMP的spf格式，结果metadata.tsv进行KO或Description差异比较
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
      ${db}/eggnog/KO.anno result/eggnog/kotab.count | \
      sed 's/^\t/Undescription\t/' > result/eggnog/kotab.count.spf
      
    head result/eggnog/kotab.count.spf


### 3.3.2 碳水化合物dbCAN2

    # 比对CAZy数据库, 用时1m; sensitive 31m
    mkdir -p temp/dbcan2
    # Error: Error opening file /db/dbCAN2/CAZyDB.07312018.dmnd，没有设置好数据库位置
    time diamond blastp --db ${db}/dbCAN2/CAZyDB.07312018 --query result/NR/protein.fa \
    	--outfmt 6 --threads 2 --max-target-seqs 1 --quiet -e 1e-5 --sensitive \
    	--out temp/dbcan2/gene_diamond.f6 
    # 整理比对数据为表格
    mkdir -p result/dbcan2
    # 提取基因对应基因家族，同一基因存在1对多，只取第一个
    cut -f 1,2 temp/dbcan2/gene_diamond.f6 | uniq | sed 's/|/\t/g' | cut -f 1,3 | \
    	cut -f 1,2 -d '_' |sed '1 i Name\tKO' > temp/dbcan2/gene_fam.list
    # 基因丰度矩阵末尾添加对应FAM编号，没注释的直接删除
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print $0,a[$1]}' temp/dbcan2/gene_fam.list \
      result/salmon/gene.count | sed '/\t$/d' > temp/dbcan2/gene_fam.count
    # 按基因家族合并
    # /anaconda2/bin/Rscript --vanilla ${db}/script/mat_gene2ko.R -i temp/dbcan2/gene_fam.count -o result/dbcan2/cazytab
    ${Rscript} ${db}/script/mat_gene2ko.R -i temp/dbcan2/gene_fam.count -o result/dbcan2/cazytab
    # 结果中添加FAM注释，spf格式用于stamp分析
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' ${db}/dbCAN2/fam_description.txt \
    	result/dbcan2/cazytab.count > result/dbcan2/cazytab.count.spf


### 3.3.3 抗生素抗性ResFam

    mkdir -p temp/resfam
    # 比对至抗生素数据库 1s, 8s
    time diamond blastp --db ${db}/resfam/Resfams-proteins --query result/NR/protein.fa \
    	--outfmt 6 --threads 2 --max-target-seqs 1 --quiet -e 1e-5 --sensitive \
    	--out temp/resfam/gene_diamond.f6
    
    mkdir -p result/resfam
    # 提取基因对应基因家族
    cut -f 1,2 temp/resfam/gene_diamond.f6 | uniq | \
      sed '1 i Name\tResGeneID' > temp/resfam/gene_fam.list
    # 基因丰度矩阵末尾添加对应FAM编号，没注释的直接删除
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' \
      temp/resfam/gene_fam.list result/salmon/gene.count | \
    	sed '/^\t/d' > result/resfam/resfam.count
    
    # 统计注释基因的比例
    wc -l result/salmon/gene.count
    wc -l result/resfam/resfam.count # 172/7734=2.2%
    
    # 结果中添加FAM注释，spf格式用于stamp分析
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$4"\t"$3"\t"$2} NR>FNR{print a[$1],$0}' \
      ${db}/resfam/Resfams-proteins_class.tsv  result/resfam/resfam.count \
      > result/resfam/resfam.count.spf


### 3.3.4 抗生素抗性CARD

    # (推荐在线，可选)




# 四、挖掘单菌基因组/分箱 Binning (MetaWRAP)——官方测试数据

    # 主页：https://github.com/bxlab/metaWRAP
    # 挖掘单菌基因组，需要研究对象复杂度越低、测序深度越大，结果质量越好。
    # 推荐单样本数据量30GB+，至少3个样本起。
    # 上面的测序数据6个样一共有70MB，适合演示分析流程但无法挖掘到单菌基因组，
    # 这里使用官方测序数据演示讲解
    # 软件和数据库布置需2-3天，演示数据分析过程超10h，标准30G样也需3-30天，由服务器性能决定。
    # 流程: https://github.com/bxlab/metaWRAP/blob/master/Usage_tutorial.md

    # 输入数据：质控后的FASTQ序列，文件名格式必须为*_1.fastq和*_2.fastq
    #           p153C_1_kneaddata_paired_1.fastq  -> p153C_1_1.fq
    #           p153C_1_kneaddata_paired_2.fastq  -> p153C_1_2.fq
    #           放置到 binning/temp/qc 目录下
    
    # 拼装获得的contig文件：result/megahit/final.contigs.fa
    #           放置到 binning/temp/megahit 目录下
    #           

## 4.0 准备数据和环境变量



    # 准备原始数据从头分析，详见公众号或官网教程
    # 这里我们从质控后数据和拼接结果开始
    cd ${wd}
    mkdir -p binning && cd binning
    mkdir -p temp && cd temp
    # 这里基于质控clean数据和拼接好的contigs，自己链接自上游分析
    # 7G质控数据，输入数据文件名格式必须为*_1.fastq和*_2.fastq
    mkdir -p qc
    cd qc
    # 方法1. 下载测序数据
    # for i in `seq 7 9`;do
    #     wget -c http://210.75.224.110/share/meta/metawrap/ERR01134${i}_1.fastq.gz
    #     wget -c http://210.75.224.110/share/meta/metawrap/ERR01134${i}_2.fastq.gz
    # done
    # gunzip *.gz # 解压文件
    # rename .fq .fastq *.fq # 批量修改扩展名
    # 方法2. 复制准备好的数据
    ln -s ${db}/metawrap/*.fastq ./
    cd ..
    # megahit拼接结果
    mkdir -p megahit
    cd megahit
    # wget -c http://210.75.224.110/share/meta/metawrap/final.contigs.fa.gz
    # gunzip *.gz
    ln -s ${db}/metawrap/*.fa ./
    
    # 加载运行环境
    cd ${wd}/binning
    conda activate metawrap


## 4.1 运行三种bin软件

    # 输入文件为contig和clean reads
    # 调用三大主流binning程序cococt, maxbin2, metabat2
    # 8线程耗时 0.5 - 2 小时
    # nohup 和 & 保证任务在后台不被中断，且记录输出内容到 nohup.out(可选)
    nohup metawrap binning -o temp/binning -t 3 -a temp/megahit/final.contigs.fa \
      --metabat2 --maxbin2 --concoct temp/qc/ERR*.fastq &
    # 用自己的文件，替换输出文件名为 *1_kneaddata_paired*.fastq 
    # 输出文件夹 temp/binning 包括3种软件结果和中间文件
    # OpenBLAS Warning : Detect OpenMP Loop and this application may hang. Please rebuild the library with USE_OPENMP=1 option.


## 4.2 Bin提纯

    # 一般要求c完整度70，x污染率5，这里数据少降低阈值保证有结果作演示
    # 8线程耗时 1 - 2 小时
    cd ${wd}/binning
    # /bin/rm -rf temp/bin_refinement
    nohup metawrap bin_refinement -o temp/bin_refinement -t 3 \
      -A temp/binning/metabat2_bins/ \
      -B temp/binning/maxbin2_bins/ \
      -C temp/binning/concoct_bins/ \
      -c 50 -x 10 &
    # 查看高质量Bin的数量
    # cat temp/bin_refinement/metawrap_bins.stats | awk '$2>50 && $3<10' | wc -l
    wc -l temp/bin_refinement/metawrap_50_10_bins.stats
    
    head temp/bin_refinement/metawrap_50_10_bins.stats
    # 结果改进程度见temp/bin_refinement/figures/目录


## 4.3 Bin定量

    # 使用salmon计算每个bin在样本中相对丰度
    # 耗时3m，系统用时10m，此处可设置线程，但salmon仍调用全部资源
    # 需要指定输出文件夹，包括4.3中的参数的输出目录
    nohup metawrap quant_bins -b temp/bin_refinement/metawrap_50_10_bins -t 8 \
      -o temp/bin_quant -a temp/megahit/final.contigs.fa temp/qc/ERR*.fastq  &bg
    # 文件名字改变
    # 结果包括bin丰度热图`temp/bin_quant/bin_abundance_heatmap.png`
    # 如果想自己画图，原始数据位于`temp/bin_quant/bin_abundance_table.tab`
    ls -l temp/bin_quant/bin_abundance_heatmap.png

## 4.4 Bin注释

    # Taxator-tk对每条contig物种注释，再估计bin整体的物种，11m (用时66 min)
    nohup metawrap classify_bins -b temp/bin_refinement/metawrap_50_10_bins \
      -o temp/bin_classify -t 20 &
      
    # The file /disk2/home/meta_ct/NCBI_NT_DB/nt.00.nhd doesnt exist, which likely means that           *****
    # *****          you havent set the correct path to your NCBI_nt database or	 havent downloaded it. Please           *****
    # *****            look and the the script and manually set the correct path to NCBI_nt, and follow these            *****
    # *****           steps to download the	  database: wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz;            *****
    # *****                                   for a in nt.*.tar.gz; do tar xzf ; done   
    # which config-metawrap
    # 注释结果见`temp/classify_bins/bin_taxonomy.tab`
    
    # 基于prokka基因注释，4m
    metaWRAP annotate_bins -o temp/bin_annotate -t 2 \
      -b temp/bin_refinement/metawrap_50_10_bins
    # 每个bin基因注释的gff文件bin_funct_annotations, 
    # 核酸ffn文件bin_untranslated_genes，
    # 蛋白faa文件bin_translated_genes


## S4附录. 其它可用分析

### S4.1 物种组成预测(kraken+krona，可选)

    # kraken物种注释, -o输出目录, -t线程数, -s抽样1M减少计算量, 质控数据, contig文件
    # 需要有256GB以上可用内存，7G演示数据，8线程耗时15m，测试服务器暂不开放，防多人运行内存耗尽死机
    # metawrap kraken -o temp/kraken -t 8 -s 1000000 \
    #   temp/final.contigs.fa temp/ERR*.fastq 
    # 结果文件夹中有注释结果文件*.kraken(每个reads的注释结果)、*.krona(注释结果分类汇总)
    # 可视化的Krona网页图表kronagram.html，链接结果到result方便查看
    # mkdir -p result
    # ln temp/kraken/kronagram.html result/


### S4.2 Bin可视化

    # 计算每个contig的GC含量和在每个样本中的丰度
    # 8线程耗时1-6h
    nohup metawrap blobology -a temp/megahit/final.contigs.fa -t 30 \
      -o temp/bloblogy --bins temp/bin_refinement/metawrap_50_10_bins \
      temp/qc/ERR*.fastq &
    # Something went wrong with making the blob file from 
    # the .bam and .megablast files；可能R包环境影响，但可视化不影响主流程
    # 结果为final.contigs.binned.blobplot，方便使用ggplot2可视化
    # 参考脚本${soft}/envs/metawrap/bin/metawrap-scripts/blobology/makeblobplot_with_bins.R


### S4.3 重组装

    # 需合并所有样本作为此步输入, 6s
    cat temp/qc/ERR*_1.fastq > temp/qc/all_1.fq
    cat temp/qc/ERR*_2.fastq > temp/qc/all_2.fq
    # 提纯的bin还可以通过再组装进一步改善结果，8核，100G内存，用时2小时
    nohup metawrap reassemble_bins -o temp/bin_reassemble \
      -1 temp/qc/all_1.fq -2 temp/qc/all_2.fq -t 10 -m 5 \
      -c 50 -x 10 -b temp/bin_refinement/metawrap_50_10_bins &
    # 结果统计见`temp/bin_reassemble/reassembled_bins.stats`，
    # `temp/bin_reassemble/reassembly_results.png`，
    # 比对重组装前后的变化，N50、这完整度和污染率均有改进。
    # `temp/bin_reassemble/reassembled_bins.png`展示CheckM对bin评估结果的可视化。
    
    
    # 其它单菌基因组分析可结合进化、COG、KEGG、GO、碳水化合物、
    # 抗生素抗性和基因簇分析(如antismash)等
    
    # 分析结束，退出虚拟环境
    conda deactivate


# 五、进化树构建和美化

    # 进入工作目录
    mkdir -p ~/meta/35Evolution
    cd ~/meta/35Evolution
    # 右侧进入工作目录，上传U盘35Evolution目录中测试数据input.zip，包括9个细菌的蛋白组
    
    # 手动加载指定位置的环境变量，如使用某用户装好的conda环境
    source deactivate
    export PATH="/home/liuyongxin/miniconda2/bin:$PATH"


## 5.1 OrthoFinder同源基因的筛选(Linux)

    # 运行Orthofinder
    # -M基因树推断方法，-A多序列比对方法，-S搜索方法默认diamond
    # -a并行任务数，-t线程数
    # -f输入目录，里面包含每个物种蛋白序列fasta文件，可上传input.ZIP测试文件夹
    # 运行orthofinder查看帮助，详见官方文档 
    # https://github.com/davidemms/OrthoFinder/blob/master/OrthoFinder-manual.pdf
    nohup time orthofinder -M msa -A mafft -S diamond -t 8 -a 8 -f input &bg
    # 9个细菌64核运行大约3小时，结果目录位于输入文件夹中OrthoFinder目录
    
    # 参数详细说明：OPTIONS:
    # -t <int>: 序列搜索调用线程数，整数，默认40，根据服务器配置和使用人数调整，推荐8 - 12
    # -a <int>: 序列分析并行任务数，整数，默认1，根据服务器配置和使用人数调整，推荐 3 -5
    # -M <txt>: 基因树推断方法，文本，默认dendroblast，可选msa
    # -S <txt>: 指定比对软件，文本，默认 : diamond，可选blast, mmseqs, blast_gz
    # -A <txt>: 多序列比对程序，文本，需要有'-M msa'参数，默认为mafft，可选muscle, mafft
    # -T <txt>: 树推断方法，文本，需要有'-M msa'参数，默认为fasttree，可选iqtree, raxml-ng, raxml
    # -s <file>: 用户自定义根物种树，文件
    # -I <int>: MCL膨胀参数，默认1.5
    # -x <file>: 输出OrthoXML格式结果
    # -p <dir>: 保存临时文件目录
    # -1: 仅进行单种方法的序列搜索
    # -X: 序列ID中不添加序列名字
    # -n <txt>: 输出目录后辍
    # -o <txt>: 设置输出文件夹，默认为输入文件夹中的OrthoFinder目录
    # -h: 输出帮助文档
    # 工作流程中断参数：WORKFLOW STOPPING OPTIONS:
    # -op: 运行准备BLAST输入文件后停止；Stop after preparing input files for BLAST
    # -og: 推断同源组后停止；Stop after inferring orthogroups
    # -os: 输出同源组序列后停止；Stop after writing sequence files for orthogroups (requires '-M msa')
    # -oa: 多序列比对后停止；Stop after inferring alignments for orthogroups (requires '-M msa')
    # -ot: 建树后停止；Stop after inferring gene trees for orthogroups
    # 工作流程继续参数：WORKFLOW RESTART COMMANDS:
    # -b  <dir>: 基于预计算的BLAST结果目录继续运行; Start OrthoFinder from pre-computed BLAST results in <dir>
    # -fg <dir>: 基于预计算的orthogroups结果目录继续运行; tart OrthoFinder from pre-computed orthogroups in <dir>
    # -ft <dir>: 基于预计算的基因树结果目录继续运行; Start OrthoFinder from pre-computed gene trees in <dir>    
    
    # -f : 输入文件工作目录
    # 需将待分析序列提前放入此目录
    # 若基因组较大，使用 diamond 代替 blast 作为比对软件，可大幅提高执行速度。
    # orthofinder 会自动生成orthomcl 格式的结果，除此之外还会生成单拷贝基因文件及其树文件。
    # 过程的程序选择比orthomcl, 结果更直接。


## 5.2 IQ-TREE快速构建ML进化树(Linux/Widows)

    # IQ-TREE识别的格式包括Phylip; fasta; nexus以及clustalw
    # 该方法适合于大数据，例如几百个OTUs、多基因的系统发育树！
    # IQ-TREE整合了ModelFinder这个模块，可以快速获得目的数据的最佳进化模型
    # 大家把iqtree.exe文件和libiomp5md.dll文件拷贝到自己电脑C盘的bin目录下，
    # 前面我们已经把这个bin目录整体放在了环境变量中，这样我们就可以直接调用该文件夹下的所有程序。
    
    # iqtree会生成大量文件
    # 新建目录并进入
    mkdir -p iqtree && cd iqtree
    # 准备输入文件，注意与上步输入目录，结果生成日期对应，需要手动修改
    ln -s ../input/OrthoFinder/Results_Jun23/MultipleSequenceAlignments/SpeciesTreeAlignment.fa ./
    
    # 用IQ-TREE进行模型的筛选
    # -s 参数对应的是序列文件；
    # -m 参数默认使用ModelFinder选模型，也可使用jModeltest或ProtTest，只需可将MF改为TEST即可；
    # -nt 对应的是CPU线程数，根据硬件配置修改（如：4线程，将AUTO改为4即可），也可以使用AUTO交给程序来选择 ....
    time iqtree -s SpeciesTreeAlignment.fa -m MF -nt AUTO
    # 运行48m，Best-fit model: LG+F+R3 chosen according to BIC，详见SpeciesTreeAlignment.fa.log
    # ModelFinder会根据不同标准推荐最佳模型，而ModelFinder默认推荐使用BIC标准。
    
    # ML系统发育树的构建
    # iqtree -h 查看帮助文档, -o 指定外群，-pre输出文件前缀
    time iqtree -s SpeciesTreeAlignment.fa -m LG+F+R3 -bb 1000 -redo -alrt 1000 -nt AUTO -pre output_Pseudomonas
    # 4m, *.contree 为一致树文件，日志见 output_Pseudomonas.log


## 5.3 进化树美化

    cd ${wd}/35Evolution/iTOL
    # 访问 http://itol.embl.de/ ，上传tree_3837_genomes/otus.nwk，
    # 再拖拽下方生成的注释方案于树上即美化
    
    # table2itol生成注释文件
    # 支持绘制iTOL domains, 中的colour strips, simple bars, 
    # gradients, binary data, heat maps以及texts.
    # 部分支持iTOL中的分支注释
    # 选择合适的数据输入格式
    # 提供精心选择的颜色向量多达40个级别和可选的将它们与符号以最大化对比度。
    
    # 方案A. 外圈颜色、形状分类和丰度方案
    # annotation.txt 菌对应物种注释、人工分组和丰度等
    # -a 找不到输入列将终止运行（默认不执行）-c 将整数列转换为factor或具有小数点的数字
    # -t 偏离提示标签时转换ID列，-w 颜色带，区域宽度等， -D输出目录，
    # -i taxon_oid列名，-l Name显示名称如种/属/科名，
    # 注意：指定的列名务必与注释文章中一致
    Rscript table2itol/table2itol.R -a -c double -D planA \
      -i taxon_oid -l Name -t %s -w 0.5 annotation.txt
    # 生成注释文件中每列为单独一个文件
    
    ## 方案B. 生成丰度柱形图注释文件
    Rscript table2itol/table2itol.R -a -d -c none \
      -D planB -b IMG_Phylum -i taxon_oid -l Name -t %s -w 0.5 annotation.txt
    
    ## 方案C. 自定义颜色
    ${Rscript} -- table2itol/table2itol.R -a \
      -C table2itol/tests/INPUT/colours_1.yml -c double -D planC \
      -i taxon_oid -l IMG_Genus -t %s -w 0.5 annotation.txt


## S5.1 进化树相关软件安装

###  orthofinder及依赖软件

    # 方法1. 推荐conda安装orthofinder及依赖关系，16.3M，-y自动安装
    conda install -y orthofinder
    orthofinder # version 2.3.3
    
    # 方法2. 手动安装    
    # OrthoFinder软件的下载地址 https://github.com/davidemms/OrthoFinder
    # 下载完成后不需要安装，只需要指定路径或者添加到环境变量就行，
    # 但需要手动安装依赖关系：MCL, FastME, diamond/MMseqs2/Blast+
    # 参考博文：https://www.jianshu.com/p/296186cde6dc    
    ## 依赖软件MCL软件的下载和安装
    wget https://www.micans.org/mcl/src/mcl-latest.tar.gz --no-check-certificate
    wget https://micans.org/mcl/src/mcl-latest.tar.gz
    tar -zxvf mcl-latest.tar.gz
    cd mcl-14-137 #（视具体情况而定）
    ./configure
    make
    make check
    sudo make install
    ## 依赖软件diamond软件的下载和安装
    #没有root权限的用户把diamond所在目录加入环境变量。
    ##关于几款比对软件Blast，Blast+，Diamond的介绍：https://www.cnblogs.com/jessepeng/p/10447373.html
    wget https://github.com/bbuchfink/diamond/releases/download/v0.9.24/diamond-linux64.tar.gz
    tar -xzvf diamond-linux64.tar.gz
    sudo cp diamond /usr/local/bin
    ## 依赖软件FastME软件的下载和安装
    下载，解压并将主程序加入环境变量：
    wget http://www.atgc-montpellier.fr/download/sources/fastme/fastme-2.1.5.tar.gz
    tar -zxvf fastme-2.1.5.tar.gz
    sudo cp fastme-2.1.5/binaries/fastme-2.1.5-linux64 /usr/local/bin/fastme
    ## 依赖软件MMseqs2软件的下载和安装https://github.com/soedinglab/MMseqs2/releases
    #下载对应版本，解压并将主程序拷贝至存在于环境变量的目录下或将其所在的目录加入环境变量：
    wget https://github.com/soedinglab/MMseqs2/releases/download/7-4e23d/MMseqs2-Linux-AVX2.tar.gz
    tar -zxvf MMseqs2-Linux-AVX2.tar.gz
    sudo cp mmseqs2/bin/mmseqs /usr/local/bin
