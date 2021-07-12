###################################################
# subject:NGS后数据整理
# email:jingzhy@shanghaitech.edu.cn
# author:荆征宇
# time:2020/08/15
###################################################

#读取并计算污染情况--------------------------------------------------------------
contamination_condition <- function(data_ori,cutoff = 99){
  marker <- data_ori$Name[which(str_detect(data_ori$Name,"Map"))-1]
  contamination <- data.frame(index = marker,ratio =
                                data_ori$Name[which(str_detect(data_ori$Name,"Map"))])
  contamination$ratio <- str_remove(contamination$ratio,"Map.*= ")

  output <- contamination %>% filter(as.numeric(str_extract(ratio,"\\d*")) < cutoff)
  return(output)
}


#对数据中的样品进行分组，加上编号--------------------------------------------------------------
marc_group <- function(data_ori){   #该函数给文件中各个结果添加上对应的序号
  for_out <- data_ori
  for_out$marker <- NA
  for_out$order <- NA
  group_num <- NA

  for (index in 1:length(for_out$Name)) {   #ifelse 可以对向量每个数都进行判断。if 只能对一个值进行判断
    ifelse(str_detect(for_out[index,1],"#"),
           group_num <- str_extract(for_out[index,1],"\\d+"),
           ifelse(!str_detect(for_out[index,1],"(Name)|(Mapping)"),
                  for_out$marker[[index]] <- group_num,
                  for_out$marker[[index]] <- NA_real_
           )
    )
  }
  #做一些规整，调整一些数据类型
  for_out %<>% filter(!is.na(marker))
  for_out$NumReads <- round(as.numeric(for_out$NumReads))
  for_out$TPM <- round(as.numeric(for_out$TPM),1)
  for_out <- as_tibble(for_out)
  for_out$Name <- factor(for_out$Name,levels = unique(for_out$Name))
  for_out$mouse_strain <- substring(for_out$Name,1,4)
  #添加gRNA顺序，自动识别长度
  for(line in unique(for_out$mouse_strain)){
    for_out$order[for_out$mouse_strain == line] <-
      1:length(unique(for_out$Name[for_out$mouse_strain == line]))
  }

  for_out <- for_out %>% dplyr::select(-Length, -EffectiveLength)

  #为了解决salmon将reads平分给相同的序列，来将1329中相同序列的gRNA的read数×2
  for_out$NumReads[str_detect(for_out$Name,'(^1329_kpnb1-1$)|(^1329_kpnb1-22$)')] <-
    for_out$NumReads[str_detect(for_out$Name,'(^1329_kpnb1-1$)|(^1329_kpnb1-22$)')] * 2

  #两者分别*2最后会差别1 2个read

  return(for_out)
}


#看每个sgRNA在总体中所占的丰度比例------------------------------------------------------
percent_marc <- function(data_group, remove_grna_index = c(1,3), remove_gene = "Cd47"){
  new_data_f <- data_group[0,]
  new_data_f$percent <- 1
  new_data_f$percent_remove_some_grna <- 1

  for(i in unique(data_group$marker)){
    single_data_f <- data_group %>% filter(marker == i)

    #计算grna所占百分比
    single_data_f$percent <- single_data_f$NumReads / sum(single_data_f$NumReads)
        #排除掉两种不要的gRNA
    s <- sum(single_data_f$NumReads[-c(remove_grna_index,
                                    which(
                                      str_detect(single_data_f$Name, regex(remove_gene, ignore_case = T))
                                      )
                                    )])

    single_data_f$percent_remove_some_grna <- single_data_f$NumReads / s
    single_data_f$percent_remove_some_grna[remove_grna_index] <- NA
    single_data_f$percent_remove_some_grna[str_detect(single_data_f$Name, remove_gene)] <- NA

    #将每个样品重新合并
    new_data_f <- rbind(new_data_f, single_data_f)
  }
  new_data_f <- new_data_f %>% dplyr::select(-TPM)
#  new_data_f$logFC <- NA
  new_data_f$marker <- factor(new_data_f$marker, levels = unique(new_data_f$marker))
  new_data_f$median_reads <- NA_real_
  new_data_f$median_FC <- NA_real_
  new_data_f$control_reads <- NA_real_
  new_data_f$besides_ratio <- NA_real_
  new_data_f$ratio_FC <- NA_real_
  new_data_f$percent_FC <- NA_real_
  new_data_f$control_percent <- NA_real_

  return(new_data_f)
}


#差异gRNA画图前数据整理----------------------------------------------------------------
make_up_number <- function(gene_list){
  data <- data.frame(gene = gene_list, number = 1)
  gene_frame <- data.frame(gene = unique(gene_list), number = 1)
  collection <- c()
  for(i in 1:length(gene_list)){
    if(gene_list[[i]] %in% collection){
      gene_frame[gene_frame[[1]] == gene_list[[i]], 2] <-
        gene_frame[gene_frame[[1]] == gene_list[[i]], 2] + 1
    }
    data[i, 2] <- gene_frame[gene_frame[[1]] == gene_list[[i]], 2]
    collection <- c(collection, gene_list[[i]])
  }
  return(letters[data[[2]]])
}

marc_DataClean_for_plot <- function(data_nor_logfc){
  data_nor_logfc$Name <- str_remove(data_nor_logfc$Name,"\\d+_")
  data_nor_logfc$Name <- factor(data_nor_logfc$Name, levels = unique(data_nor_logfc$Name))
  data_nor_logfc$serial_number <- NA
  data_nor_logfc$gene <- NA


  #将1329和100-mer分开处理
  sub_1329 <- data_nor_logfc %>% filter(mouse_strain == 1329)
  sub_other <- data_nor_logfc %>% filter(mouse_strain != 1329)

  #先对1329进行处理
  if(length(sub_1329$Name) != 0){
    sub_1329$type <- case_when(
      str_detect(sub_1329$Name,"pol") | str_detect(sub_1329$Name,"kpnb") ~ "down",
      str_detect(sub_1329$Name,"p53") | str_detect(sub_1329$Name,"pten") ~ "up",
      str_detect(sub_1329$Name,"cd19") | str_detect(sub_1329$Name,"cd45") ~ "control",
      str_detect(sub_1329$Name,"cd47") ~ "first gRNA")

    grna_table <- data.frame(name = as.character(unique(sub_1329$Name)))
    sub_1329$gene <- as.data.frame(str_split(grna_table$name,"-",simplify = TRUE))[[1]] %>%
      as.character() %>% rep(length(unique(sub_1329$marker)))

    #对gRNA加上编号
    for(gene in unique(sub_1329$gene)){
        sub_1329$serial_number <- rep(make_up_number(filter(sub_1329, marker == marker[[1]])$gene),
                                      length(unique(sub_1329$marker)))

      #对kpnb1做特殊处理，将两个一样的gRNA改成一样的名字
      if(gene == "kpnb1"){
        sub_1329$serial_number[sub_1329$serial_number == "c" & sub_1329$gene == "kpnb1"] <- "a"
      }
    }
  }

  if(length(sub_other$Name) != 0){
    #对100-mer进行处理
    sub_other$type <- case_when(
      sub_other$order == 1 ~ "first gRNA",
      str_detect(sub_other$Name,"NC|Nc|GFP") ~ "control",
      TRUE ~ "unknown"
    )
    sub_other$gene <- sub_other$Name

    sub_other$gene <- str_remove(sub_other$gene, "-\\d+")
    sub_other$gene <- sub_other$Name %>% as.character()
    sub_other$gene[str_detect(sub_other$Name,"NC5|Nc5")] <- "NC"
    sub_other$gene[str_detect(sub_other$Name,"GFP")] <- "GFP"

    for(i in unique(sub_other$mouse_strain)){
      sub_other$serial_number[sub_other$mouse_strain == i] <-
        rep(make_up_number(sub_other %>%
                             filter(mouse_strain == i) %>%
                             filter(marker == marker[[1]]) %>%
                             dplyr::select(gene) %>%
                             unlist()),
            length(unique(sub_other$marker[sub_other$mouse_strain == i])))
    }

  }
  output <- rbind(sub_1329,sub_other) %>% arrange(marker, mouse_strain)

  return(output)
}

#添加样品信息----------------------------------------------------------------------------------------------
data_clean_for_analysis <- function(data_for_plot, sample_path = "./sample_info.xlsx"){

  new_data_f <- data_for_plot
  sample_info <- xlsx::read.xlsx(sample_path, sheetIndex = 1, encoding = "UTF-8")
  new_data_f$sample <- sample_info$sample[match(paste0("BZ",new_data_f$marker), sample_info$index_num)]
  #new_data_f$sample <- str_replace(new_data_f$sample,"\\s+","_")
  new_data_f$organ <- str_remove(new_data_f$sample,".+?-")
  new_data_f$mouse_index <- str_extract(new_data_f$sample,"^[^-]+")
  new_data_f$mouse_index <- factor(new_data_f$mouse_index, levels = unique(new_data_f$mouse_index))

  return(new_data_f)
}

#在后期使用双阳来给样品进行计算对照--------------------------------------------------------------------------
gm_mean = function(x, na.rm = TRUE, zero_propagate = T){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero_propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  }else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

find_nearest_sample <- function(control_data, strain, first_nc_index, sample_ratio, y){
  #计算重组程度最相近的几个数据
  control_first <-
    control_data %>%
    filter(mouse_strain == strain, !is.na(median_reads)) %>%
    group_by(sample, unique_num) %>%
    summarise(first = percent[1],
              n1 = percent[1] / percent[first_nc_index]) %>%
    mutate(diff = abs((sample_ratio / n1) - 1)) %>%  #使用样本的比值来进行标准化使其接近1，再减一后取绝对值拿到最近似0的
    dplyr::select(sample, first, unique_num, n1, diff) %>%
    arrange(diff)

  if(length(control_first[[1]]) == 0){ #去除没有对照的品系剩余计算步骤
    return(NULL)
  }
  six_nearest_data <- control_first$unique_num[1:6]
  within_ten <- control_first$unique_num[control_first$diff < 0.2]
  #比较这两种哪个多就使用哪些
  if(length(six_nearest_data) >= length(within_ten)){
    control_index <- six_nearest_data
  }else{
    control_index <- within_ten
  }
  control_read <-
    control_data %>%
    filter(unique_num %in% control_index) %>%
    group_by(Name, order) %>%
    summarise(value = mean(get(y)),
              percentage = mean(percent)) %>%
    arrange(order)

  return(control_read)
}

normalization_with_median <-
  function(data, control_l = T, zero = F){

    new_data <- data[0,]
    for_out <- data[0,]
    #计算出中位数标准化用到的的那个平均数
    for(line in unique(data$mouse_strain)){
      assign(paste0("median",line),
             data %>% filter(mouse_strain == line) %>%
               group_by(Name) %>% summarise(mean = gm_mean(NumReads, zero_propagate = zero)))
    }

    #开始循环对每个样本进行处理
    for(i in unique(data$unique_num)){

      single_data <- data %>% filter(unique_num == i)
      for_median_data <- get(paste0("median", single_data$mouse_strain[[1]]))  #获得该样本对应的每个gRNA的平均值

      #对每个样品计算标准化因子
      if(control_l){  #判断是否只使用control来进行中位数计算
        #提取用来计算的reads数
        if(single_data$mouse_strain[[1]] == 1329){          #对于1329只使用cd19
          for_median <- unlist(single_data %>% filter(gene == "Cd19") %>% dplyr::select(NumReads))
          names(for_median) <- unlist(single_data %>% filter(gene == "Cd19") %>% dplyr::select(Name))
        }else{
          for_median <- unlist(single_data %>% filter(type == "control") %>% dplyr::select(NumReads))
          names(for_median) <- unlist(single_data %>% filter(type == "control") %>% dplyr::select(Name))
        }
      }else{
        for_median <- unlist(single_data %>% filter(order != 1) %>% dplyr::select(NumReads))
        names(for_median) <- unlist(single_data %>% filter(order != 1) %>% dplyr::select(Name))
      }
      #计算中位数的那个size factor
      size_factor <- median(for_median / for_median_data$mean[match(names(for_median),for_median_data$Name)])
      #进行标准化
      if(size_factor > 0.01){
        single_data$median_reads <-
          as.numeric(unlist(
            (single_data %>% dplyr::select(NumReads)) /
              (size_factor)
          ))

        single_data$median_reads <- round(single_data$median_reads)
      }

      new_data <- rbind(new_data, single_data)
    }

    #计算差异倍数
    control_data <- filter(new_data, cas == "no", !str_detect(organ,"thymus|heart."))
    #对每个数据进行计算
    for(i in unique(data$unique_num)){
      single_data <- new_data %>% filter(unique_num == i)

      #抽出作为对照的数据
      sample_ratio <- single_data$percent[1] / single_data$percent[single_data$type == "control"][[1]]
      control_read <- find_nearest_sample(control_data, single_data$mouse_strain[1],
                                          first_nc_index = single_data$order[single_data$type == "control"][[1]],
                                          sample_ratio = sample_ratio, "median_reads")

      if(is.null(control_read)){
        for_out <- rbind(for_out, single_data)
        next
      }

      #用标准化后的值计算差异倍数
      single_data$median_FC <-
        single_data$median_reads /
        control_read$value
      #增加一列直接使用百分比来除
      single_data$percent_FC <-
        single_data$percent /
        control_read$percentage

      single_data$control_reads <- control_read$value
      single_data$control_percent <- control_read$percentage

      for_out <- rbind(for_out, single_data)
    }
    return(for_out)
  }


find_flanking_gRNA <- function(site, all_num, flank_num = 4, last = F){
  if(!last){
    all_num <- all_num[c(-1,-length(all_num))]
  }
  distance <- data.frame(order = all_num, dis = abs(site - all_num))
  distance <- distance[order(distance$dis),]
  distance <- distance[distance$dis != 0,]

  return(distance$order[1:flank_num])
}

normalization_with_beside <-
  function(data, flank_num = 4){

    ratio_data <- data[0,]
    for_out <-  data[0,]

    for(i in unique(data$unique_num)){  #对每个样品挨个进行处理
      single_data <- data %>% filter(unique_num == i)
      #判断最后一个gRNA是否参与计算
      if(unique(single_data$mouse_strain) %in% c(1287,1329)){
        last <- T
      }else{
        last <- F
      }
      #对所有的gRNA进行计算
      for(n in 1:length(single_data$Name)){
        single_data$besides_ratio[n] <- (single_data$NumReads[n] + 1e-20) /
          median(single_data$NumReads[find_flanking_gRNA(n, single_data$order, flank_num, last)] + 1e-20)

      }
      ratio_data <- rbind(ratio_data, single_data)
    }

    #计算差异倍数
    #抽出作为对照的数据
    control_data <- filter(ratio_data, cas == "no", !str_detect(organ,"thymus|heart."))
    #对每个数据进行计算
    for(i in unique(data$unique_num)){
      single_data <- ratio_data %>% filter(unique_num == i)
      #抽出作为对照的数据
      sample_ratio <- single_data$percent[1] / single_data$percent[single_data$type == "control"][[1]]

      control_ratio <- find_nearest_sample(control_data, single_data$mouse_strain[1],
                                          first_nc_index = single_data$order[single_data$type == "control"][[1]],
                                          sample_ratio = sample_ratio, "besides_ratio")

      if(is.null(control_ratio)){
        for_out <- rbind(for_out, single_data)
        next
      }
      #用标准化后的值计算差异倍数
      single_data$ratio_FC <-
        single_data$besides_ratio /
        control_ratio$value

      for_out <- rbind(for_out, single_data)
    }
    return(for_out)
  }


#合并所有的数据,并输出到包外----------------------------------------------------------------------------------------
#' Title
#'
#' @param all_data
#'
#' @return
#' @export
#'
#' @examples
normalization <- function(all_data){
  timestart <- Sys.time()
  all_data <- normalization_with_median(all_data)
  all_data <- normalization_with_beside(all_data)

  timeend <- Sys.time()
  runningtime <- timeend-timestart
  print(runningtime)
  return(all_data)
}

#' Title
#'
#' @param all_data
#' @param control
#'
#' @return
#' @export
#'
#' @examples
formalize_gene_name <- function(all_data, control){
  control_data <- data.frame(mouse_index = str_extract(control, "\\w?\\d*"),
                             mouse_strain = str_extract(control, "\\d+$"))
  #添加唯一识别符号
  times <- all_data %>% mutate(temp = paste0(all_data$marker,all_data$batch)) %>% group_by(temp) %>% dplyr::count(temp)
  all_data$unique_num <- rep(1:length(unique(paste0(all_data$marker,all_data$batch))),
                             times = times$n[match(unique(paste0(all_data$marker,all_data$batch)), times$temp)])
  #对对照添加识别信息
  all_data$cas <- NA
  for(strain in unique(control_data$mouse_strain)){
    sub_control <- control_data$mouse_index[control_data$mouse_strain == strain]
    all_data$cas[all_data$mouse_strain == strain] <- case_when(
      all_data$mouse_index[all_data$mouse_strain == strain] %in% sub_control ~ "no",
      T ~ "yes"
    )
  }


  all_data$Name <- str_replace(all_data$Name, "(?i)^pol2", "Polr2a")
  all_data$gene <- str_replace(all_data$gene, "(?i)^pol2", "Polr2a")
  all_data$Name <- str_replace(all_data$Name, "(?i)^p53", "Tp53")
  all_data$gene <- str_replace(all_data$gene, "(?i)^p53", "Tp53")

  all_data$Name <- paste(toupper(substr(all_data$Name, 1, 1)), substr(all_data$Name, 2, 50), sep="")
  all_data$gene <- paste(toupper(substr(all_data$gene, 1, 1)), substr(all_data$gene, 2, 50), sep="")

  all_data$organ <- str_replace(all_data$organ, "naïve", "naive")
  all_data$sample <- str_replace(all_data$sample, "naïve", "naive")

  all_data$Name <- factor(all_data$Name, levels = unique(all_data$Name))
  return(all_data)
}

#' Title
#'
#' @param name
#' @param sample_path
#'
#' @return
#' @export
#'
#' @examples
data_clean_for_ngs <- function(name, sample_path = "./sample_info.xlsx"){
  data_ori <- read.table(name, fill = T, sep = "\t", comment.char = "", stringsAsFactors = F) %>%
    as_tibble()
  colnames(data_ori) <- data_ori[3,]

  print(contamination_condition(data_ori, 85))

  data_group <- marc_group(data_ori)
  data_per <- percent_marc(data_group, remove_grna_index = c(1:3), remove_gene = "cd47")
  #data_nor_logfc <- marc_normalization(data_per, gm = T, control_l = T, two_sides = T,
  #control_index = control_index)

  data_for_plot <- marc_DataClean_for_plot(data_per)
  data_for_analysis <- data_clean_for_analysis(data_for_plot, sample_path = sample_path)
  data_for_analysis$batch <- str_extract(name,"HP.*(?=\\.)")
  return(data_for_analysis)
}




