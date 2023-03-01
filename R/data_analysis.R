#指定对照和样本，来计算特定样本对中的差异------------------------------------------------------------------
#提取用于批次比较的样本
extract_samples <- function(data_for_analysis, samples_info = NULL, unique_num_info = NULL){
  out_data <- data_for_analysis[0,]
  out_data$numerator  <- NA_character_
  out_data$group <- NA_character_

  if(!is.null(samples_info)){
    unique_num_info <- data.frame(exp_index = 0, con_index = 0, group = 0)

    for (i in 1:dim(samples_info)[[1]]) {
      unique_num_info[i, 1] <-
        data_for_analysis %>%
        dplyr::filter(batch == samples_info[i, 1], marker == samples_info[i, 2], order == 1) %>%
        select(unique_num)

      unique_num_info[i,2] <-
        data_for_analysis %>%
        dplyr::filter(batch == samples_info[i, 3], marker == samples_info[i, 4], order == 1) %>%
        select(unique_num)

      unique_num_info[i, 3] <- samples_info[i, 5]
    }
    unique_num_info <- unique_num_info[order(unique_num_info$group),]
  }

  for(i in 1:dim(unique_num_info)[[1]]){
    #将等下分析作为分子的样本放进去
    temp_df <-
      data_for_analysis %>%
      dplyr::filter(unique_num == unique_num_info[i, 1]) %>%
      mutate(group = unique_num_info[i, 3], numerator = "yes")
    out_data <- rbind(out_data, temp_df)
    #将等下分析作为分母的样本放进去
    temp_df <-
      data_for_analysis %>%
      dplyr::filter(unique_num == unique_num_info[i, 2]) %>%
      mutate(group = unique_num_info[i, 3], numerator = "no")
    out_data <- rbind(out_data, temp_df)
  }
  return(out_data)
}


#获取用于两两比较的数据的unique_num
find_sample_by_regex <-
  function(single_mouse_data, exp_regex, con_regex, group, limit = NULL, except = NULL){
    if (!is.null(limit)) {
      single_mouse_data <- single_mouse_data %>% dplyr::filter(str_detect(sample, limit))
    }
    if (!is.null(except)) {
      single_mouse_data <- single_mouse_data %>% dplyr::filter(!str_detect(sample, except))
    }

    output <- data.frame(exp_index = 0, con_index = 0, group = 0)
    exp <- unique(single_mouse_data$unique_num[str_detect(single_mouse_data$sample, exp_regex)])
    con <- unique(single_mouse_data$unique_num[str_detect(single_mouse_data$sample, con_regex)])

    #防止同一个老鼠有多组样本的情况
    for(i in 1:length(exp)){
      output[i, 1] <- exp[i]
      if(length(con) == 1){
        output[i, 2] <- con[1]
      }else{
        output[i, 2] <- con[i]
      }
      output[i, 3] <- group
    }
    return(output)
  }

#' Title
#'
#' @param data_for_analysis
#' @param mouse_index_list
#' @param strain
#'
#' @return
#' @export
#'
#' @examples
get_unique_num_table <-
  function(data_for_analysis, mouse_index_list, strain, search_table = NULL){
    unique_num_output <- data.frame(exp_index = 0, con_index = 0, group = 0)[0,]

    if(is.null(search_table)){
      search_table <-
        tibble::tribble(
          ~first_sample, ~second_sample, ~lable, ~limit, ~except,
          #肿瘤中cd4 pd1阳性与阴性比较
          "(?i)cd4.*pd-?1_(p|hi)", "(?i)cd4.*pd-?1_(n|lo)", "cd4_pd-1P-N", "(?i)tumor", NULL,
          #肿瘤中cd8 pd1阳性与阴性比较
          "(?i)cd8.*pd-?1_(p|hi)", "(?i)cd8.*pd-?1_(n|lo)", "cd8_pd-1P-N", "(?i)tumor", NULL,
          #肿瘤中cd4 pd1阳性与外周naive比较
          "(?i)(?<!spleen.{0,200}?)cd4.*pd-?1_(p|hi)(?!.{0,200}?spleen)",
          "(?i)cd4.*naive", "cd4_pd-1P-naive", NULL, NULL,
          #肿瘤中cd8 pd1阳性与外周naive比较
          "(?i)(?<!spleen.{0,200}?)cd8.*pd-?1_(p|hi)(?!.{0,200}?spleen)",
          "(?i)cd8.*naive", "cd8_pd-1P-naive", NULL, NULL,
          #cd8细胞的TNF等分泌物
          "(?i)tnf_(p|hi)", "(?i)tnf_(n|lo)", "cd8-TNF_P-N", "(?i)cd8", NULL,
          #cd8细胞的107a等分泌物
          "(?i)cd107a_(p|hi)", "(?i)cd107a_(n|lo)", "cd8-cd107a_P-N", "(?i)cd8", NULL,
          #肿瘤中cd4阳性与外周naive比较
          "(?i)(?<=tumor.{0,200}?)cd4|cd4(?=.{0,200}?(tumor))",
          "(?i)(?<!tumor.{0,200}?)cd4.*naive(?!.{0,200}?tumor)",
          "cd4_tumor-naive", NULL, "(?i)pd-?1",
          #肿瘤中cd8阳性与外周naive比较
          "(?i)(?<=tumor.{0,200}?)cd8|cd8(?=.{0,200}?(tumor))",
          "(?i)(?<!tumor.{0,200}?)cd8.*naive(?!.{0,200}?tumor)",
          "cd8_tumor-naive", NULL, "(?i)pd-?1",
          #外周cd4 effector与naive比较
          "(?i)cd4.*(?<!cen.{0,10})effector", "(?i)cd4.*naive",
          "cd4_effector-naive", "(?i)spleen|ln", NULL,
          #外周cd8 effector与naive比较
          "(?i)cd8.*(?<!cen.{0,10})effector", "(?i)cd8.*naive", "cd8_effector-naive", "(?i)spleen", NULL,
          #外周cd8 CM与naive比较
          "(?i)cd8.*(CM|central)", "(?i)cd8.*naive", "cd8_CM-naive", "(?i)spleen", NULL,
          #巨噬细胞里tnf高表达与低表达之间比较
          "(?i)(?<!spleen.{0,200}?)tnf_(p|hi)(?!.{0,200}?spleen)",
          "(?i)tnf_(n|lo)", "MAC_TNF_P-N", "(?i)(F4/80_P|cd11b_P)", "(?i)spleen",
          #巨噬细胞里INOS高表达与低表达之间比较
          "(?i)(?<!spleen.{0,200}?)inos_(p|hi)(?!.{0,200}?spleen)",
          "(?i)inos_(n|lo)", "MAC_INOS_P-N", "(?i)(F4/80_P|cd11b_P)", "(?i)spleen",
          #巨噬细胞里arg高表达与低表达之间比较
          "(?i)(?<!spleen.{0,200}?)arg1?_(p|hi)(?!.{0,200}?spleen)",
          "(?i)arg1?_(n|lo)", "MAC_ARG1_P-N", "(?i)(F4/80_P|cd11b_P)", "(?i)spleen",
          #巨噬细胞里tgf高表达与低表达之间比较
          "(?i)(?<!spleen.{0,200}?)tgf_(p|hi)(?!.{0,200}?spleen)",
          "(?i)tgf_(n|lo)", "MAC_TGFβ1_P-N", "(?i)(F4/80_P|cd11b_P)", "(?i)spleen",
          #巨噬细胞里tnf高表达与arg1高表达之间比较
          "(?i)(?<!spleen.{0,200}?)tnf_(p|hi)(?!.{0,200}?spleen)",
          "(?i)arg1?_(p|hi)", "MAC_TNF_P-ARG1_P", "(?i)(F4/80_P|cd11b_P)", "(?i)spleen",
          #nk中cd107a高表达与低表达之间比较
          "(?i)cd107a_(p|hi)", "(?i)cd107a_(n|lo)", "NK-cd107a_P-N", "(?i)NK", NULL,
          #nk中tnf高表达与低表达之间比较
          "(?i)tnf_(p|hi)", "(?i)tnf_(n|lo)", "NK_TNF_P-N", "(?i)NK", NULL,
          #nk中pd1阳性与阴性之间比较
          "(?i)pd.?1_(p|hi)", "(?i)pd.?1_(n|lo)", "NK_pd-1P-N", "(?i)NK", NULL,
          #nk是否在肿瘤中之间比较
          "(?i)(?<=tumor.{0,200}?)nk|nk(?=.{0,200}?(tumor))",
          "(?i)(?<!tumor.{0,200}?)nk(?!.{0,200}?tumor)",
          "NK_tumor-periphery", "(?i)NK", "(?i)pd.?1|cd107a|tnf_"
        )
      search_table[4] <- as.character(search_table[[4]])
      search_table[5] <- as.character(search_table[[5]])
      search_table <- as.data.frame(search_table)
    }

    for(mouse in mouse_index_list){
      single_mouse_data <-
        data_for_analysis %>%
        filter(mouse_index == mouse, mouse_strain == strain, !is.na(median_reads)) %>%
        select(sample, mouse_index, unique_num) %>%
        arrange(sample)
      #各老鼠中循环查找表中信息
      for(i in 1:dim(search_table)[1]){

        unique_num_output <-
          rbind(unique_num_output,
                find_sample_by_regex(single_mouse_data = single_mouse_data,
                                     exp_regex = search_table[i, 1],
                                     con_regex = search_table[i, 2],
                                     group = search_table[i, 3],
                                     limit = unlist(search_table[i, 4]),
                                     except = unlist(search_table[i, 5])
                )
          )
      }
    }

    unique_num_output <-
      unique_num_output %>%
      dplyr::filter(!is.na(exp_index), !is.na(con_index)) %>%
      arrange(group)

    return(unique_num_output)
  }
#计算用于决定上下游的阈值
calculation_of_threshold <- function(sub_data_f = sub_data_f, top = top, up_lfc = up_lfc, down_lfc = down_lfc){
  sub_data_f <- sub_data_f %>% group_by(Name) %>% mutate(up_num = sum(lfc >= 0), down_num = sum(lfc < 0))

  if(is.null(top)){
    up_grna <- sub_data_f %>% dplyr::filter(lfc > up_lfc)
    down_grna <- sub_data_f %>% dplyr::filter(lfc < down_lfc)

    up_grna <- up_grna %>% mutate(for_out = paste0(sample, ": ", round(lfc,2)))
    down_grna <- down_grna %>% mutate(for_out = paste0(sample, ": ", round(lfc,2)))
  }else{
    up_grna <-
      sub_data_f %>%
      group_by(sample) %>%
      dplyr::filter(lfc > quantile(lfc, 1-top, na.rm = T))
    down_grna <-
      sub_data_f %>%
      group_by(sample) %>%
      dplyr::filter(lfc < quantile(lfc, top, na.rm = T))

    up_grna <- up_grna %>% mutate(for_out = paste0(sample, ": ", round(lfc,2)))
    down_grna <- down_grna %>% mutate(for_out = paste0(sample, ": ", round(lfc,2)))
  }

  return(list(up_grna, down_grna))
}
#初始化计算各种数值
get_lfc <- function(ana_data = ana_data, con_data = con_data, ana_index = ana_index, con_index = con_index,
                    batch = batch, pair = pair, ana = ana, con = con, fc = fc){
  if(batch){
    ana_data <- ana
    con_data <- con
  }else{
    for(n in 1:length(ana_index)){
      ana_data <- rbind(ana_data,data_for_analysis[data_for_analysis$unique_num == ana_index[n],])
      con_data <- rbind(con_data,data_for_analysis[data_for_analysis$unique_num == con_index[n],])
    }
  }

  #根据是否配对来获得不同的计算数值
  if(pair){
    if(length(ana_index) != length(con_index)){
      stop("长度不等！")
    }
    ana_data$ana_vs_con <- paste0(ana_data$sample," / ",con_data$sample)
    ana_data$num_for_view <- paste0(ana_data$unique_num,"/",con_data$unique_num)
    ana_data$ratio_FC <- (ana_data$besides_ratio + 1e-20) / (con_data$besides_ratio + 1e-20)
    ana_data$median_FC <- (ana_data$median_reads + 1e-20) / (con_data$median_reads + 1e-20)
    ana_data$con_median_read <- con_data$median_reads
  }else{
    new_a <- ana_data %>% group_by(Name) %>% summarise(ratio = mean(besides_ratio), median = mean(median_reads))
    new_c <- con_data %>% group_by(Name) %>% summarise(ratio = mean(besides_ratio), median = mean(median_reads))

    ana_data$ana_vs_con <- paste0(paste0(unique(ana_data$sample), collapse = ","), " / ", paste0(unique(con_data$sample), collapse = ","))
    ana_data$ratio_FC <- (new_a$ratio + 1e-20) / (new_c$ratio + 1e-20)
    ana_data$median_FC <- (new_a$median + 1e-20) / (new_c$median + 1e-20)
    ana_data$con_median_read <- new_c$median
  }

  if(fc == "median_FC"){
    ana_data$lfc <- log2(ana_data$median_FC)
  }
  if(fc == "ratio_FC"){
    ana_data$lfc <- log2(ana_data$ratio_FC)
  }

  return(ana_data)
}
#生成输出矩阵
generate_output_matrix <- function(ana_data = ana_data, sub_data_f = sub_data_f, up_grna = up_grna,
                                   down_grna = down_grna, pair = pair, passed_list = passed_list){
  table_collected <-
    tibble(Name = NA_character_, mean_reads = NA_real_, mean_ratio = NA_real_,
           median_medianLFC = NA_real_, median_ratioLFC = NA_real_, sd_ratioLFC = 0,
           order = 0, sample_num = "", sample_reads = NA_character_,
           control_reads = NA_character_, reads_without_cas = NA_character_,
           trend = NA_character_, sample = NA_character_, worked_mice = NA_character_,
           unique_num = NA_character_)

  i <- 1
  for (f in list(up_grna, down_grna)) {
    grna_list <- vector()
    for (name in f$Name) {
      if(!name %in% grna_list & name %in% passed_list){
        #对列表开始填充
        temp_data_f <- f %>% dplyr::filter(Name == name) %>% ungroup()
        other_temp <- ana_data %>% dplyr::filter(Name == name) %>% ungroup()

        table_collected[i,"Name"] <- name

        table_collected[i,"mean_reads"] <- temp_data_f %>% summarise(mean_read = mean(median_reads))

        table_collected[i,"mean_ratio"] <- temp_data_f %>%
          summarise(temp = mean(besides_ratio)) %>%
          round(3)

        table_collected[i,"median_medianLFC"] <- temp_data_f %>%
          summarise(temp = log2(median(median_FC))) %>%
          round(2)

        table_collected[i,"median_ratioLFC"] <- temp_data_f %>%
          summarise(mean_fc = log2(median(ratio_FC))) %>%
          round(2)

        table_collected[i,"sd_ratioLFC"] <- temp_data_f %>%
          summarise(temp = sd(ratio_FC)) %>%
          round(2)

        table_collected[i,"order"] <- temp_data_f %>% select(order) %>% unique()

        table_collected[i,"sample_num"] <- temp_data_f %>%
          summarise(temp = paste0(n(),"/",length(unique(sub_data_f$unique_num))))

        table_collected[i,"sample_reads"] <- other_temp %>%
          summarise(sample = paste0(round(median_reads,0), collapse = "; "))

        table_collected[i,"control_reads"] <- other_temp %>%
          summarise(sample = paste0(round(con_median_read,0), collapse = "; "))

        table_collected[i,"reads_without_cas"] <- other_temp %>%
          summarise(sample = paste0(round(control_reads,0), collapse = "; "))

        table_collected[i,"trend"] <-
          temp_data_f %>%
          dplyr::filter(order(unique_num) == 1) %>%
          summarise(paste0("more: ", up_num, "; lesser: ", down_num))

        table_collected[i,"sample"] <- temp_data_f %>% summarise(sample = paste0(ana_vs_con, collapse = "; "))

        table_collected[i,"worked_mice"] <-
          paste0(unlist(str_extract_all(temp_data_f$ana_vs_con, "^\\w?\\d+")), collapse = "; ")

        table_collected[i,"unique_num"] <- temp_data_f %>% summarise(sample = paste0(num_for_view, collapse = "; "))

        i = i + 1
        grna_list <- c(grna_list, name)
      }
    }
  }
  table_collected <- table_collected %>%  arrange(desc(sample_num), desc(abs(median_medianLFC)))
  return(table_collected)
}

#单次来看gRNA在两组样品间的变化
variational_grna_specified_sample <-
  function(data_for_analysis, ana_index = c(), con_index = c(), min_reads = 0, min_sample = 0, batch = T,
           ana = NULL, con = NULL, pair = T, top = NULL, up_lfc = 0.58, down_lfc = -1, keys = NULL, line = 1610,
           fc = "median_FC", remove_nc = F){
    #初始化并计算一下
    ana_data <- data_for_analysis[0,]
    con_data <- data_for_analysis[0,]
    sub_data_f <- data_for_analysis[0,]
    if(remove_nc){
      data_for_analysis <- data_for_analysis %>% dplyr::filter(type != "control")
    }

    #计算惰性的参数
    ana_index_temp = ana_index
    con_index_temp = con_index
    min_reads_temp = min_reads
    min_sample_temp = min_sample
    batch_temp = batch
    ana_temp = ana
    con_temp = con
    pair_temp = pair
    top_temp = top
    up_lfc_temp = up_lfc
    down_lfc_temp = down_lfc
    keys_temp = keys
    line_temp = line
    fc_temp = fc

    #获得各种lfc的数值，用于后续计算
    ana_data <- get_lfc(ana_data = ana_data, con_data = con_data, ana_index = ana_index,
                        con_index = con_index, batch = batch, pair = pair, ana = ana, con = con, fc = fc)

    #如果有key就把带有key的提取出来
    if(!is.null(keys)){
      for(key in keys){
        temp_data_f <- ana_data %>% dplyr::filter(mouse_strain == line,
                                           str_detect(ana_vs_con,regex(key, ignore_case = T)))
        sub_data_f <- rbind(sub_data_f,temp_data_f)
      }
    }else{
      sub_data_f <- ana_data
    }

    #对reads较低的进行一次过滤
    filter_determination <- data_for_analysis %>% dplyr::filter(unique_num %in% c(ana$unique_num, con$unique_num))
    temp_sum <- filter_determination %>% group_by(Name) %>% summarise(n = sum(NumReads > min_reads))
    passed_list <- temp_sum$Name[temp_sum$n >=  min_sample] %>% as.character()

    #决定阈值
    threshold <- calculation_of_threshold(sub_data_f = sub_data_f, top = top, up_lfc = up_lfc, down_lfc = down_lfc)
    up_grna <- threshold[[1]]
    down_grna <- threshold[[2]]

    #将每个样本中同样的gRNA合并起来
    table_collected <- generate_output_matrix(ana_data = ana_data, sub_data_f = sub_data_f, up_grna = up_grna,
                                              down_grna = down_grna, pair = pair, passed_list = passed_list)

    #输出交集与独有的数据
    if(!batch){
      venn_info <- table(table_collected$sample)
      names(venn_info) <- str_replace_all(names(venn_info), "/", "_")
      names(venn_info) <- str_extract_all(names(venn_info), "\\b\\w?\\d+")
      print(venn_info)
    }

    return(table_collected)
  }



#' Title
#'
#' @param data_for_analysis
#' @param samples_csv_file
#' @param unique_num_info
#' @param remove_grna
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
variational_grna_specified_sample_batch <-
  function(data_for_analysis, samples_csv_file = NULL, unique_num_info = NULL, remove_grna = c(), ...){
    if(!is.null(samples_csv_file)){
      samples_info <- read.csv(samples_csv_file, header = T)
    }else{
      samples_info <- NULL
    }

    data_for_analysis[0,]
    all_results <- list()
    all_for_compare_data <-
      extract_samples(data_for_analysis, samples_info = samples_info,
                      unique_num_info = unique_num_info)
    data_for_analysis <- data_for_analysis %>% dplyr::filter(!Name %in% remove_grna)
    i = 1

    for(group_num in sort(unique(all_for_compare_data$group))){
      temp_data <- all_for_compare_data[all_for_compare_data$group == group_num,]
      ana = temp_data[temp_data$numerator == "yes",]
      con = temp_data[temp_data$numerator == "no",]

      if(unique(ana$mouse_strain != con$mouse_strain)){
        browser()
        break
      }

      temp_results <- variational_grna_specified_sample(data_for_analysis,
                                                        ana = ana,
                                                        con = con, ...)
      if(!is.na(temp_results$Name[1])){
        all_results[[paste0(group_num, "_", ana$mouse_strain[[1]])]] <- temp_results
      }
      i = i + 1
    }
    return(all_results)
  }


