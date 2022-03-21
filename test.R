#加载所有函数 用于测试
devtools::load_all()

load("D:/research/my_lab/MARC_bioinformatics/amplicon/result/RData/HP7-65 all data.RData")
rm(list = ls()[!ls() %in% c("data_for_analysis")])


mouse_list <-
  data_for_analysis %>%
  filter(mouse_strain == 1913, order == 1, cas == "yes", str_detect(sample, "(?i)nk"), batch == "HP-65") %>%
  select(mouse_index) %>%
  unique() %>%
  unlist()

unique_info <- get_unique_num_table(data_for_analysis, mouse_list, 1913)

tcell <-
  variational_grna_specified_sample_batch(data_for_analysis, remove_grna = "Cd47",
                                          unique_num_info = unique_info,
                                          up_lfc = 0.7, down_lfc = -0.7)

hp_65 <-
  data_clean_for_ngs("D:/research/my_lab/MARC_bioinformatics/amplicon/result/2022-02-06 HP-65.tsv",
                     sample_path = "D:/research/my_lab/MARC_bioinformatics/amplicon/result/sample_info/sample_info HP-65.xlsx")
