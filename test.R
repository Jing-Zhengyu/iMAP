#加载所有函数 用于测试
devtools::load_all()

load("E:/my_lab/MARC_bioinformatics/amplicon/result/RData/HP7-64 all data.RData")
rm(list = ls()[!ls() %in% c("data_for_analysis", "all_data")])




immune_data <-
  data_for_analysis %>%
  dplyr::filter(mouse_index %in% c(8,264,391,385,212) &
                  !unique_num %in% c() &   #517,1162,1029,1160,700,721,697,1156,851,852,954,959
                  !str_detect(organ, "thymus|epididymis_|NA|brown_fat|gland|bone|muscle|ovary"),
                str_detect(organ, "(?i)spleen.|ln.|blood")
  ) %>%
  group_by(sample, Name) %>%
  dplyr::filter(order(unique_num, decreasing = T) == 1)
compare_table <- tribble(
  ~first_sample, ~second_sample, ~label, ~limit, ~except,
  "B_cell", "CD\\d|treg", "B-T", NULL, NULL,
  "8_e", "8_n", "cd8E-cd8N", NULL, NULL,
  "8_cm", "8_n", "cd8CM-cd8N", NULL, NULL,
  "4_e", "4_n", "cd4E-cd4n", NULL, NULL,
  "8_e", "4_e", "cd8e-cd4e", NULL, NULL,
  "8_n", "4_n", "cd8n-cd4n", NULL, NULL,
  "treg", "4_n", "treg-cd4n", NULL, NULL,
  "treg", "4_e", "treg-cd4e", NULL, NULL,
  "blood", "B_cell", "blood-B", NULL, NULL,
  "blood", "CD\\d|treg", "blood-T", NULL, NULL
)
compare_table <- as.data.frame(compare_table)
info_1610 <- get_unique_num_table(as.data.frame(immune_data),
                                  mouse_index_list = c(8,264,391,385,212),
                                  strain = 1610, search_table = compare_table)
tcell_1610 <-
  variational_grna_specified_sample_batch(immune_data, remove_grna = "Cd47",
                                          unique_num_info = info_1610,
                                          up_lfc = 0.5, down_lfc = -0.5)
