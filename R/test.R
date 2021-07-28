unique_info_1913 <-
  get_unique_num_table(data_for_analysis,
                       c("K221", "K226", "K222", "K219", "K134", "K291", "K141", "K273", "K278",
                         "K308", "K140", "K129", "G85", "G206", "K52", "K65", "K67"), 1913)
tcell_1913 <-
  variational_grna_specified_sample_batch(data_for_analysis, remove_grna = "Cd47",
                                          unique_num_info = unique_info_1913,
                                          up_lfc = 0.7, down_lfc = -0.7)





single_mouse_data <- data_for_analysis %>%
       dplyr::filter(mouse_index == "K273", !is.na(median_reads)) %>%
       dplyr::select(sample, mouse_index, unique_num) %>%
       arrange(sample)


find_sample_by_regex(single_mouse_data, "(?i)tnf_(p|hi)", "(?i)tnf_(n|lo)",
                     "NK_TNF_P-N", limit = "(?i)NK", except = NULL
)
