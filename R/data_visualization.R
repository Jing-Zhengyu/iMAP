#' Title
#'
#' @param data_for_analysis
#' @param plot_index
#' @param grna_index
#' @param smooth
#'
#' @return
#' @export
#'
#' @examples
marc_distribution_inturn <-
  function(data_for_analysis, plot_index, grna_index = 2:100, smooth = F){

    s <- geom_smooth(aes(group = 1),size = 0.8,color = "#BE9063",span = 0.22,se = F)

    #循环出图
    for (plot_num in plot_index) {
      if(plot_num %in% data_for_analysis$unique_num){
        data_for_plot <- filter(data_for_analysis, unique_num == plot_num, order %in% grna_index)
        data_for_plot$Name <- forcats::fct_inorder(as.character(data_for_plot$Name))
        p <- ggplot(data_for_plot,
                    aes(Name,percent,fill=type)) +
          geom_bar(stat = "identity") +
          theme(axis.text.x  = element_text(angle=45, vjust=0.5,face = "bold",size = 8)) +
          xlab("gRNA") +
          ggtitle(paste0("num.",plot_num,"-",unique(data_for_plot$sample)))

        #判断是否有拟合曲线
        if(smooth){
          p1 <- p + s
        }else{
          p1 <- p
        }
        print(p1)
      }
    }
  }

#' Title
#'
#' @param data_for_analysis
#' @param plot_index
#' @param y
#' @param group
#' @param logy
#' @param order
#' @param order_index
#' @param grna_index
#'
#' @return
#' @export
#'
#' @examples
marc_distribution_comparative <-
  function(data_for_analysis, plot_index, y = "percent", group = "sample",
           logy = F, order = F, order_index = NULL, grna_index = 2:100){
    if(!order){
      data_for_analysis$Name <- forcats::fct_inorder(as.character(data_for_analysis$Name))
    }else{
      for_order <- data_for_analysis %>% filter(unique_num == order_index)
      name <- for_order$Name[order(for_order[y])]
      data_for_analysis$Name <- factor(data_for_analysis$Name, levels = name)
    }
    plot_data <- filter(data_for_analysis, marker %in% plot_index, order %in% grna_index)

    if(logy){
      p <- ggplot(plot_data) +
        geom_line(aes(Name, log2(get(y)),color = get(group),group = get(group)),size = 1) +
        geom_point(data = filter(plot_data, type == "control"),
                   aes(Name,log2(get(y))), color = "black", size = 1.5) +
        scale_color_d3() +
        theme_classic() +
        theme(axis.text.x  = element_text(angle=90, vjust=0.5,face = "bold",size = 8)) +
        xlab("gRNA") +
        ylab(y)
    }else{
      p <- ggplot(plot_data) +
        geom_line(aes(Name, get(y),color = get(group), group = get(group)),size = 1) +
        geom_point(data = filter(plot_data, type == "control"),
                   aes(Name,get(y)), color = "black", size = 1.5) +
        scale_color_d3() +
        theme_classic() +
        theme(axis.text.x  = element_text(angle=90, vjust=0.5,face = "bold",size = 8)) +
        xlab("gRNA") +
        ylab(y)
    }

    print(p)

  }

#' Title
#'
#' @param data_for_analysis
#' @param plot_index
#' @param FC
#' @param grna_index
#' @param logy
#'
#' @return
#' @export
#'
#' @examples
marc_logfc_inturn <-
  function(data_for_analysis, plot_index,  FC = "ratio_FC", grna_index = 2:100,
           logy = T){

    data_for_analysis <- data_for_analysis %>% filter(!is.na(get(FC)), order %in% grna_index)
    index_FC <- which(colnames(data_for_analysis) == FC)

    for (plot_num in plot_index) {
      if(plot_num %in% data_for_analysis$unique_num){
        single_data <- filter(data_for_analysis,unique_num == plot_num)
        single_data <- single_data[order(single_data[index_FC]),] %>% as_tibble()
        single_data$Name <- factor(single_data$Name,levels = single_data$Name)

        single_data$fill <- case_when(
          single_data$type == "control" ~ "#00b0f0",
          single_data$type == "unknown" ~ "lightgray"
        )
        if(logy){   #看是否要取对数
          single_data[FC] <- log2(single_data[FC])
        }

        #计算展现多少区域的数据
        if(min(single_data[index_FC]) < -8){
          lim <- range(single_data[index_FC])
          lim[1] <- -8
          lim[2] <- lim[2] + 0.5
        }else{
          lim <- range(single_data[index_FC]) + c(-0.5,0.5)
        }
        ratio <- NA

        #判断不同品系小鼠的标签高度
        if(unique(single_data$mouse_strain == 1329)){
          ratio <- 0.035 * (lim[[2]] - lim[[1]])
        }else{
          ratio <- 0.09 * (lim[[2]] - lim[[1]])
        }

        #计算label的坐标
        single_data$for_label_y <- case_when(
          single_data[index_FC] < 0 ~ ratio,
          single_data[index_FC] >= 0 ~ -1 * ratio)

        #计算纵坐标标签个数和大小
        top_section_percent <- lim[[2]] / (lim[[2]] - lim[[1]])
        num_top_break <- round(top_section_percent * 6)
        top_break <- seq(0,round(lim[[2]]),length.out = num_top_break + 1)
        bottom_break <- seq(0,round(lim[[1]]),length.out = (6 - num_top_break) + 1)
        break_y <- c(rev(bottom_break),top_break[-1]) %>% round(1)

        #获取样品信息当做标题
        plot_title <- paste0(plot_num,"_",unique(single_data$sample))


        #复杂的画图函数
        #判断哪个品系，做些微调
        if(unique(single_data$mouse_strain == 1329)){
          p <-
            ggplot(single_data,aes(Name,get(FC))) +
            geom_bar(aes(fill=gene),stat = "identity",color = "black", alpha = 1) +
            geom_text(aes(x = Name, y = for_label_y,label = serial_number),size = 4) +
            #               geom_point(aes(Name,0,color=gene),size=3) +
            coord_cartesian(ylim=lim) +
            scale_y_continuous(breaks = break_y, labels = break_y) +
            ylab("FC") +
            theme_classic() +
            scale_fill_manual(values = c("#FF1493","#FFA54F","#7e8aa2",
                                         "grey","#436EEE","#00BFFF"),
                              limits = c("Polr2a","Kpnb1","Cd19","Cd45","Tp53","Pten"), #修改图例顺序
                              labels = c("Polr2a","Kpnb1","Cd19","Cd45","Tp53","Pten")) +
            #               scale_y_continuous(limits = c()) +
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.line.x = element_blank(),
                  axis.title.y = element_text(size = 18),
                  axis.text.y = element_text(size = 15),
                  legend.position = c(0.45,0.9),
                  legend.title =  element_blank(),
                  legend.direction = "horizontal",
                  legend.box = "vertical",
                  legend.text = element_text(
                    margin = margin(r = 5, unit = "pt")),
                  plot.title = element_text(size = 25,hjust = 0.5)
            ) +
            guides(fill = guide_legend(nrow = 1)) +
            ggtitle(plot_title)
        }else{
          #计算每个标签的翻转方向
          angle <- rep(c(90,270),c(sum(single_data[index_FC] < 0),
                                   sum(single_data[index_FC] >= 0)))
          p <-
            ggplot(single_data,aes(Name,get(FC))) +
            geom_bar(stat = "identity",color = "white",alpha = 1, fill = single_data$fill) +
            geom_bar(data = filter(single_data, gene == "Polr2a"), aes(Name,get(FC)),
                     fill = "#ff7f0e", stat = "identity",color = "#ff7f0e", alpha = 1) +
            geom_text(aes(x = Name, y = for_label_y,label = Name),
                      size = 3, angle = angle) +
            coord_cartesian(ylim=lim) +
            ylab("FC") +
            scale_y_continuous(breaks = break_y, labels = break_y) +
            theme_classic() +
            theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.line.x = element_blank(),
                  axis.title.y = element_text(size = 18),
                  axis.text.y = element_text(size = 15),
                  legend.position = c(0.45,0.9),
                  legend.title =  element_blank(),
                  legend.direction = "horizontal",
                  legend.box = "vertical",
                  legend.text = element_text(
                    margin = margin(r = 5, unit = "pt")),
                  plot.title = element_text(size = 25,hjust = 0.5)
            ) +
            guides(fill = guide_legend(nrow = 1)) +
            ggtitle(plot_title)
        }
        print(p)
      }
    }
  }

#' Title
#'
#' @param data
#' @param unique_n
#' @param control_unique_n
#' @param y
#' @param show_name
#' @param grna_index
#' @param remove_order
#' @param limit
#' @param seq
#'
#' @return
#' @export
#'
#' @examples
line_with_control <-
  function(data, unique_n, control_unique_n = 0, y = "percent", show_name = T,
           grna_index = c(1:1000), remove_order = c(), limit = 0.5, seq = 0.1){
    single_data <- data %>% filter(unique_num == unique_n)
    single_data$geom <- "line"

    if(control_unique_n == 0){
      first_data <- data %>% filter(unique_num == unique_n, order == 1)
      c_data <- data %>% filter(mouse_strain == unique(single_data$mouse_strain), cas =="no",
                                !str_detect(sample,"heart.|thymus"))
      c_data <- c_data %>% filter(order == 1) %>% group_by(sample, unique_num) %>%
        summarise(d = abs(percent[[1]] - first_data$percent[[1]])) %>%
        arrange(d)
      control_data <- data[data$unique_num == c_data$unique_num[[1]], ]
    }else{
      control_data <- data[data$unique_num == control_unique_n, ]
    }
    control_data$geom <- "area"


    plot_data <- rbind(control_data, single_data)
    plot_data %<>% filter(!(order %in% remove_order), order %in% grna_index)
    #更改名字用于x轴坐标
    if(show_name){
      plot_data$name <- paste0((plot_data$order - 1),"-",plot_data$gene)
    }else{
      plot_data$name <- plot_data$order - 1
    }
    plot_data$name <- factor(plot_data$name,levels = unique(plot_data$name))
    plot_data$sample <- factor(plot_data$sample,
                               levels = c(unique(plot_data$sample[plot_data$geom == "area"]),
                                          unique(plot_data$sample[plot_data$geom == "line"])))


    p <- ggplot() +
      geom_area(data = filter(plot_data, geom == "area"),
                aes(name, get(y), group = 1), fill = "#8fbc94", alpha = 1) +
      geom_line(data = filter(plot_data, geom == "line"),
                aes(name, get(y), group = sample), color = "black",
                size = 0.8, linetype = "dotted", alpha = 0.8) +
      geom_point(data = filter(plot_data, geom == "line"),
                 aes(name, get(y),color = type), size = 5) +
      ylab("percentage") +
      theme_classic() +
      scale_color_manual(values = c("#ff7f0e","#1f77b4", "gray"),
                         limits = c("control","unknown", "first gRNA"), #修改图例顺序
                         labels = c("NC","screening genes", "first gRNA")) +
      scale_x_discrete(expand = c(0.01,0)) +
      scale_y_continuous(expand = c(0,0), limits = c(0,limit), breaks = seq(0,limit,seq),
                         labels = paste0(seq(0, limit * 100, seq * 100))) +
      theme(axis.title.x = element_blank(),
            axis.text.x  = element_text(angle=270, hjust=0.05,
                                        vjust=0.3,  size = 15),
            axis.title.y = element_text(size = 20),
            axis.text.y = element_text(size = 20, family="Arial", colour = "black"),
            legend.position = "none"
      )
    print(p)
    return(p)
  }

