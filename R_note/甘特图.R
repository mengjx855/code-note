####info####
# created date: 2022-8-10

# 甘特图(Gantt chart) 又称为横道图、条状图(Bar chart)。
# 其通过条状图来显示项目、进度和其他时间相关的系统进展的内在关系随着时间进展的情况。
# 一般情况下，横轴表示时间，纵轴表示各个项目，线条表示期间计划和实际完成的情况.


####ggplot2####
library(tidyverse)
library(ggtext)
library(hrbrthemes)
library(wesanderson)
library(LaCroixColoR)
library(RColorBrewer)
library(ggsci)

# 构建数据
data <- data.frame(name = c('A', 'B', 'C', 'D','E','F'), 
                   start = c(4, 7, 12, 16,10,8),
                   end = c(12, 11, 8, 10,13,14),
                   shift_type = c('Early', 'Mid_day', 'Mid_day', 'Late','Mid_day','Early'))
# 数据可视化
ggplot(data, aes(x = start, xend = end, y = name, yend = name, color = shift_type)) +
  geom_segment(size = 8, lineend = "round") +
  ggsci::scale_color_nejm() +
  labs(title = "Example of ggplot2::Gantt Chart ",
    subtitle = "processed charts with geom_segment()",
    caption = "Visualization by DataCharm") +
  #hrbrthemes::theme_ft_rc(base_family = "Roboto Condensed")  +
  hrbrthemes::theme_ipsum(base_family = "Roboto Condensed")  +
  theme(plot.title = element_markdown(hjust = 0.5, vjust = .5, color = "black",size = 20, margin = margin(t = 1, b = 12)),
        plot.subtitle = element_markdown(hjust = 0, vjust = .5, size = 15),
        plot.caption = element_markdown(face = 'bold', size = 12))


####ganttrify####
# R-ganttrify包专门绘制绘制甘特图，其默认的排版和格式都是符合一般的审美的，其还有很多用于定制化的参数.
library("ganttrify")
ganttrify(project = ganttrify::test_project,
          project_start_date = "2021-03",
          font_family = "Roboto Condensed") +
  labs(title = "Example of ganttrify::Gantt Chart ",
       subtitle = "processed charts with geom_segment()",
       caption = "Visualization by DataCharm") +
  theme(plot.title = element_markdown(hjust = 0.5, vjust = .5, color = "black",size = 20, margin = margin(t = 1, b = 12)),
        plot.subtitle = element_markdown(hjust = 0, vjust = .5, size = 15),
        plot.caption = element_markdown(face = 'bold', size = 12))

ganttrify(project = ganttrify::test_project,
          project_start_date = "2021-03",
          font_family = "Roboto Condensed") +
  labs(title = "Example of ganttrify::Gantt Chart ",
       subtitle = "processed charts with geom_segment()",
       caption = "Visualization by DataCharm") +
  theme(plot.title = element_markdown(hjust = 0.5, vjust = .5, color = "black",size = 20, margin = margin(t = 1, b = 12)),
        plot.subtitle = element_markdown(hjust = 0, vjust = .5, size = 15),
        plot.caption = element_markdown(face = 'bold', size = 12))
