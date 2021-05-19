ligrec <- compile_ligrec_descr()

# ligrecx <- list("CellChatDB" = ligrec$CellChatDB)
#
# ligrec <- ligrec %>%
#     ligrec_decomplexify
#
# ligrec_decomplex$CellChatDB$interactions %>%
#     distinct_at(.vars = c("target", "source"))


ligrec_olap <- ligrec %>%
    ligrec_decomplexify %>%
    ligrec_overlap %>%
    summarize_overlaps %T>%
    total_unique_bar

ligrec_olap$interactions

# ligrec_olap$interactions %>%
    # filter(name == 'total')









ligrec_uniq <- ligrec_olap %>%
    map2(
        names(.),
        function(data, label){

            data %<>%
                filter(name != 'total') %>%
                mutate(
                    resource = factor(
                        resource,
                        levels = res_order,
                        ordered = TRUE
                    ),
                    name = factor(
                        name,
                        levels = c('n_shared', 'n_unique'),
                        ordered = TRUE
                    )
                )

            return(data)

        }
    )


data <- ligrec_uniq %>%
    bind_rows(.id="entity") %>%
    mutate(entity = str_to_title(entity)) %>%
    group_by(entity, resource) %>%
    mutate(perc = round((value/sum(value))*100))


ggplot(data, aes(y = resource, x = value, fill = name)) +
    geom_col() +
    scale_fill_manual(
        values = c('#B3C5E9', '#4268B3'),
        label = c(n_shared = 'Shared', n_unique = 'Unique'),
        guide = guide_legend(title = '')
        ) +
    xlab('') +
    ylab('Resources') +
    facet_grid(.~entity, scales='free_x')  +
    geom_text(
        aes(label = if_else(perc > 5, str_glue("{round(perc)}%"), "")),
        position = position_stack(vjust = 0.6), size = 6, color = "white"
    ) +
    theme_bw() +
    theme(
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 17),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size=18),
        legend.text = element_text(size=16),
        strip.text.x = element_text(size=17),
        legend.key.size = unit(17, 'mm')
    )










total_unique_bar <- function(ligrec_olap){

    log_success('Drawing overlap barplots.')

    res_order <-
        ligrec_olap$interactions %>%
        filter(name == 'total') %>%
        arrange(value) %>%
        pull(resource) %>%
        unique

    ligrec_olap %<>%
        map2(
            names(.),
            function(data, label){

                data %<>%
                    filter(name != 'total') %>%
                    mutate(
                        resource = factor(
                            resource,
                            levels = res_order,
                            ordered = TRUE
                        ),
                        name = factor(
                            name,
                            levels = c('n_shared', 'n_unique'),
                            ordered = TRUE
                        )
                    )
            }
        )  %>%
        bind_rows(.id="entity") %>%
        mutate(entity = str_to_title(entity)) %>%
        group_by(entity, resource) %>%
        mutate(perc = round((value/sum(value))*100))


    ligrec_olap %<>% rename(typ)


     ggplot(ligrec_olap, aes(y = resource, x = value, fill = typ)) +
         geom_col() +
            scale_fill_manual(
                values = c('#B3C5E9', '#4268B3'),
                label = c(n_shared = 'Shared', n_unique = 'Unique'),
                guide = guide_legend(title = '')
            ) +
            xlab('') +
            ylab('Resources') +
            facet_grid(.~entity, scales='free_x')  +
            geom_text(
                aes(label = if_else((perc > 3 & typ!="n_shared"), str_glue("{round(perc)}%"), "")),
                position = position_stack(vjust = 0.6, reverse = TRUE),
                size = 8, color = "black",
                inherit.aes = TRUE, check_overlap = TRUE
            ) +
            theme_bw() +
            theme(
                axis.text.x = element_text(size = 16),
                axis.text.y = element_text(size = 17),
                axis.title.y = element_text(size = 20),
                panel.grid.major = element_blank(),
                panel.background = element_blank(),
                axis.ticks = element_blank(),
                legend.title = element_text(size=18),
                legend.text = element_text(size=16),
                strip.text.x = element_text(size=17),
                legend.key.size = unit(17, 'mm')
            )

    cairo_pdf(
        figure_path('size_overlap_combined.pdf'),
        width = 19,
        height = 10,
        family = 'DINPro'
    )

    print(p)

    dev.off()

}

ligrec_olap <- ligrec %>%
    ligrec_decomplexify %>%
    ligrec_overlap %>%
    summarize_overlaps %T>%
    total_unique_bar
