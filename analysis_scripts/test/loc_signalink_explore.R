ligrecx <- compile_ligrec_descr()


ligrec <- ligrecx %>%
    ligrec_decomplexify %T>%
    ligrec_overheats %>%
    ligrec_overlap

ligrec %>%  ligrec_classes_bar_enrich('OP-L', location)



ligrec_signal <- ligrec_classes_bar_enrich(ligrec, 'SignaLink_pathway', pathway, NULL)

ligrec_signal$interactions

ligrec_signal_counts <- ligrec_signal$interactions %>%
    filter(resource=="OmniPath") %>%
    group_by(target, pathway)  %>%
    summarise(n = n()) %>%
    mutate(perc = n / sum(n)) %>%
    ungroup()


#
var <- sym("pathway")
var <- ensym(var)

resource = sym("SignaLink_pathway")

var <- enquo(var)
legend_title <- sprintf(
    '%s (%s)',
    var %>% quo_text %>% str_to_title,
    resource %>% str_replace('_.*$', '')
)


# Ligrec axis
attr <- sym("pathway")

attr_str <- quo_text(attr)
attr_src <- sprintf('%s_source', attr_str) %>% sym
attr_tgt <- sprintf('%s_target', attr_str) %>% sym

annot <-
    import_omnipath_annotations(resource = "SignaLink_pathway", wide = TRUE)


ligrec$interactions %<>%
    annotated_network(annot = annot, !!attr) %>%
    filter(!!attr_src == !!attr_tgt) %>%
    select(-!!attr_tgt) %>%
    rename(!!attr := !!attr_src)



lig_data <- ligrec$interactions %>%
    filter(!is.na(!!var)) %>%
    order_by_group_size(resource) %>%
    mutate(
        !!var := factor(
            !!var,
            levels = sort(unique(!!var)),
            ordered = TRUE
        )
    )

p <- ggplot(lig_data, aes(x = resource, fill = !!var)) +
    geom_bar() +
    stat_count() +
    scale_fill_manual(
        values = .palette1,
        guide = guide_legend(title = legend_title)
    ) +
    xlab('Resources') +
    ylab(str_to_title("entity")) +
    theme_bw() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.text = element_text(size = 7),
        legend.key.size = unit(3, 'mm')
    )
