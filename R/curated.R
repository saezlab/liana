
LR_RESOURCES <- c(
    'CellPhoneDB',
    'Cellinker',
    'CellTalkDB',
    'CellChatDB',
    'CellCall',
    'connectomeDB2020',
    'Guide2Pharma',
    'Baccin2019',
    'Kirouac2010',
    'Ramilowski2015',
    'scConnect',
    'talklr',
    'ICELLNET',
    'EMBRACE',
    'LRdb',
    'iTALK',
    'SignaLink',
    'HPMR'
)

CURATED_LR_RESOURCES <- c(
    'Guide2Pharma',
    'HPMR',
    'ICELLNET',
    'Kirouac2010',
    'CellTalkDB',
    'CellChatDB',
    'connectomeDB2020',
    'SignaLink',
    'Cellinker',
    'CellPhoneDB',
    'talklr'
)


#' Statistics about literature curated L-R interactions
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom rlang set_names
#' @importFrom tidyr unnest_wider
#' @importFrom tibble tibble
#' @export
curated_stats <- function(){

    # NSE vs. R CMD check workaround
    resource <- NULL

    LR_RESOURCES %>%
    set_names(
        map(., stats_one_resource),
        .
    ) %>%
    tibble(
        resource = names(.),
        data = .
    ) %>%
    unnest_wider(col = data)

}


#' Statistics about curated interactions in one resource
#'
#' @importFrom magrittr %>% inset2
#' @importFrom purrr map
#' @importFrom rlang set_names
#' @noRd
stats_one_resource <- function(resource){

    c('total', 'literature', 'curated') %>%
    set_names(
        map(
            .,
            function(field){
                get(sprintf('%s_one_resource', field))(resource)
            }
        ),
        .
    ) %>%
    map(nrow) %>%
    inset2(
        'literature_pct',
        .$literature / .$total * 100L
    ) %>%
    inset2(
        'curated_pct',
        .$curated / .$total * 100L
    )

}


#' Curated interactions from one resource
#'
#' @param resource Character: name of a single resource
#'
#' @noRd
curated_one_resource <- function(resource){

    # some resources have their dedicated function:
    suffixes <- c('curated', 'ligrec')

    for(suf in suffixes){

        func_name <- sprintf('%s_%s', tolower(resource), suf)

        func <- tryCatch(
            get(func_name, envir = asNamespace('OmnipathR')),
            error = function(cond){NULL}
        )

        if(!is.null(func)){

            return(func())

        }

    }

    # at the end we fall back to the generic method:
    # interactions with literature references from the resource
    literature_one_resource(resource = resource)

}


#' All interactions from one resource
#'
#' @importFrom magrittr %<>%
#' @importFrom OmnipathR import_post_translational_interactions
#' @noRd
total_one_resource <- function(resource){

    resource %<>% network_resource_name

    import_post_translational_interactions(resources = resource)

}


#' Interactions with literature references from one resource
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr filter
#' @importFrom stringr str_detect
#' @importFrom OmnipathR import_post_translational_interactions
#' @importFrom OmnipathR with_references
#' @noRd
literature_one_resource <- function(resource){

    # NSE vs. R CMD check workaround
    references <- NULL

    resource %<>% network_resource_name

    import_post_translational_interactions(resources = resource) %>%
    with_references(resources = resource)

}


#' Network resource name
#'
#' @importFrom magrittr %>%
#' @noRd
network_resource_name <- function(resource){

    resource %>%
    {`if`(. == 'SignaLink', 'SignaLink3', .)}

}
