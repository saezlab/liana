require(OmnipathR)
require(dplyr)
require(stringr)
require(assertthat)

omni_list <- list(
  'CellChatDB',
  'CellPhoneDB',
  'Ramilowski2015',
  'Baccin2019',
  'LRdb',
  'Kirouac2010',
  'ICELLNET',
  'iTALK',
  'EMBRACE',
  'HPMR',
  'Guide2Pharma',
  'connectomeDB2020',
  'talklr',
  'CellTalkDB',
  'OmniPath'
)

# Collect all ligrec interactions from all resources
all_ligrec_interactions <- import_intercell_network()

# Remove complexes, check for duplicates, and give all interactions a unique label to enable UpSet plot:
OmniPath_ligrec_interactions <- all_ligrec_interactions %>% 
  dplyr::filter(entity_type_intercell_source != 'complex' & entity_type_intercell_target != 'complex') %>%
  distinct() %>%
  tibble::rowid_to_column()

#' create_lr_networks
#'
#' @param omni_list The list of resources for which to retrieve ligand-receptor
#'   interactions, including 'OmniPath' if desired.
#'
#' @return A list of dataframes with the ligand-receptor interactions contained
#'   in each resource
#' @export
#'
#' @examples
create_lr_networks <- function(resources){
  # Check resource names are valid:
  stopifnot("One or more invalid resource names"=all(resources %in% c(intersect(get_interaction_resources(), get_intercell_resources()), 'OmniPath')))
  
  get_network <- function(resource){
    cat("Collecting interactions for: ", resource, "\n")
    if(resource == 'OmniPath'){
      return(OmniPath_ligrec_interactions)
    }
    # Try statement because of bug in import_intercell_network
    # 'MatrixDB', 'Adhesome', "CellChatDB" throw error (0 intercellular communication role records).
    try({
      intercell_network <- import_intercell_network(
        interactions_param = list(resources = resource)
      )
      # Remove complexes and any duplicates, remove cols database_intercell_source, database_intercell_target
      intercell_network <- intercell_network %>% 
        dplyr::filter(entity_type_intercell_source != 'complex' & entity_type_intercell_target != 'complex') %>%
        distinct() %>%
        dplyr::select(-c(database_intercell_source, database_intercell_target))
      # Retrieve unique id amongst all OmniPath interactions to enable UpSet plot:
      intercell_network_id <-  suppressMessages(dplyr::inner_join(OmniPath_ligrec_interactions, intercell_network))#, by=merge_columns)
      # Assert that no interactions are removed in this process:
      assertthat::validate_that(nrow(intercell_network) == nrow(intercell_network_id), 
                                msg="Some interactions may have been removed when assigning ID")
      return(intercell_network_id)
    })
  }
  intercell_networks <- lapply(resources, get_network)
  names(intercell_networks) <- resources
  return(intercell_networks)
}

networks <- create_lr_networks(omni_list)
saveRDS(networks, "./R/networks.RData")
