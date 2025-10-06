struct SmallTreeModel
  _id::String # model encoding id, example 3-0-0-JC
  model::GroupBasedPhylogeneticModel
  model_type::String
end


group_based_phylogenetic_model(stm::SmallTreeModel) = stm.model
phylogenetic_model(stm::SmallTreeModel) = phylogenetic_model(group_based_phylogenetic_model(stm))
