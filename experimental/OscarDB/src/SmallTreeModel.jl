struct SmallTreeModel
  _id::String
  model::GroupBasedPhylogeneticModel{PhylogeneticTree{QQFieldElem}}
  tree::String
  model_type::String
end


group_based_phylogenetic_model(stm::SmallTreeModel) = stm.model
phylogenetic_model(stm::SmallTreeModel) = phylogenetic_model(group_based_phylogenetic_model(stm))
