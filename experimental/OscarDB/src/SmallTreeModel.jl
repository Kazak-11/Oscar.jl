struct SmallTreeModel
  _id::String # model encoding id, example 3-0-0-JC
  model::GroupBasedPhylogeneticModel
  model_type::String
  n_leaves::Int

  function SmallTreeModel(model_encoding::String, pm::GroupBasedPhylogeneticModel)
    n_leaves, _, _, model_type = split(model_encoding, "-")
    new(model_encoding, pm, String(model_type), parse(Int, n_leaves))
  end
end

group_based_phylogenetic_model(stm::SmallTreeModel) = stm.model
phylogenetic_model(stm::SmallTreeModel) = phylogenetic_model(group_based_phylogenetic_model(stm))

