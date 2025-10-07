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
model_type(stm::SmallTreeModel) = stm.model_type
n_leaves(stm::SmallTreeModel) = stm.n_leaves
graph(stm::SmallTreeModel) = graph(group_based_phylogenetic_model(stm))

function Base.show(io::IO, stm::SmallTreeModel)
  println(io, "Small tree phylogenetic model")
  println(io, "Model type: $(model_type(stm))")
  println(io, "Number of leaves: $(n_leaves(stm))")
end

