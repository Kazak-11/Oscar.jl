function _overlattice_orbits(L::ZZLat; even=true)
  d = ZZ(det(L))
  D = discriminant_group(L)
  idD = hom(D,D,gens(D))
  G,iG = image_in_Oq(L)
  orders = [i for i in divisors(d) if divides(d,i^2)[1]]
  result = ZZLat[]
  for ord in orders 
    #@show ord, D
    b, l, p = is_prime_power_with_data(ord)
    if b && is_elementary(D, p)
      sg = first.(first.(_isotropic_subspaces_representatives_and_stabilizers_elementary(D, iG, valuation(ord,p);do_stab=false)))
    else 
      # slooow
      sg = domain.(first.(_subgroups_orbit_representatives_and_stabilizers(idD, G, ord)))
    end
    for S in sg 
      M = cover(S)
      if !is_integral(M) || (even && !is_even(M))
        continue
      end 
      push!(result,M)        
    end
  end
  return result
end

function root_overlattices(n::Int)
  result = ZZLat[]
  for R in root_lattices(n)
    for S in _overlattice_orbits(R)
      # only add the ones not adding new roots 
      RS = root_sublattice(S)
      if RS == R # this is terribly inefficient
        push!(result,S)
      end 
      if S != RS
        #@show "new"
      end
    end 
  end 
  return result
end

# Dev comment: initial algorithm in paper initializes Ga as directed edge-vertex weighted graph using adjenct matrix Ga_matrix
# It looks that it's better to transform it directly from initial adjenct matrix Ga_matrix to the resulting graph,
# as there is no good function in Oscar to create a weighted graph from adjenct matrix
function _get_edge_labeled_graph(char_vectors_set, gram_matrix)
    B = Matrix(ZZ, char_vectors_set)
    Ga_matrix = B*gram_matrix*transpose(B)
    a = 1 + maximum(Ga_matrix)
    b = a + 1
    p = number_of_columns(B)
    graph_dict = Dict()
    tmp = 0
    for j = [1; p+2]
        for i = [1; j]
            if i < j && j <= p
                tmp = Ga_matrix[i,j]
            elseif i <= j
                if j == p + 1
                    tmp = Ga_matrix[i,i]
                else 
                    tmp = a
                end
            elseif i == p + 1 && j == p + 2
                tmp = b
            end
            push!(graph_dict, (i,j)=>tmp)
        end
    end

    return graph_from_labeled_edges(Undirected, graph_dict, name = :weight) # T1(Ga)
end

function _lift_canonical_ordering(i, j, w, can_order) 
    can_i = find(can_order .== (1+(w-1)*(i-1)))[0] # returnes index of min Si = (i,0) in canonical ordering 
    can_j = find(can_order .== (1+(w-1)*(j-1)))[0]
    return can_i < can_j
end

function _get_canonical_form(A, char_vectors_set, canonical_ordering)
    can_char_vectors_set = sort(char_vectors_set, canonical_ordering)
    (H, U) = hnf_with_transform(can_char_vectors_set)
    return U*A*transpose*U
end


function canonical_form(L::ZZLat)
    gram = gram_matrix(L)
    n = dim(ambient_space(L))
    if n>=2 && n<5
        char_vectors_set = characteristic_vectors(L)
    elseif n>=5 #maybe it will be removed or condition made higher due to performance
        char_vectors_set = _vor_vectors_set()
    end 
    T1_graph = _get_edge_labeled_graph(char_vectors_set, gram) # transform from adjenctcy matrix A to edge-vertex weighted graph Ga, then to edge weighted graph T1(Ga)
    T2_graph =  _edge_label_to_vertex_label(T1_graph, label = :weight) # transform from T1(Ga) to vertex weighted graph T2(T1(Ga))
    can_order = _canonical_perm(T2_graph, label = :weight) # get canonical ordering of T2(T1(Ga))
    lift_canonical_ordering = (i, j) => _lift_canonical_ordering(i, j, w, can_order)
    return _get_canonical_form(A, char_vectors_set, lift_canonical_ordering)
end