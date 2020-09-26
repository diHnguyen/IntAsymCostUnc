#USES TREE-SA BUT KEEP A LIMIT BETWEEN M_MAX AND M_MIN
#Ready to upload
using JuMP 
using Gurobi
using LightGraphs
using DataFrames
using CSV
using TimerOutputs

#myInstance = string(parse(Int64,ARGS[1]))
myInstance = "5"
myRun = Dates.format(now(), "HH:MM:SS")
epsilon = 1e-3
outfile = "./Output/MulticutMultipartitionScore"*myInstance*"_"*myRun*".txt"
gurobi_env = Gurobi.Env()

to = TimerOutput()
myFile = "./Instances/20NodesCenter_"*myInstance*".jl"
include(myFile)
setparams!(gurobi_env, Heuristics=0.0, Cuts = 0, OutputFlag = 0)

h1 = Model(solver = GurobiSolver(gurobi_env)) # If we want to add # in Gurobi, then we have to turn of Gurobi's own Cuts 

@variable(h1, 1 >= y_h[1:Len]>=0)
#@variable(h, q[1:Len]>=0)

#Setting constraint for start node
leaving = 1 .* (edge[j=1:Len,1] .== origin)
@constraint(h1, sum(y_h[k]*leaving[k] for k=1:Len) == 1)

#Setting constraint for other nodes

for i in all_nodes
    if i != destination
        if i != origin
            incoming = 1 .* (edge[j=1:Len,2] .== i)
            leaving = -1 .* (edge[j=1:Len,1] .== i)
            @constraint(h1, sum(y_h[k]*leaving[k] + y_h[k]*incoming[k] for k=1:Len) == 0)
        end
    end
end
#@objective(h, Min, sum((cL_orig[i]+d[i]*x_now[i])*y_h[i] + q[i] for i=1:Len))
#print(h1)
#solve(h1)



function what_arc(T, y, e)
    arc_type = -1 #type 0 = basic-SP, 1 = basic-nonSP, 2 = non-basic
    if T[e] == 1
        if y[e] == 1
            arc_type = 0
        else
            arc_type = 1
        end
    else
        arc_type = 2
    end
end

function classify_edges_costs(T,c)
    basic = []
    basic_cost = []
    i_basic = []
    nonbasic = []
    nonbasic_cost = []
    i_nonbasic = []
    
    for e = 1:Len
        #arc_type = what_arc(T,y,e)
        if T[e] != 1
            nonbasic = vcat(nonbasic, edge[e,:]')
            push!(i_nonbasic, e)
            push!(nonbasic_cost, c[e])
        else 
            basic = vcat(basic, edge[e,:]')
            push!(i_basic, e)
            push!(basic_cost, c[e])
        end
    end
    #We need to index back to the original edge values for reference purposes
    return basic, basic_cost, i_basic, nonbasic, nonbasic_cost, i_nonbasic
end

function getSortedBasicArcList(pred, basic, i_basic)
#     succ = pred.*0
    sorted_basic = copy(basic)
    sorted_i_basic = zeros(Int64,length(i_basic))
#     println("pred = ", pred)
#     println("basic = ", basic)
#     println("i_basic = ", i_basic)
    
    index_sorting = length(i_basic)
#     println(index_sorting)
    S = [origin]
#     println("S = ",S)
    while isempty(S) == false
#         println("S = ",S)
        tailNode = S[1]
        deleteat!(S, 1)  
        for i = 1:length(pred)
            #origin vs destination
            if i != tailNode & pred[i] == tailNode
                u, v = tailNode, i
#                 println("(",u , ", ", v,")")
                myArc = transpose([u, v])
                temp_index = find(all(basic .== myArc,2))
                if length(temp_index) > 0
                    sorted_i_basic[index_sorting] = i_basic[temp_index[1]]
                    index_sorting = index_sorting - 1
                    push!(S,i)
                end
                
            end
        end        
    end
    return sorted_i_basic
#     println("sorted_i_basic = ", sorted_i_basic)
end

function find_affected_node(T, pred, e, edge) #NEED TO CONSIDER CYCLES
    S = Int[]
    N = Int[]
    push!(S, edge[e,2])
    push!(N, edge[e,2])
    visit = collect(1:last_node)
    l = length(S)
    visit[edge[e,2]] = 0
    while l > 0 
        node = S[1]
        for i = 1:last_node 
            if pred[i] == node && visit[i] != 0
                push!(S, i)
                push!(N, i)
                visit[i] = 0
            end
            
        end
        S[1], S[length(S)] = S[length(S)], S[1]
        pop!(S)
        l = length(S)
    end
    return N
end

function Tree_basic(N, nonbasic, nonbasic_cost, label_temp, temp)
    delta_pos = []
    delta_pos_ind = []
    δ_pos_ind = 0
    val_pos_temp = 0
    delta_neg = []
    delta_neg_ind = []
    δ_neg_ind = 0
    val_neg_temp = 0
    
    if length(N) > 0
        ##println("nonbasic = ", nonbasic)
        for j in N
            head = find(nonbasic[:,2] .== j)
            if length(head) > 0
                ##println("head = ", head)
                for i in head
#                     println("Arc #", i)
#                     println("Arc name = ", nonbasic[i,:])
                    condd = nonbasic[i,1] in N
                    if condd == false
#                         println("π[",nonbasic[i,2],"] + Δ] - π[",nonbasic[i,1],"] - c_ij ")
#                         println("= ", label_temp[nonbasic[i,2]], " - ", label_temp[nonbasic[i,1]], " - ", nonbasic_cost[i])
                        val_pos_temp = label_temp[nonbasic[i,1]] + nonbasic_cost[i] - label_temp[nonbasic[i,2]]
#                         println("= ", val_pos_temp)
                        push!(delta_pos, val_pos_temp)
                        push!(delta_pos_ind, i)
                    end
                end
            end
        end
#         println("delta_pos = ", delta_pos)
        if length(delta_pos) > 0
            δ_pos_min = minimum(delta_pos)
            δ_pos_ind = find(delta_pos .== δ_pos_min)[1]
#             for i = 1:length(delta_pos)
#                 if i == 1
#                     δ_pos_min = delta_pos[1]
#                     δ_pos_ind = delta_pos_ind[1]
#                 end
#                 if delta_pos[i] < δ_pos_min
#                     δ_pos_min = delta_pos[i]
#                     δ_pos_ind = delta_pos_ind[i]
#                 end
#             end
#             temp = temp - δ_pos_min
        else
            δ_pos_min = -1.0
        end
#         print("δ_pos_min ", δ_pos_min)
        for j in N
            tail = find(nonbasic[:,1] .== j)
            if length(tail) > 0
                ##println("tail = ", tail)
                for i in tail
#                     println("Arc #", i)
#                     println("Arc name = ", nonbasic[i,:])
                    condd = nonbasic[i,2] in N
                    if condd == false
                        ##println("FALSE THEN")
                        val_neg_temp = label_temp[nonbasic[i,1]] + nonbasic_cost[i] - label_temp[nonbasic[i,2]]
#                         println("= ", val_neg_temp)
                        push!(delta_neg, val_neg_temp)
                        push!(delta_neg_ind, i)
                    end
                end
            end
        end
        
    end
    #IF DELTA DOESNT RETURN ANY VALUE AKA EMPTY: LET MIN DELTA = -1 THEN REMOVE IT LATER
    ###println("GOOD HERE")
#     println("delta_neg = ", delta_neg)
    if length(delta_neg) > 0
        δ_neg_min = minimum(delta_neg)
        δ_neg_ind = find(delta_neg .== δ_neg_min)[1]
#         for i = 1:length(delta_neg)
#             if i == 1
#                 δ_neg_min = delta_neg[1]
#                 δ_neg_ind = delta_neg_ind[1]
#             end
#             if delta_neg[i] < δ_neg_min
#                 δ_neg_min = delta_neg[i]
#                 δ_neg_ind = delta_neg_ind[i]
#             end
#         end
        temp = temp - δ_neg_min
    else
        δ_neg_min = -1.0
    end
#     print("δ_neg_min ", δ_neg_min)
#     print(temp)
    return δ_pos_min, δ_pos_ind, δ_neg_min, δ_neg_ind, temp
end

function find_basic_Delta(myNode, Delta_Final_Pos, Delta_Final_Neg, sign, basic, i_basic, nonbasic, i_nonbasic)
    myNonTreeArcList = []
    myTreeArcList = []
    Delta_List = []
    delta = -1
#     println("myNode = ", myNode)
    if sign == "POS"
        myNonTreeArcList = find(all(nonbasic[:,2].== myNode,2))
        Delta_List = Delta_Final_Pos
    end
    if sign == "NEG"
        myNonTreeArcList = find(all(nonbasic[:,1].== myNode,2))
        Delta_List = Delta_Final_Neg
    end
    
    #NONBASIC ARCS
    if length(myNonTreeArcList) > 0
#         println("Nonbasic")
#         println("myNonTreeArcList = ", myNonTreeArcList)
        
        e = i_nonbasic[myNonTreeArcList[1]]
#         println("Arc ", e)
        delta = Delta_Final_Neg[e]
#         println("New delta = ", delta)
        for i = 2 : length(myNonTreeArcList)
            e = i_nonbasic[myNonTreeArcList[i]]
#             println("Arc ", e)
            if delta > Delta_Final_Neg[e] && Delta_Final_Neg[e] >= 0
                delta = Delta_Final_Neg[e]
#                 println("New delta = ", delta)
            end
        end
    end
    
    
    #BASIC ARCS
    myTreeArcList = find(all(basic[:,1].== myNode,2))
#     println("basic = ", basic)
#     println("myTreeArcList = ", myTreeArcList)
#     println("l = ",length(myTreeArcList))
    
#     println("myTreeArcList = ", myTreeArcList)
    if length(myTreeArcList) > 0
#         println("Basic")
        if delta == -1
            e = i_basic[myTreeArcList[1]]
            delta = Delta_List[e]
#             println("New delta = ", delta)
        end
        for i = 1 : length(myTreeArcList)
            e = i_basic[myTreeArcList[i]]
#             println("Arc ", e)
            if delta > Delta_List[e] && Delta_List[e] >= 0
                delta = Delta_List[e]
#                 println("New delta = ", delta)
            end
        end
    end
    return delta 
end


function find_SA_Delta(myNode, Delta_Final_Pos, Delta_Final_Neg, sign, basic, i_basic, nonbasic, i_nonbasic)
    myNonTreeArcList = []
    myTreeArcList = []
    Delta_List = []
    delta = maximum(Delta_Final_Neg)
    entering_arc = 0
    found = false
#     println("myNode = ", myNode)
    if sign == "POS"
        myNonTreeArcList = find(all(nonbasic[:,2].== myNode,2))
        Delta_List = Delta_Final_Pos
    end
    if sign == "NEG"
        myNonTreeArcList = find(all(nonbasic[:,1].== myNode,2))
        Delta_List = Delta_Final_Neg
    end
    
    #NONBASIC ARCS
    if length(myNonTreeArcList) > 0
#         println("Nonbasic")
#         println("myNonTreeArcList = ", myNonTreeArcList)
        
#         e = i_nonbasic[myNonTreeArcList[1]]
#         println("Arc ", e)
#         delta = Delta_Final_Neg[e]
#         println("New delta = ", delta)
        for i = 1 : length(myNonTreeArcList)
            e = i_nonbasic[myNonTreeArcList[i]]
#             println("Arc ", e, "Delta ", Delta_Final_Neg[e])
            if delta >= Delta_Final_Neg[e] && Delta_Final_Neg[e] >= 0
                delta = Delta_Final_Neg[e]
                entering_arc = e
#                 println("New delta = ", delta)
            end
        end
    end
    
    
    #BASIC ARCS
    myTreeArcList = find(all(basic[:,1].== myNode,2))
#     println("basic = ", basic)
#     println("myTreeArcList = ", myTreeArcList)
#     println("l = ",length(myTreeArcList))
    
#     println("myTreeArcList = ", myTreeArcList)
    if length(myTreeArcList) > 0
#         println("Basic")
        if delta == -1
#             e = i_basic[myTreeArcList[1]]
            delta = maximum(Delta_List)
#             println("New delta = ", delta)
        end
        for i = 1 : length(myTreeArcList)
            e = i_basic[myTreeArcList[i]]
#             println("Arc ", e, "Delta ", Delta_List[e])
            if delta >= Delta_List[e] && Delta_List[e] >= 0
                delta = Delta_List[e]
                entering_arc = e
#                 println("New delta = ", delta)
            end
        end
    end
#     println("entering arc inside function: ", entering_arc)
    return delta, entering_arc
end

function findNonBasicArc_Loop(label, T, c_g, e, path, y, pred, basic, i_basic)
    label_temp = copy(label)
    Δ = 0
    δ = []
    δ_ind = 0
    T_temp = copy(T)
    pred_temp = copy(pred)
    c_temp = copy(c_g)
    edge_num = copy(e)
    N = Int[]
    cur_edge = Int[]
    path_temp = copy(path)
    
    new_SP = false
    
    
    u = edge[e,1]
    v = edge[e,2]
    t = pred[v]
    e_temp = [t v]
    ###println(f,"(t v) = (",t, " ", v, ")")
    B = find(all(basic .== e_temp,2)) #find a tree edge to be replaced by nonbasic e
    i = B[1]
#     println("Entering arc = ", e )
#     println("Exiting arc = ", i)
    i = i_basic[i]
    val = label_temp[u] - label_temp[v] + c_temp[e]
    push!(δ, val)
    ###println(f,"δ = ", δ)
    δ_min = minimum(δ)
    c_temp[e] = c_temp[e] - δ_min
    Δ = Δ + δ_min
    pred_temp[v] = u
    ###println(f, "Tree before replace: ", T_temp)
    T_temp[e] = 1
    T_temp[i] = 0
    delta = -0.5
    while new_SP == false
        basic, basic_cost, i_basic, nonbasic, nonbasic_cost, i_nonbasic = classify_edges_costs(T_temp,c_temp)
        #loop = loop+1
        ###println("Loop ", loop)
        ##println( "Current T = ", T_temp)
        ##println( "Current y = ", y)
        ##println("Current label = ", label_temp)
        arc_type = what_arc(T_temp, y, e)
#         println("Double check: Arc ", e, " is now type ", arc_type)
        if arc_type == 0 #THIS STATEMENT IS WORKING
            new_SP = true
            delta = Δ
        end
        if arc_type == 1
            N = find_affected_node(T_temp, pred_temp, e, edge)
            ##println("N = ", N)
            δ_min, δ_ind, temp = SA_basic(N, nonbasic, nonbasic_cost, label_temp, c_temp[e])
            ##println("Δ = ", Δ)
            ##println("δ_min = ", δ_min)
            if δ_min < 0
                delta = -0.5 #δ_min
                new_SP = true
            else
                ##println("nonbasic arcs: ", nonbasic)
                ##println("index of arc moving into the basis: ", δ_ind)
                u = nonbasic[δ_ind, 1]
                v = nonbasic[δ_ind, 2]
                t = pred_temp[v]
                e_temp = [t v]
                ##println("Basic arc becoming nonbasic: ", e_temp)
                ##println("Nonbasic arc becoming basic: ", [u v])
                B = find(all(basic .== e_temp,2)) #find a tree edge to be replaced by nonbasic e
                i = B[1]
                i = i_basic[i]
                pred_temp[v] = nonbasic[δ_ind, 1]
                c_temp[e] = temp
                T_temp[i_nonbasic[δ_ind]] = 1
#                 println("Entering arc ", i_nonbasic[δ_ind])
#                 println("Exiting arc ", i)
                T_temp[i] = 0 ################
                Δ = Δ + δ_min
                delta = Δ
#                 println("current delta = ", delta)
                ##println("Pred = ", pred_temp)
            end
            
        end 
        #CHECK FOR LOOP:


        if y.*T_temp != y
            new_SP = true
        else
            loop = false
            if pred_temp[origin] != 0
                loop = true
            else
                v = edge[e,2]
                visit = collect(1:last_node)
                ##println("visit = ", visit)
                while v != origin && visit[v] != 0
                    ##println("v = ", v)
                    ##println("pred[v] = ", pred_temp[v] )
                    visit[v] = 0
                    v = pred_temp[v]

                    ##println("visit = ", visit)

                end

                if v != origin
                    loop = true
                    
                    ##println("HERE")
                end
            end
            if loop == true
                new_SP = true

            end
        end
        ##println("new_SP = ", new_SP)

    end

    return delta
end

function SA_basicSP(N, nonbasic, nonbasic_cost, label_temp, temp)
    delta = []
    delta_ind = []
    δ_ind = 0
    val_temp = 0
    if length(N) > 0
        ##println("nonbasic = ", nonbasic)
        for j in N
            head = find(nonbasic[:,2] .== j)
            if length(head) > 0
                ##println("head = ", head)
                for i in head
                    ##println("Arc #", i)
                    ##println("Arc name = ", nonbasic[i,:])
                    condd = nonbasic[i,1] in N
                    if condd == false
                        val_temp = label_temp[nonbasic[i,1]] + nonbasic_cost[i] - label_temp[nonbasic[i,2]]
                        push!(delta, val_temp)
                        push!(delta_ind, i)
                    end
                end
            end
        end
    end
    #IF DELTA DOESNT RETURN ANY VALUE AKA EMPTY: LET MIN DELTA = -1 THEN REMOVE IT LATER
    ###println("GOOD HERE")
    ##println("delta = ", delta)
    if length(delta) > 0
        for i = 1:length(delta)
            if i == 1
                δ_min = delta[1]
                δ_ind = delta_ind[1]
            end
            if delta[i] < δ_min
                δ_min = delta[i]
                δ_ind = delta_ind[i]
            end
        end
        temp = temp - δ_min
    else
        δ_min = -1.0
    end

    return δ_min, δ_ind, temp

end

function SA_basic(N, nonbasic, nonbasic_cost, label_temp, temp)
    delta = []
    delta_ind = []
    δ_ind = 0
    val_temp = 0
    if length(N) > 0
        ##println("nonbasic = ", nonbasic)
        for j in N
            tail = find(nonbasic[:,1] .== j)
            if length(tail) > 0
                ##println("tail = ", tail)
                for i in tail
                    ##println("Arc #", i)
                    ##println("Arc name = ", nonbasic[i,:])
                    condd = nonbasic[i,2] in N
                    if condd == false
                        ##println("FALSE THEN")
                        val_temp = label_temp[nonbasic[i,1]] + nonbasic_cost[i] - label_temp[nonbasic[i,2]]
                        push!(delta, val_temp)
                        push!(delta_ind, i)
                    end
                end
            end
        end
    end
    
    ##println("delta = ", delta)
    if length(delta) > 0
        ##println("DELTA > 0")
        for i = 1:length(delta)
            if i == 1
                δ_min = delta[1]
                δ_ind = delta_ind[1]
            end
            if delta[i] < δ_min
                δ_min = delta[i]
                δ_ind = delta_ind[i]
            end
        end
        temp = temp - δ_min
    else
        δ_min = -1.0
    end

    return δ_min, δ_ind, temp
end

function gx_bound(c_L, c_U, c, c_g, x_now, edge)
    
    #println(f,"Current CELL's LB = ", c_L)
    #println(f,"Current CELL's UB = ", c_U)
    #println(f, "Current interdiction x = ", x_now)
    start_node = edge[:,1]
    end_node = edge[:,2]

    no_node = max(maximum(start_node), maximum(end_node) )
    no_link = length(start_node)


    function getShortestX(state, start_node, end_node, origin, destination)
        _x = zeros(Int, length(start_node))
        _path = enumerate_paths(state, destination)

        for i=1:length(_path)-1
            _start = _path[i]
            _end = _path[i+1]

            for j=1:length(start_node)
                if start_node[j]==_start && end_node[j]==_end
                _x[j] = 1
                break
                end
            end

        end
        _x
    end


    graph = Graph(no_node)
    distmx = Inf*ones(no_node, no_node)

    # Adding links to the graph
    for i=1:no_link
        add_edge!(graph, start_node[i], end_node[i])
        distmx[start_node[i], end_node[i]] = c_g[i]
    end

    # Run Dijkstra's Algorithm from the origin node to all nodes
    state = dijkstra_shortest_paths(graph, origin, distmx)
    label = state.dists
    pred = state.parents
    b_arc = ""
    
    for i = 1: length(state.parents)
        if state.parents[i] != 0 
            b_arc = string(b_arc, "(", state.parents[i], ",", i, ")")
        end
    end
    
    # Retrieving the shortest path
    path = enumerate_paths(state, destination)
    
    #parents = LightGraphs.DijkstraState(state, destination)

    # Retrieving the 'x' variable in a 0-1 vector
    y = getShortestX(state, start_node, end_node, origin, destination)
    #println(f,"y vector:", y)
    
    gx = sum(c_g[i]*y[i] for i = 1:no_link)    
    SP = sum(c[i]*y[i] for i = 1:length(c))
    T = Int64[]
    for i = 1:Len
        if pred[edge[i,2]] == edge[i,1]
            push!(T, 1)
        else
            push!(T, 0)
        end
    end

    y_index = find(y .== 1)

    #println("Minimum Spanning Tree: ", T)
    #println("Shortest path y = ", y)
    #println("Indices of shortest path (edge) = ", y_index)
    #println("Nodes visited =", path)
    #println("Minimum Spanning Tree =", pred)
    #println("Node label =" ,label)

    return y, gx, SP, T, pred, label, path
end

function hx_bound(h, c_L, c_U, c, d, x_now, edge)
    
    c = (c_L + c_U)/2
    #println(f,"Current CELL's LB = ", c_L)
    #println(f,"Current CELL's UB = ", c_U)
    #println(f,"Current interdiction x = ", x_now)
    M = c_U - c_L

    h2 = copy(h1)

    y2 = h2[:y_h]

    @variable(h2, q[1:Len]>=0)
    for i = 1:Len
        @constraint(h2, q[i] >= c[i] - c_L[:,1][i] - M[i]*(1-y2[i])) #_h[i]))
    end

    @objective(h2, Min, sum((c_L[i]+d[i]*x_now[i])*y2[i] + q[i] for i=1:Len))

    #print(f,h2)
    solve(h2)
    hx = getobjectivevalue(h2)
    
    return getvalue(y2), hx
end
function corner_g(y, pred, c_L, c_U, d, x_now, edge, origin, destination)
    c_corner = zeros(length(edge[:,1]))
    
    for i = 1:length(edge[:,1])
        if y[i] == 1
            c_corner[i] = c_U[i] + x_now[i]*d[i]
        else
            c_corner[i] = c_L[i] + x_now[i]*d[i]
        end
    end
    π = zeros(length(edge[:,1])) + sum(c_U[k] + x_now[k]*d[k] for k = 1:length(edge[:,1])) + 1
    u = origin
    y_index = find(y.==1)
    π[u] = 0
    S = [u]
    while u != destination
        found_node = false
        count = 0
        while found_node == false 
            count = count + 1
            i = y_index[count]
            if edge[i,1] == u
                found_node = true
                v = edge[i,2]
                π[v] = π[u] + c_corner[i]
                u = v
                push!(S, u)
            end
        end
    end
    B = copy(S) + 0
    
    opt = true
    while opt == true && length(B) > 0
        u = B[1]
        B[1], B[length(B)] = B[length(B)], B[1]
        pop!(B)
        myArcs = find(edge[:,1].==u)
        if length(myArcs) > 0
            for i in myArcs
                if y[i] != 1
                    v = edge[i,2]
                    if π[v] > π[u] + c_corner[i]
                        π[v] = π[u] + c_corner[i]
                        if length(find(B.==v)) == 0
                            push!(B, v)
                        end
                        if length(find(S.==v)) > 0
                            opt = false
                        end
                    end
                end
            end
        else
            succ = find(pred.==u)
            for i = 1:length(succ)
                if length(find(B.==succ[i])) > 0
                    push!(B, succ[i])
                end
            end
        end
    end
    return opt, c_corner
end
# f = open(outfile, "w")
#MAIN PROGRAM:
df = DataFrame(CELL = Int[],ARCS_USED = Array[], LB = Array[], UB = Array[], PROB = Float64[], SP = Float64[])
# df2 = DataFrame(ITER = Int[], X = Array[])
# Cell_List = [1]
#=
weighted_unc = []
w_u = sum(cU_orig[i] - cL_orig[i] for i = 1:Len)
push!(weighted_unc, w_u)
weighted_unc_temp = copy(weighted_unc)
=#
last_x = zeros(Len)
# checkedGs = Int64[]
# checkedHs = Int64[]
checkedCells = Int64[]
x_count = 0
MP_obj = 0.0

push!(df, (1, yy, cL_orig, cU_orig, 1, SP_init))

##println(f,"MASTER PROBLEM==========================================================================================")

m = Model(solver = GurobiSolver(gurobi_env)) # If we want to add # in Gurobi, then we have to turn of Gurobi's own Cuts 

@variable(m,  x[1:Len], Bin)
#====For 15NodesPrelim Only ======#
#@variable(m, 1e6 >= z[1:200000] >= 0)
#@constraintref constr[1:200000]

@variable(m, 1e6 >= z[1:5000000] >= 0)
@constraintref constr[1:5000000]
@constraint(m, sum(x[i] for i=1:Len) <= b) #attack budget = 2
constr[1] = @constraint(m, z[1] <= SP_init + sum(yy[i]*x[i]*d[i] for i=1:Len) )
@objective(m, Max, sum(p[i]*z[i] for i = 1:length(p)) )#w - sum(s[k] for k=1:length(s))/length(s) )
##println(f,p)
##println(f,z[1:length(p)])

con_num = 1
stopping_cond = delta2 + 0
total_time = 0.0
iter = 0

while stopping_cond >= delta2 #===Both Timing & While Loop===#
    tic()
    stopping_cond = 0
    iter = iter + 1
#         println("\nIter ", iter)
        
    solve(m) 
#         println(m)
    MP = getvalue(sum(p[i]*z[i] for i = 1:length(p)))
    MP_obj = copy(MP)
    x_now = getvalue(x) + 0.0
#         println("Interdiction: ", find(x_now.==1) )
#         println("Number of Cells: ", length(p))
    split_Candidates = collect(1:length(p))
    if x_now != last_x 
        checkedCells = Int64[]
#             checkedGs = Int64[]
#             checkedHs = Int64[]
        x_count = 0
    else
        x_count = x_count + 1
    end

#         f = open(outfile, "a")
#             println(f,"\nIter : ", iter)
#             println(f, "Obj = ", MP_obj)
#             println(f, "Interdiction ", find(x_now.==1))
#             println(f, "# Cell = ", length(p))
#         close(f)
    z_now = getvalue(z)+ 0.0 
    k = 0
    myLength = length(p)

    if x_count >= 0 #!= 10
        zCutAdded = false  
#             println("checkedCells = ", checkedCells)
#         println("checkedGs = ", checkedGs)
        for k = 1:myLength #===1===#
            if (k in checkedCells) == false
#             if (k in checkedGs) == false
                row = find(df[:CELL] .== k)[end]
                c_L = df[:LB][row] 
                c_U = df[:UB][row] 
                y = df[:ARCS_USED][row]
                SP = df[:SP][row]
                c = (c_U + c_L)/2
                M = zeros(Len)
                for i = 1:Len
                    M[i] = c_U[i] - c_L[i]
                end
                c_g = c + d.*x_now
#BELONGS TO THE ORIGINAL CODE
                if x_now != last_x
                    y, gx, SP = gx_bound(c_L, c_U, c, c_g, x_now, edge)
                    g[k] = gx
                end
                if z_now[k] - g[k] > epsilon    
                    stopping_cond = delta2 + 1 
                    push!(df, (k, y, c_L + zeros(length(c_L)), c_U + zeros(length(c_U)), p[k], SP))
                    con_num = con_num + 1 
                    constr[con_num] = @constraint(m, z[k] <= sum(c[i]*y[i] + d[i]*y[i]*x[i] for i = 1:length(c_L)))
                    y_h, hx = hx_bound(h, c_L, c_U, c, d, x_now, edge)
                    h[k] = hx
                    if g[k]-h[k] <= delta2
                        push!(checkedCells, k)
                    end
                    if z_now[k] - g[k] > delta1                                
                        zCutAdded = true
                    end  
                end
            end
        end #===End of 1===#
#             println("zCutAdded = ", zCutAdded)
        ghCutAdded = false
        if zCutAdded == false #===2===#
            for k = 1:myLength #===2.1===#
                if (k in checkedCells) == false #===2.1.1===#                 
                    myLength = length(p)
                    row = find(df[:CELL] .== k)[1]
                    c_L = df[:LB][row] 
                    c_U = df[:UB][row] 

                    r_SA = 0.5*(c_U - c_L)
                    c = (c_U + c_L)/2
                    M = zeros(Len)
                    for i = 1:Len
                        M[i] = c_U[i] - c_L[i]
                    end
                    c_g = c + d.*x_now

                    y_h, hx = hx_bound(h, c_L, c_U, c, d, x_now, edge)
                    h[k] = hx

                    if g[k]-h[k] > delta2 #===2.1.1.1===#
#                                 println("k = ", k)
#                                 println("g[k] = ", g[k])
#                                 println("h[k] = ", h[k])
                    #SENSITIVITY ANALYSIS HERE
                        Delta_Final_Pos = ones(Len).*-1
                        Delta_Final_Neg = ones(Len).*-1

                        M_max = maximum(M)
                        temp = find(M .>0)
                        M_min = M_max
                        for e in temp
                            if M[e] <= M_min
                                M_min = M[e]
                            end
                        end
                        if M_max/M_min > 10
#                             println("Pick Max")
                            split_Candidates[k] = find(M .== maximum(M))[1]
                        end
                        if M_max/M_min <= 10
#                             println("SA")
                            #FIND TREE, PRED
                            y, gx, SP, T, pred, label, path = gx_bound(c_L, c_U, c, c_g, x_now, edge)
                            g[k] = gx
                            basic, basic_cost, i_basic, nonbasic, nonbasic_cost, i_nonbasic = classify_edges_costs(T,c_g)
                            orderedArcList = getSortedBasicArcList(pred, basic, i_basic)

                            label_temp = copy(label)
                            c_temp = copy(c_g)

                            for e in i_nonbasic
                                u = edge[e,1]
                                v = edge[e,2]
                                Delta_Final_Neg[e] = label_temp[u] - label_temp[v] + c_temp[e]
                            end

                            for i = 1 : length(orderedArcList)
                                e = orderedArcList[i]
                                myNode = edge[e,2]
                                if y[e] == 1 #if arc is on the shortest path y
                                    Delta_Final_Pos[e] = find_basic_Delta(myNode, Delta_Final_Pos, Delta_Final_Neg, "POS", basic, i_basic, nonbasic, i_nonbasic) 
                                    #Do not have to find Delta_Final_Neg : all set to (-1)
                                else
                                    Delta_Final_Neg[e] = find_basic_Delta(myNode, Delta_Final_Pos, Delta_Final_Neg, "NEG", basic, i_basic, nonbasic, i_nonbasic)
                                end
                            end          

                            arcs_neg = find(Delta_Final_Neg .>= 0)
                            arcs_pos = find(Delta_Final_Pos .>= 0)
                            found_split_arc = false
                            arcs_delta = []
                            myMax = maximum(Delta_Final_Neg)

                            arcs = find(M.>0)                        
                            if myMax < maximum(Delta_Final_Pos)
                                myMax = maximum(Delta_Final_Pos)
                            end

                            myScore = myMax/M_min
                            for i in arcs
#                                 println("i = ", i)
                                if Delta_Final_Pos[i] >= 0 && Delta_Final_Pos[i]/M[i] < myScore
                                    myScore = Delta_Final_Pos[i]/M[i]
                                    split_Candidates[k] = i
                                    found_split_arc = true
                                end
                                if Delta_Final_Neg[i] >= 0 && Delta_Final_Neg[i]/M[i] < myScore
                                    myScore = Delta_Final_Neg[i]/M[i]
                                    split_Candidates[k] = i
                                    found_split_arc = true
                                end
                            end
                            if found_split_arc == false
                                split_Candidates[k] = find(M .== maximum(M))[1]
                            end
                        end                        
                                                         
                        stopping_cond = delta2 + 1
                        found = true
                        arc_split = split_Candidates[k]
#                                 arc_split = find(M .== maximum(M))[1]
#                                     print(" arc split ", arc_split)
                        gap = (c_U[arc_split] - c[arc_split])/2
                        newCell = maximum(df[:CELL])+1

                        p_temp = p[k]/2
                        p[k] = p_temp
                        push!(p,p_temp)

                        last_row = length(df[:CELL])

                        c_Mid_k = copy(c_U)
                        c_Mid_k[arc_split] = (c_U[arc_split]+c_L[arc_split])/2
                        c_Mid_newCell = copy(c_L)
                        c_Mid_newCell[arc_split] = c_Mid_k[arc_split]

                        c_Ltemp = (c_L + c_Mid_k)/2
                        c_g_Ltemp = c_Ltemp + d.*x_now

                        c_Utemp = (c_U + c_Mid_newCell)/2
                        c_g_Utemp = c_Utemp + d.*x_now

                        y_low, gx_low, SP_low = gx_bound(c_L, c_Mid_k, c_Ltemp, c_g_Ltemp, x_now, edge)
                        y_high, gx_high, SP_high = gx_bound(c_Mid_newCell, c_U, c_Utemp, c_g_Utemp, x_now, edge)
                        g[k] =  gx_low
                        push!(g, gx_high)
                        push!(h, 0)
#                                     println("y_low ", find(y_low.==1))
#                                     println("y_high ", find(y_high.==1))
                        y_low_exists = false
                        y_high_exists = false

                        constr_of_k = find(df[:CELL].== k)

                        for i in constr_of_k 
                            df[:PROB][i] = p_temp #Probability must change
                            ARCS_USED_i = df[:ARCS_USED][i] #Get path y in that row/constraint
                            newCell_RHS = df[:SP][i] #Get a temporary RHS for newCell
                            df[:UB][i] = c_Mid_k
                            #ONLY UPDATE RHS IF ARC SPLIT IS ON THAT PATH           
                            if arc_split in find(ARCS_USED_i.==1)
                                k_new_RHS = df[:SP][i] - gap #UPDATE RHS OF CELL k
                                newCell_RHS = df[:SP][i] + gap #UPDATE RHS OF NEWCELL
                                df[:SP][i] = k_new_RHS #FIX DF INFORMATION
                                JuMP.setRHS(constr[i], k_new_RHS)
                            end

                            #COPY THE CONSTRAINT FROM k to NEWCELL
                            con_num = con_num + 1
                            constr[con_num] = @constraint(m, z[newCell] <= 
                    sum(d[i]*x[i]*ARCS_USED_i[i] for i = 1:length(c_L)) + newCell_RHS)

                            push!(df, (newCell, ARCS_USED_i, c_Mid_newCell, c_U, p_temp, newCell_RHS))

                            if ARCS_USED_i == y_low
                                y_low_exists = true
                            end
                            if ARCS_USED_i == y_high
                                y_high_exists = true
                            end

                        end

                        if y_low_exists == false
                            push!(df, (k, y_low, c_L, c_Mid_k, p_temp, SP_low)) 
                            con_num = con_num + 1
                            constr[con_num] = @constraint(m, z[k] <= 
                            sum(d[i]*x[i]*y_low[i] for i = 1:length(c_L)) + SP_low)
                        end
                        if y_high_exists == false
                            push!(df, (newCell, y_high, c_Mid_newCell, c_U, p_temp, SP_high)) 
                            con_num = con_num + 1
                            constr[con_num] = @constraint(m, z[newCell] <= 
                            sum(d[i]*x[i]*y_high[i] for i = 1:length(c_L)) + SP_high)
                        end
                        ghCutAdded = true
                    else
                        push!(checkedCells, k)
                    end #===End of 2.1.1.1===#

#                         end #End of if rule100 > 1.0
                end #===End of 2.1.1===#
            end #===End of 2.1===#
            @objective(m, Max, sum(p[i]*z[i] for i = 1:length(p)) )
        end #===End of 2===#
    end
#     println("x_now end "= )
#         println("z = ", z_now[1:length(p)])
#         println("g = ", g)
#         println("h = ", h)
    time_lapse = toq()    
    total_time = total_time + time_lapse
    last_x = x_now
#         if total_time > 3600
#             total_time = ">3600"
# #             println("total time" > 3600")
#             break
#         end
end #===End While Loop===#


opt_gap = MP_obj - sum(p[i]*(h[i]) for i = 1:length(p))
# @timeit to "CALC OPT GAP" opt_gap = sum(p[i]*(g[i] - h[i]) for i = 1:length(p))

# f = open(outfile, "w")
# println(f,"Optimality gap = ", opt_gap)
# println(f,"MP_obj = ", MP_obj)
# println(f,"g = ", sum(p[i]*(g[i]) for i = 1:length(p)))
# println(f,"h = ", sum(p[i]*(h[i]) for i = 1:length(p)))
# println(f,"Overall runtime = ", total_time)
# println(f,"No. iterations = ", iter) #length(df[:CELL]))
temp_g = sum(p[i]*g[i] for i = 1:length(p))
temp_h = sum(p[i]*h[i] for i = 1:length(p))
timesFile = open("/Users/din/Documents/July2020Paper1/Summary/local_MulticutMultipartitionScore.txt", "a")
println("Laptop ", Grp, "; laptopIns ", myInstance, "; Time ", total_time, "; MP ", MP_obj, "; g ", temp_g,"; h ",temp_h,"; Opt ",opt_gap,"; Cells ", length(p))
println(find(last_x.==1))
close(timesFile)



