
#Ready for Upload

using JuMP 
using Gurobi
using LightGraphs
using DataFrames
using CSV
using TimerOutputs

#myInstance = string(parse(Int64,ARGS[1]))
myInstance = "2"
myRun = Dates.format(now(), "HH:MM:SS")
epsilon = 0.001
outfile = "./Output/Multicut"*myInstance*"_"*myRun*".txt"
gurobi_env = Gurobi.Env()

to = TimerOutput()
myFile = "./Instances/15NodesPrelim_"*myInstance*".jl"
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

#     println("Minimum Spanning Tree: ", T)
    #println("Shortest path y = ", y)
    #println("Indices of shortest path (edge) = ", y_index)
    #println("Nodes visited =", path)
    #println("Minimum Spanning Tree =", pred)
    #println("Node label =" ,label)

    return y, gx, SP
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

MP_obj = 0.0

push!(df, (1, yy, cL_orig, cU_orig, 1, SP_init))

##println(f,"MASTER PROBLEM==========================================================================================")

m = Model(solver = GurobiSolver(gurobi_env)) # If we want to add # in Gurobi, then we have to turn of Gurobi's own Cuts 

@variable(m,  x[1:Len], Bin)
@variable(m, 1e6 >= z[1:200000] >= 0)
@constraintref constr[1:200000]
@constraint(m, sum(x[i] for i=1:Len) <= b) #attack budget = 2
constr[1] = @constraint(m, z[1] <= SP_init + sum(yy[i]*x[i]*d[i] for i=1:Len) )
@objective(m, Max, sum(p[i]*z[i] for i = 1:length(p)) )#w - sum(s[k] for k=1:length(s))/length(s) )
##println(f,p)
##println(f,z[1:length(p)])

con_num = 1
stopping_cond = delta2 + 0
total_time = 0.0
iter = 0
while stopping_cond >= delta2 
    tic()
    stopping_cond = 0
    iter = iter + 1
#         println(f,"\n\nITERATION ",length(df[:CELL]),": ========================================================================================== ")
#         println("\nIter : ", iter)
#         println(f, "p = ", p)
#         println(f,"Cell List = ", Cell_List)
#         if iter == 90
#             println(m)
#         end
    solve(m) 

    MP = getvalue(sum(p[i]*z[i] for i = 1:length(p)))
    MP_obj = copy(MP)

    x_now = getvalue(x) + 0.0    
    last_x = x_now

    z_now = getvalue(z)+ 0.0 
    gap = 0.0
    k = 0
    found = false

    myLength = length(p)
    zCutAdded = false
    for k = 1:myLength

        row = find(df[:CELL] .== k)[1]
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
        y, gx, SP = gx_bound(c_L, c_U, c, c_g, x_now, edge)
        g[k] = gx

        if z_now[k] - g[k] > epsilon    
            stopping_cond = delta2 + 1 
            push!(df, (k, y, c_L + zeros(length(c_L)), c_U + zeros(length(c_U)), p[k], SP))
            con_num = con_num + 1 
            constr[con_num] = @constraint(m, z[k] <= sum(c[i]*y[i] + d[i]*y[i]*x[i] for i = 1:length(c_L)))
            if z_now[k] - g[k] > delta1                             
                zCutAdded = true
            end  
        end
    end
    @objective(m, Max, sum(p[i]*z[i] for i = 1:length(p)) )
        
    if zCutAdded == false
#             found = false
        k=0
        while found == false && k < myLength
            k=k+1
            row = find(df[:CELL] .== k)[1]
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

            y_h, hx = hx_bound(h, c_L, c_U, c, d, x_now, edge)

            ##println(f,"k = ", k)
            if length(h) >= k 
                h[k] = hx
            else
                h = push!(h,hx)
            end
     
            if g[k]-h[k] > delta2
#                     println("Cell ", k)
#                     println("g[k] = ", g[k])
#                     println("h[k] = ", h[k])
#                    println("gx hx")
#                     println("Add g-h cut at cell ", k)   
                stopping_cond = delta2 + 1
                found = true
                arc_split = find(M .== maximum(M))[1]
                gap = (c_U[arc_split] - c[arc_split])/2
                newCell = maximum(df[:CELL])+1

                p_temp = p[k]/2
                p[k] = p_temp

                push!(p,p_temp)
#                         push!(Cell_List, newCell)
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

                y_low_exists = false
                y_high_exists = false

                constr_of_k = find(df[:CELL].== k)

                for i in constr_of_k #LOOP THRU ALL ROWS IN DF ASSOCIATED WITH k
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
                        
            end #END OF if gx-hx > delta2
        end #END OF while found == false && k < myLength
            @objective(m, Max, sum(p[i]*z[i] for i = 1:length(p)) )
    end #if zCutAdded == false
    time_lapse = toq()    
    total_time = total_time + time_lapse

end #while stopping_cond >= delta2 
 
# print(m)
# CSV.write("/Users/din/Desktop/Research_S2019/Base_V2_"*myInstance*".cvs", df)  
##CSV.write("/Users/din/Desktop/Research/InterdicionX_Ins1_BaseClever", df2) 

opt_gap =MP_obj - sum(p[i]*(h[i]) for i = 1:length(p))
# f = open(outfile, "w")
# println(f,"Optimality gap = ", opt_gap)
# println(f,"MP_obj = ", MP_obj)
# println(f,"g = ", sum(p[i]*(g[i]) for i = 1:length(p)))
# println(f,"h = ", sum(p[i]*(h[i]) for i = 1:length(p)))
# println(f,"Overall runtime = ", total_time)
# println(f,"No. iterations = ", iter) #length(df[:CELL]))
temp_g = sum(p[i]*g[i] for i = 1:length(p))
    temp_h = sum(p[i]*h[i] for i = 1:length(p))
timesFile = open("/Users/din/Documents/July2020Paper1/Summary/local_Multicut.txt", "a")
println(timesFile, "Laptop ", Grp, "; laptopIns ", myInstance, "; Time ", total_time, "; MP ", MP_obj, "; g ", temp_g,"; h ",temp_h,"; Opt ",opt_gap,"; Cells ", length(p))
close(timesFile)

#show(to; allocations = false)


