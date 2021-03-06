edge = [1 2; 1 5; 1 12; 1 13; 2 9; 2 13; 2 15; 3 2; 3 5; 3 15; 4 16; 5 10; 5 12; 5 18; 7 18; 8 12; 8 15; 8 18; 9 2; 9 5; 9 8; 10 3; 10 8; 10 9; 10 18; 12 2; 12 5; 12 8; 12 10; 13 4; 13 9; 13 12; 13 20; 15 2; 15 7; 15 16; 16 4; 16 5; 16 7; 16 13; 16 20; 18 2; 18 5; 18 9]
cL_orig = Any[6.0, 44.0, 105.0, 124.0, 64.0, 108.0, 122.0, 13.0, 17.0, 118.0, 119.0, 42.0, 63.0, 130.0, 102.0, 41.0, 73.0, 101.0, 73.0, 33.0, 12.0, 69.0, 18.0, 9.0, 84.0, 100.0, 64.0, 38.0, 10.0, 89.0, 37.0, 5.0, 67.0, 135.0, 74.0, 14.0, 120.0, 100.0, 88.0, 35.0, 31.0, 163.0, 133.0, 95.0]
cU_orig = Any[20.0, 44.0, 105.0, 124.0, 76.0, 122.0, 142.0, 13.0, 33.0, 118.0, 119.0, 60.0, 75.0, 130.0, 118.0, 41.0, 73.0, 101.0, 73.0, 45.0, 12.0, 79.0, 18.0, 9.0, 84.0, 100.0, 78.0, 38.0, 28.0, 89.0, 37.0, 5.0, 67.0, 135.0, 86.0, 14.0, 120.0, 116.0, 88.0, 35.0, 41.0, 163.0, 133.0, 95.0]
d = Any[10.0, 8.0, 17.0, 7.0, 19.0, 17.0, 12.0, 18.0, 19.0, 3.0, 20.0, 5.0, 2.0, 9.0, 14.0, 20.0, 15.0, 3.0, 17.0, 6.0, 16.0, 10.0, 12.0, 18.0, 18.0, 9.0, 1.0, 13.0, 15.0, 18.0, 15.0, 5.0, 20.0, 4.0, 18.0, 6.0, 8.0, 10.0, 1.0, 9.0, 18.0, 4.0, 2.0, 15.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)
p = [1.0]
g = [SP_init]
h = [0.0]
origin = 1
destination = 20
last_node = maximum(edge)
all_nodes = collect(1:last_node)
M_orig = zeros(Len)
for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end
delta1 = 1.0
delta2 = 3.0
b = 2
last_node = maximum(edge)
Grp = "-5Nodes"
