edge = [1 15; 1 19; 2 14; 2 17; 3 6; 3 20; 4 5; 4 10; 5 2; 5 9; 5 19; 6 3; 6 11; 7 13; 8 7; 8 11; 8 16; 8 18; 9 2; 9 12; 10 4; 10 14; 11 13; 12 5; 12 6; 12 14; 13 5; 13 14; 14 4; 14 17; 14 20; 15 11; 16 2; 16 3; 16 8; 16 9; 17 3; 17 4; 17 13; 17 15; 18 13; 19 7; 19 8; 19 18]
cL_orig = Any[137.0, 184.0, 115.0, 146.0, 22.0, 165.0, 7.0, 63.0, 26.0, 31.0, 136.0, 32.0, 38.0, 65.0, 5.0, 30.0, 74.0, 99.0, 75.0, 28.0, 46.0, 45.0, 12.0, 62.0, 56.0, 24.0, 73.0, 12.0, 104.0, 28.0, 57.0, 36.0, 145.0, 119.0, 75.0, 66.0, 142.0, 133.0, 43.0, 17.0, 50.0, 121.0, 113.0, 5.0]
cU_orig = Any[151.0, 184.0, 115.0, 146.0, 38.0, 179.0, 21.0, 63.0, 26.0, 41.0, 136.0, 32.0, 54.0, 65.0, 5.0, 30.0, 88.0, 99.0, 75.0, 40.0, 66.0, 45.0, 22.0, 74.0, 70.0, 24.0, 87.0, 12.0, 104.0, 28.0, 57.0, 36.0, 145.0, 133.0, 75.0, 66.0, 142.0, 133.0, 43.0, 17.0, 50.0, 121.0, 113.0, 5.0]
d = Any[8.0, 8.0, 8.0, 1.0, 16.0, 16.0, 10.0, 10.0, 14.0, 17.0, 20.0, 15.0, 8.0, 10.0, 1.0, 20.0, 11.0, 4.0, 10.0, 6.0, 18.0, 16.0, 7.0, 12.0, 14.0, 1.0, 16.0, 5.0, 11.0, 9.0, 7.0, 8.0, 2.0, 15.0, 3.0, 4.0, 13.0, 6.0, 3.0, 18.0, 19.0, 3.0, 7.0, 4.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
delta2 = 4.5
b = 2
last_node = maximum(edge)
Grp = "0.1Density"
