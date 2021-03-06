edge = [1 4; 1 5; 1 12; 1 14; 1 15; 1 19; 2 9; 3 10; 3 12; 4 10; 4 11; 4 15; 5 12; 5 14; 6 2; 6 5; 6 18; 7 4; 8 7; 8 9; 8 15; 9 5; 9 13; 10 7; 10 20; 11 3; 11 8; 11 9; 11 10; 11 12; 11 13; 12 7; 12 10; 12 14; 12 16; 12 17; 13 2; 13 6; 13 11; 14 10; 14 15; 15 5; 15 8; 15 11; 15 20; 16 12; 16 15; 16 19; 17 15; 18 6; 19 10; 19 15]
cL_orig = Any[35.0, 39.0, 109.0, 128.0, 129.0, 179.0, 69.0, 75.0, 88.0, 49.0, 62.0, 105.0, 65.0, 91.0, 41.0, 3.0, 118.0, 32.0, 10.0, 5.0, 72.0, 41.0, 32.0, 35.0, 89.0, 75.0, 22.0, 22.0, 8.0, 14.0, 12.0, 52.0, 16.0, 21.0, 39.0, 45.0, 106.0, 69.0, 23.0, 28.0, 12.0, 104.0, 71.0, 29.0, 39.0, 37.0, 8.0, 34.0, 25.0, 115.0, 95.0, 37.0]
cU_orig = Any[35.0, 39.0, 109.0, 138.0, 147.0, 179.0, 69.0, 75.0, 88.0, 65.0, 82.0, 105.0, 65.0, 91.0, 41.0, 23.0, 128.0, 32.0, 10.0, 19.0, 72.0, 41.0, 44.0, 35.0, 109.0, 89.0, 32.0, 22.0, 8.0, 14.0, 30.0, 52.0, 16.0, 21.0, 39.0, 45.0, 122.0, 69.0, 23.0, 42.0, 12.0, 104.0, 71.0, 47.0, 51.0, 37.0, 8.0, 34.0, 25.0, 125.0, 95.0, 37.0]
d = Any[14.0, 8.0, 12.0, 5.0, 15.0, 10.0, 16.0, 12.0, 20.0, 7.0, 1.0, 3.0, 3.0, 19.0, 14.0, 14.0, 7.0, 16.0, 2.0, 13.0, 7.0, 7.0, 17.0, 14.0, 7.0, 11.0, 1.0, 19.0, 3.0, 18.0, 6.0, 13.0, 13.0, 7.0, 5.0, 5.0, 3.0, 3.0, 7.0, 14.0, 16.0, 11.0, 13.0, 8.0, 2.0, 16.0, 14.0, 8.0, 9.0, 18.0, 3.0, 11.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
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
delta2 = 5.0
b = 2
last_node = maximum(edge)
Grp = "0.1Density"
