edge = [1 2; 1 7; 1 9; 1 11; 2 6; 2 19; 3 4; 3 15; 4 18; 4 19; 5 15; 6 5; 6 13; 6 17; 7 10; 8 4; 8 7; 8 15; 8 19; 9 6; 9 12; 10 7; 10 12; 10 13; 11 3; 11 5; 11 7; 11 8; 11 14; 11 15; 11 20; 12 15; 13 14; 13 16; 13 17; 14 3; 14 10; 14 16; 15 2; 15 7; 15 9; 16 3; 16 4; 16 11; 16 17; 17 2; 17 10; 17 11; 17 14; 17 19; 18 2; 18 4; 18 5; 18 8; 18 12; 18 14; 19 6; 19 9; 19 15]
cL_orig = [0.0, 63.0, 69.0, 100.0, 39.0, 174.0, 0.0, 115.0, 132.0, 147.0, 95.0, 5.0, 67.0, 112.0, 26.0, 39.0, 13.0, 73.0, 115.0, 34.0, 25.0, 18.0, 12.0, 35.0, 83.0, 58.0, 37.0, 25.0, 29.0, 45.0, 94.0, 22.0, 15.0, 27.0, 40.0, 112.0, 43.0, 23.0, 131.0, 82.0, 55.0, 133.0, 123.0, 44.0, 11.0, 146.0, 69.0, 58.0, 33.0, 24.0, 159.0, 144.0, 132.0, 104.0, 55.0, 42.0, 128.0, 99.0, 45.0]
cU_orig = [28.0, 63.0, 95.0, 100.0, 39.0, 174.0, 26.0, 115.0, 146.0, 147.0, 95.0, 5.0, 67.0, 112.0, 26.0, 39.0, 13.0, 73.0, 115.0, 34.0, 45.0, 48.0, 28.0, 35.0, 83.0, 58.0, 37.0, 25.0, 29.0, 45.0, 94.0, 40.0, 15.0, 27.0, 40.0, 112.0, 43.0, 23.0, 131.0, 82.0, 55.0, 133.0, 123.0, 64.0, 11.0, 146.0, 69.0, 58.0, 33.0, 24.0, 159.0, 144.0, 132.0, 104.0, 65.0, 42.0, 128.0, 99.0, 45.0]
d = [10.0, 11.0, 8.0, 2.0, 4.0, 9.0, 17.0, 19.0, 8.0, 8.0, 3.0, 14.0, 19.0, 12.0, 2.0, 6.0, 6.0, 12.0, 4.0, 16.0, 12.0, 13.0, 5.0, 10.0, 20.0, 8.0, 10.0, 6.0, 18.0, 20.0, 20.0, 13.0, 19.0, 17.0, 8.0, 1.0, 13.0, 19.0, 15.0, 13.0, 20.0, 1.0, 20.0, 14.0, 17.0, 4.0, 15.0, 4.0, 5.0, 18.0, 12.0, 6.0, 6.0, 2.0, 1.0, 3.0, 3.0, 8.0, 8.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
Grp = center
last_node = maximum(edge)
