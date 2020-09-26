edge = [1 2; 1 11; 1 14; 1 19; 2 3; 2 19; 3 6; 4 2; 4 3; 4 9; 4 11; 4 14; 4 17; 5 3; 5 16; 6 9; 6 17; 6 20; 7 4; 7 13; 8 5; 8 14; 8 18; 9 2; 9 11; 9 13; 9 17; 10 12; 10 13; 10 19; 11 8; 11 14; 11 15; 12 9; 12 13; 12 19; 13 4; 13 8; 13 12; 13 14; 13 17; 13 19; 14 9; 14 13; 14 20; 15 4; 15 9; 15 16; 15 17; 15 18; 16 5; 16 14; 16 15; 17 3; 17 7; 17 10; 17 14; 18 19; 19 12; 19 13; 19 14; 19 16; 19 17; 19 18]
cL_orig = [10.0, 98.5, 127.0, 174.0, 10.0, 168.0, 34.0, 20.0, 11.0, 52.0, 67.0, 89.0, 125.0, 23.0, 105.0, 32.0, 101.0, 137.0, 30.0, 48.5, 22.0, 63.0, 100.0, 70.0, 24.0, 35.0, 81.0, 21.0, 33.0, 94.0, 23.5, 25.0, 37.0, 33.0, 10.0, 71.0, 90.0, 48.0, 13.0, 9.0, 39.0, 56.0, 47.0, 15.0, 59.5, 107.0, 64.0, 5.0, 18.0, 31.0, 105.0, 24.0, 9.0, 145.0, 95.5, 71.0, 21.5, 5.0, 62.0, 60.0, 55.0, 32.0, 14.0, 15.0]
cU_orig = [10.0, 109.5, 127.0, 182.0, 10.0, 168.0, 34.0, 20.0, 11.0, 52.0, 67.0, 103.0, 125.0, 23.0, 105.0, 32.0, 111.0, 137.0, 30.0, 65.5, 32.0, 63.0, 100.0, 76.0, 24.0, 45.0, 81.0, 21.0, 33.0, 94.0, 32.5, 29.0, 37.0, 33.0, 10.0, 71.0, 100.0, 48.0, 17.0, 9.0, 39.0, 66.0, 47.0, 15.0, 68.5, 107.0, 64.0, 5.0, 18.0, 31.0, 105.0, 24.0, 11.0, 145.0, 106.5, 71.0, 34.5, 5.0, 72.0, 60.0, 55.0, 32.0, 22.0, 15.0]
d = [16.0, 13.0, 4.0, 1.0, 6.0, 1.0, 8.0, 10.0, 18.0, 3.0, 15.0, 12.0, 8.0, 8.0, 19.0, 14.0, 19.0, 11.0, 4.0, 8.0, 17.0, 14.0, 5.0, 3.0, 5.0, 5.0, 9.0, 15.0, 7.0, 16.0, 20.0, 14.0, 17.0, 14.0, 20.0, 14.0, 17.0, 11.0, 18.0, 20.0, 5.0, 11.0, 4.0, 13.0, 14.0, 11.0, 19.0, 2.0, 13.0, 13.0, 2.0, 10.0, 4.0, 4.0, 14.0, 7.0, 19.0, 19.0, 13.0, 13.0, 13.0, 1.0, 11.0, 18.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
Grp = "+10%Uncertainties"