edge = [1 5; 1 9; 1 15; 1 19; 2 14; 2 17; 2 20; 3 6; 3 14; 3 18; 3 20; 4 5; 4 10; 4 20; 5 2; 5 9; 5 19; 6 3; 6 11; 7 5; 7 8; 7 11; 7 13; 7 17; 8 6; 8 7; 8 11; 8 16; 8 18; 9 2; 9 12; 9 14; 10 3; 10 4; 10 8; 10 14; 11 6; 11 13; 11 19; 12 5; 12 6; 12 11; 12 14; 12 17; 13 5; 13 14; 13 16; 13 20; 14 3; 14 4; 14 17; 14 20; 15 11; 15 12; 16 2; 16 3; 16 5; 16 8; 16 9; 17 3; 17 4; 17 5; 17 8; 17 13; 17 15; 18 11; 18 13; 19 7; 19 8; 19 18]
cL_orig = [32.0, 74.0, 137.0, 184.0, 115.0, 146.0, 174.0, 30.0, 105.0, 154.0, 172.0, 7.0, 63.0, 162.0, 26.0, 31.0, 136.0, 32.0, 46.0, 25.0, 6.0, 35.0, 65.0, 101.0, 20.0, 5.0, 30.0, 81.0, 99.0, 75.0, 28.0, 48.0, 72.0, 46.0, 17.0, 45.0, 48.0, 12.0, 76.0, 68.0, 56.0, 8.0, 24.0, 48.0, 80.0, 12.0, 32.0, 59.0, 107.0, 104.0, 28.0, 57.0, 36.0, 31.0, 145.0, 119.0, 111.0, 75.0, 66.0, 142.0, 133.0, 125.0, 85.0, 43.0, 17.0, 69.0, 50.0, 121.0, 113.0, 5.0]
cU_orig = [46.0, 90.0, 151.0, 184.0, 115.0, 146.0, 188.0, 30.0, 105.0, 154.0, 172.0, 21.0, 63.0, 162.0, 26.0, 41.0, 136.0, 32.0, 46.0, 25.0, 20.0, 35.0, 65.0, 101.0, 20.0, 5.0, 30.0, 81.0, 99.0, 75.0, 40.0, 48.0, 72.0, 66.0, 17.0, 45.0, 48.0, 22.0, 76.0, 68.0, 70.0, 8.0, 24.0, 48.0, 80.0, 12.0, 32.0, 75.0, 119.0, 104.0, 28.0, 57.0, 36.0, 31.0, 145.0, 133.0, 111.0, 75.0, 66.0, 142.0, 133.0, 125.0, 85.0, 43.0, 17.0, 69.0, 50.0, 121.0, 113.0, 5.0]
d = [14.0, 20.0, 8.0, 8.0, 8.0, 1.0, 14.0, 16.0, 17.0, 18.0, 16.0, 10.0, 10.0, 16.0, 14.0, 17.0, 20.0, 15.0, 8.0, 14.0, 2.0, 7.0, 10.0, 8.0, 7.0, 1.0, 20.0, 11.0, 4.0, 10.0, 6.0, 15.0, 20.0, 18.0, 13.0, 16.0, 20.0, 7.0, 12.0, 12.0, 14.0, 18.0, 1.0, 7.0, 16.0, 5.0, 16.0, 10.0, 19.0, 11.0, 9.0, 7.0, 8.0, 15.0, 2.0, 15.0, 8.0, 3.0, 4.0, 13.0, 6.0, 4.0, 18.0, 3.0, 18.0, 6.0, 19.0, 3.0, 7.0, 4.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
delta1 = 2.0
delta2 = 4.5
b = 2
last_node = maximum(edge)
Grp = "+1delta1"
