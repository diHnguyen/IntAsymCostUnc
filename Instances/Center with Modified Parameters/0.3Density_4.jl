edge = [1 3; 1 18; 2 4; 2 6; 2 14; 2 15; 2 16; 2 17; 3 5; 3 6; 3 7; 3 16; 4 2; 4 7; 4 8; 4 12; 4 15; 5 9; 5 13; 5 14; 5 19; 6 2; 6 11; 7 2; 7 12; 7 15; 8 6; 8 10; 8 11; 8 18; 9 6; 9 10; 9 16; 9 20; 10 11; 10 15; 11 17; 11 20; 12 4; 12 9; 12 13; 12 20; 13 18; 14 15; 15 9; 15 12; 15 13; 15 14; 16 3; 16 12; 16 13; 17 12; 17 14; 17 15; 18 2; 18 3; 18 12; 18 13; 19 3; 19 5; 19 7; 19 9; 1 2; 1 6; 2 5; 3 15; 4 10; 5 17; 6 3; 6 5; 6 7; 6 10; 6 12; 6 14; 8 2; 8 7; 8 14; 8 16; 9 14; 9 17; 9 19; 10 5; 11 4; 11 18; 12 10; 12 11; 13 3; 13 9; 13 19; 14 2; 14 6; 14 11; 15 7; 15 19; 16 9; 17 9; 17 10; 18 5; 18 6; 19 17; 19 20]
cL_orig = [21.0, 165.0, 20.0, 34.0, 116.0, 126.0, 131.0, 148.0, 6.0, 22.0, 40.0, 127.0, 22.0, 28.0, 43.0, 75.0, 100.0, 40.0, 76.0, 90.0, 139.0, 40.0, 46.0, 46.0, 38.0, 77.0, 16.0, 21.0, 35.0, 104.0, 31.0, 15.0, 68.0, 98.0, 10.0, 51.0, 53.0, 86.0, 85.0, 31.0, 12.0, 82.0, 49.0, 14.0, 57.0, 34.0, 20.0, 12.0, 134.0, 39.0, 32.0, 48.0, 30.0, 21.0, 145.0, 136.0, 63.0, 39.0, 165.0, 138.0, 125.0, 103.0, 8.0, 48.0, 26.0, 120.0, 56.0, 125.0, 32.0, 15.0, 7.0, 43.0, 64.0, 76.0, 59.0, 5.0, 60.0, 76.0, 47.0, 85.0, 101.0, 50.0, 75.0, 72.0, 23.0, 9.0, 101.0, 41.0, 63.0, 124.0, 79.0, 27.0, 81.0, 40.0, 72.0, 80.0, 65.0, 133.0, 125.0, 20.0, 9.0]
cU_orig = [21.0, 165.0, 20.0, 46.0, 116.0, 126.0, 147.0, 148.0, 26.0, 40.0, 40.0, 127.0, 22.0, 28.0, 43.0, 95.0, 118.0, 40.0, 76.0, 90.0, 139.0, 40.0, 46.0, 46.0, 56.0, 77.0, 16.0, 21.0, 35.0, 104.0, 31.0, 15.0, 68.0, 114.0, 10.0, 51.0, 63.0, 86.0, 85.0, 31.0, 12.0, 82.0, 49.0, 14.0, 57.0, 34.0, 20.0, 12.0, 134.0, 39.0, 32.0, 48.0, 30.0, 21.0, 165.0, 154.0, 63.0, 51.0, 165.0, 138.0, 125.0, 103.0, 8.0, 48.0, 26.0, 120.0, 56.0, 125.0, 32.0, 15.0, 7.0, 43.0, 64.0, 76.0, 59.0, 5.0, 60.0, 76.0, 47.0, 85.0, 101.0, 50.0, 75.0, 72.0, 23.0, 9.0, 101.0, 41.0, 63.0, 124.0, 79.0, 27.0, 81.0, 40.0, 72.0, 80.0, 65.0, 133.0, 125.0, 20.0, 9.0]
d = [16.0, 11.0, 20.0, 5.0, 12.0, 19.0, 12.0, 5.0, 16.0, 17.0, 8.0, 4.0, 9.0, 15.0, 20.0, 15.0, 10.0, 3.0, 9.0, 2.0, 15.0, 1.0, 3.0, 19.0, 4.0, 5.0, 18.0, 20.0, 4.0, 10.0, 19.0, 4.0, 20.0, 13.0, 11.0, 12.0, 13.0, 19.0, 5.0, 4.0, 11.0, 13.0, 17.0, 13.0, 6.0, 9.0, 15.0, 18.0, 13.0, 9.0, 8.0, 5.0, 20.0, 8.0, 10.0, 2.0, 10.0, 4.0, 11.0, 6.0, 2.0, 10.0, 13.0, 10.0, 10.0, 18.0, 16.0, 9.0, 15.0, 13.0, 5.0, 16.0, 17.0, 5.0, 12.0, 11.0, 14.0, 5.0, 4.0, 1.0, 18.0, 11.0, 8.0, 15.0, 16.0, 2.0, 2.0, 8.0, 16.0, 18.0, 9.0, 11.0, 10.0, 10.0, 17.0, 9.0, 13.0, 17.0, 12.0, 8.0, 15.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
delta2 = 3.5
b = 2
last_node = maximum(edge)
Grp = "0.3Density"
