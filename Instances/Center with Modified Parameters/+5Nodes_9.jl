edge = [1 2; 1 5; 1 12; 1 13; 2 9; 2 13; 2 14; 2 15; 3 2; 3 5; 3 11; 3 15; 4 16; 5 10; 5 11; 5 12; 5 18; 6 5; 6 20; 7 6; 7 18; 8 12; 8 15; 8 18; 8 19; 9 2; 9 5; 9 8; 9 17; 9 19; 10 3; 10 8; 10 9; 10 17; 10 18; 11 5; 11 12; 12 2; 12 5; 12 8; 12 10; 13 4; 13 9; 13 12; 13 20; 14 2; 14 4; 14 10; 14 16; 14 20; 15 2; 15 7; 15 16; 16 4; 16 5; 16 7; 16 13; 16 20; 17 2; 17 5; 17 8; 18 2; 18 5; 18 9; 18 11; 18 19; 19 3; 19 5; 19 7; 19 11; 19 15; 19 16; 19 18; 1 23; 3 23; 3 25; 4 23; 5 25; 6 22; 10 22; 10 25; 13 24; 13 25; 14 23; 17 22; 17 24; 17 25; 18 21; 18 22; 19 22; 19 24; 21 6; 21 9; 21 12; 22 4; 22 11; 23 3; 24 4; 24 12; 24 13; 24 16; 25 10; 25 11; 25 12; 25 15; 25 16; 25 21; 25 23]
cL_orig = [6.0, 44.0, 105.0, 124.0, 64.0, 108.0, 120.0, 132.0, 13.0, 17.0, 85.0, 118.0, 119.0, 42.0, 55.0, 69.0, 130.0, 12.0, 140.0, 9.0, 102.0, 41.0, 73.0, 101.0, 108.0, 73.0, 33.0, 12.0, 77.0, 105.0, 69.0, 18.0, 9.0, 66.0, 84.0, 57.0, 3.0, 100.0, 71.0, 38.0, 19.0, 89.0, 37.0, 5.0, 67.0, 124.0, 91.0, 39.0, 25.0, 46.0, 135.0, 74.0, 14.0, 120.0, 100.0, 88.0, 35.0, 31.0, 150.0, 116.0, 87.0, 163.0, 133.0, 95.0, 67.0, 13.0, 161.0, 139.0, 123.0, 85.0, 37.0, 32.0, 12.0, 224.0, 201.0, 224.0, 193.0, 205.0, 159.0, 121.0, 153.0, 112.0, 122.0, 86.0, 48.0, 66.0, 83.0, 34.0, 35.0, 25.0, 46.0, 152.0, 125.0, 92.0, 177.0, 107.0, 205.0, 201.0, 117.0, 111.0, 82.0, 155.0, 145.0, 126.0, 98.0, 86.0, 35.0, 23.0]
cU_orig = [20.0, 44.0, 105.0, 124.0, 76.0, 122.0, 120.0, 132.0, 13.0, 33.0, 85.0, 118.0, 119.0, 60.0, 67.0, 69.0, 130.0, 12.0, 140.0, 9.0, 118.0, 41.0, 73.0, 101.0, 108.0, 73.0, 45.0, 12.0, 77.0, 105.0, 79.0, 18.0, 9.0, 66.0, 84.0, 57.0, 23.0, 100.0, 71.0, 38.0, 19.0, 89.0, 37.0, 5.0, 67.0, 124.0, 105.0, 39.0, 25.0, 64.0, 135.0, 86.0, 14.0, 120.0, 116.0, 88.0, 35.0, 41.0, 150.0, 116.0, 87.0, 163.0, 133.0, 95.0, 67.0, 13.0, 161.0, 139.0, 123.0, 85.0, 37.0, 32.0, 12.0, 224.0, 201.0, 224.0, 193.0, 205.0, 159.0, 121.0, 153.0, 112.0, 122.0, 86.0, 48.0, 66.0, 83.0, 34.0, 35.0, 25.0, 46.0, 152.0, 125.0, 92.0, 177.0, 107.0, 205.0, 201.0, 117.0, 111.0, 82.0, 155.0, 145.0, 126.0, 98.0, 86.0, 35.0, 23.0]
d = [10.0, 8.0, 17.0, 7.0, 19.0, 17.0, 11.0, 12.0, 18.0, 19.0, 11.0, 3.0, 20.0, 5.0, 13.0, 2.0, 9.0, 2.0, 3.0, 13.0, 14.0, 20.0, 15.0, 3.0, 20.0, 17.0, 6.0, 16.0, 6.0, 10.0, 10.0, 12.0, 18.0, 19.0, 18.0, 14.0, 4.0, 9.0, 1.0, 13.0, 15.0, 18.0, 15.0, 5.0, 20.0, 20.0, 4.0, 3.0, 16.0, 16.0, 4.0, 18.0, 6.0, 8.0, 10.0, 1.0, 9.0, 18.0, 20.0, 12.0, 8.0, 4.0, 2.0, 15.0, 16.0, 16.0, 6.0, 6.0, 7.0, 8.0, 13.0, 1.0, 12.0, 12.0, 8.0, 2.0, 18.0, 1.0, 9.0, 16.0, 20.0, 8.0, 11.0, 5.0, 8.0, 19.0, 10.0, 8.0, 1.0, 10.0, 6.0, 7.0, 15.0, 11.0, 9.0, 12.0, 14.0, 8.0, 7.0, 15.0, 16.0, 4.0, 16.0, 1.0, 10.0, 6.0, 6.0, 1.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
Grp = "+5Nodes"
