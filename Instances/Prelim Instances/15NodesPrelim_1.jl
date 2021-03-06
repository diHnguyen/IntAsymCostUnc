edge = [1 3; 1 7; 2 6; 3 2; 3 5; 3 14; 4 12; 5 3; 5 9; 5 15; 6 3; 6 5; 6 13; 7 3; 7 4; 7 12; 7 13; 8 10; 8 15; 9 2; 9 4; 9 8; 9 12; 10 4; 10 7; 10 12; 11 2; 11 5; 11 7; 11 13; 12 3; 12 4; 12 11; 12 14; 12 15; 13 2; 13 5; 13 11; 14 7; 14 9; 14 15]
cL_orig = [19.0, 58.0, 45.0, 5.0, 25.0, 100.0, 76.0, 11.0, 36.0, 103.0, 26.0, 6.0, 72.0, 41.0, 30.0, 54.0, 57.0, 15.0, 74.0, 73.0, 50.0, 8.0, 31.0, 58.0, 24.0, 21.0, 95.0, 57.0, 38.0, 22.0, 87.0, 76.0, 12.0, 9.0, 35.0, 106.0, 80.0, 18.0, 67.0, 35.0, 11.0]
cU_orig = [19.0, 58.0, 45.0, 5.0, 25.0, 118.0, 76.0, 31.0, 36.0, 103.0, 26.0, 6.0, 72.0, 41.0, 30.0, 54.0, 57.0, 29.0, 74.0, 73.0, 50.0, 8.0, 31.0, 58.0, 36.0, 21.0, 95.0, 73.0, 38.0, 22.0, 87.0, 76.0, 12.0, 25.0, 35.0, 106.0, 80.0, 18.0, 81.0, 55.0, 11.0]
d = [17.0, 19.0, 13.0, 8.0, 12.0, 4.0, 2.0, 2.0, 12.0, 5.0, 17.0, 1.0, 15.0, 17.0, 1.0, 15.0, 9.0, 14.0, 14.0, 3.0, 7.0, 3.0, 4.0, 13.0, 4.0, 15.0, 14.0, 2.0, 12.0, 10.0, 16.0, 3.0, 15.0, 20.0, 17.0, 5.0, 10.0, 16.0, 17.0, 10.0, 13.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
SP_init = sum(yy[i]*c_orig[i] for i = 1:Len)
p = [1.0]
g = [SP_init]
h = [0.0]
origin = 1
destination = 15
last_node = maximum(edge)
all_nodes = collect(1:last_node)
M_orig = zeros(Len)
for i = 1:Len
    M_orig[i] = cU_orig[i] - cL_orig[i]
end
delta1 = 1.0
delta2 = 3.0
last_node = maximum(edge)
b=2
Grp = "15NodesPrelim"