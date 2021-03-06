edge = [1 8; 1 10; 1 15; 1 16; 2 5; 2 7; 2 10; 2 14; 2 16; 3 4; 3 13; 3 15; 3 18; 4 3; 4 10; 4 14; 4 16; 5 2; 5 6; 5 7; 5 8; 5 13; 5 14; 6 3; 6 10; 6 12; 7 8; 7 9; 7 13; 7 20; 8 6; 8 11; 8 19; 8 20; 9 10; 9 16; 10 5; 10 7; 10 9; 10 11; 10 16; 11 18; 11 20; 12 3; 12 6; 12 10; 12 20; 13 2; 13 3; 13 5; 13 7; 13 16; 13 20; 14 7; 14 8; 14 13; 14 20; 15 6; 15 11; 15 12; 15 13; 16 2; 16 6; 16 7; 16 11; 17 14; 17 19; 17 20; 18 4; 18 11; 18 20; 19 11; 19 15; 19 17]
cL_orig = [67.0, 87.0, 140.0, 149.0, 24.0, 40.5, 80.0, 123.0, 145.0, 3.0, 97.0, 115.0, 149.5, 13.0, 52.0, 91.5, 115.0, 21.5, 15.0, 19.0, 30.0, 80.0, 90.5, 25.0, 41.0, 50.0, 15.0, 11.0, 60.0, 134.0, 17.0, 27.0, 109.0, 117.0, 15.0, 73.0, 51.0, 25.0, 6.0, 15.0, 54.0, 74.0, 94.0, 94.0, 62.0, 21.0, 75.0, 108.0, 99.0, 73.5, 62.0, 29.5, 74.0, 73.0, 60.0, 7.0, 59.0, 93.0, 40.0, 32.0, 17.0, 135.0, 101.0, 87.0, 53.0, 33.0, 23.0, 29.0, 143.0, 75.0, 11.0, 79.5, 42.0, 17.0]
cU_orig = [67.0, 97.0, 140.0, 157.0, 36.0, 55.5, 80.0, 123.0, 145.0, 11.0, 97.0, 115.0, 156.5, 13.0, 60.0, 98.5, 115.0, 28.5, 15.0, 19.0, 30.0, 90.0, 99.5, 25.0, 41.0, 62.0, 15.0, 25.0, 60.0, 134.0, 33.0, 27.0, 109.0, 123.0, 15.0, 73.0, 51.0, 25.0, 6.0, 15.0, 62.0, 74.0, 94.0, 94.0, 66.0, 21.0, 75.0, 108.0, 99.0, 84.5, 62.0, 36.5, 74.0, 73.0, 60.0, 7.0, 59.0, 93.0, 40.0, 32.0, 17.0, 135.0, 101.0, 87.0, 53.0, 33.0, 23.0, 37.0, 143.0, 75.0, 19.0, 90.5, 42.0, 17.0]
d = [7.0, 4.0, 9.0, 5.0, 9.0, 13.0, 17.0, 20.0, 8.0, 7.0, 19.0, 9.0, 2.0, 9.0, 19.0, 5.0, 5.0, 12.0, 1.0, 20.0, 13.0, 4.0, 18.0, 8.0, 10.0, 6.0, 10.0, 5.0, 17.0, 7.0, 17.0, 17.0, 7.0, 8.0, 8.0, 3.0, 4.0, 5.0, 13.0, 20.0, 16.0, 14.0, 15.0, 10.0, 7.0, 13.0, 15.0, 15.0, 13.0, 17.0, 9.0, 14.0, 7.0, 13.0, 2.0, 2.0, 9.0, 14.0, 4.0, 5.0, 19.0, 9.0, 5.0, 2.0, 20.0, 13.0, 14.0, 11.0, 7.0, 9.0, 12.0, 8.0, 15.0, 17.0]
Len = length(cL_orig)
c_orig = 0.5*(cL_orig+cU_orig)
yy = zeros(Len)
yy = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
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
Grp = "+10%Uncertainties"
