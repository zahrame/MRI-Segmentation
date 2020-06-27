%This function converts a 2D matrix to 1D matrix
function y=conv2Dto1D(x)
[row,col]=size(x);
w(row*col)=0;
for i=1:row
    for j=1:col
        w((i-1)*col+j)=x(i,j);
    end
end
y=w';