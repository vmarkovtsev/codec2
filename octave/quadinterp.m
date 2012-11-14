% quadinterp.m
% David Rowe 10/3/11

x_1 = -10; y_1 = 1;
x0 = 0; y0 = 3;
x1 = 2; y1 = 2;

plot([ x_1 x0 x1 ], [y_1 y0 y1], "r+");

c = y0
a = (y_1*x1 - y1*x_1 + c*x_1 - c*x1)/(x_1*x_1*x1 - x1*x1*x_1)
b = (y1 -a*x1*x1 - c)/x1
polyfit([x_1 x0 x1],[y_1 y0 y1],2)

hold on;
for x = -11:0.25:3;
    y = a*x*x + b*x + c;
    plot(x,y,'o');
end
hold off;

