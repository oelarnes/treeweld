clear; close all; clc

load edges_sub.dat;
edges = edges_sub;

% n0 is the number of apparent edges
n0 = edges(1, 1);

sub_div = edges(1, 2); % subdivisions between apparent nodes.

% n is the number of welding edges
n = n0 * sub_div;

edges = edges(2:end,:);

K = size(edges, 1) / (n + 2 * n0 + 1); % n0 for degree data, n for expanded lamination

% m is the number of non-welding subdivisions. Should be 1.
m = 1;
% M is for normalization
M = 10**5;

% Z is a vector of vertex data.
Z = [2 * n * m];

% c is the second complex moment of the tree
c = [];

% the welding map is constructed out of the two conformal maps below. 
% slit_map welds an interval to angle alpha, alpha_moeb resets three given points
% to the positions where slip_map will weld them up.

function value = slit_map(z, alpha)
    value = exp(alpha * log(z + i * alpha) + (1 - alpha) * log(z + i * (alpha - 1)));
end

function value = alpha_moeb(z, x1, x2, x3, alpha)
	% alpha_moeb resets the points on the imaginary axis so that
	% the points x1, x2, x3 map to i(1-alpha), 0, -i*alpha. 
	
    value = -i * (alpha - 1) * alpha * (x1 - x3) * (z - x2) ./ ...
        (alpha * (x1 - x3) * (z - x2) + (x3 - x2) * (z - x1));
end

% Used in vertex improvement. W == 0 if and only if the v are Shabat nodes.
function W = error_fun(v, vd)
    n = size(v, 1) - 1;
    W = v - v.';
    W += (W == 0);
    W = log(W);
    W = sum(W, 2);
    W += log(vd) - log(4 * n);
    W -= pi * i * round(imag(W) / pi);
    W(n+1) = vd'*v;
end

% For Newton's method.
function DW = error_fun_derivative(v, vd)
    n = size(v, 1) - 1;
    DW = v - v.' + eye(n+1);
    DW = - ones(size(DW))./DW + eye(n+1);
    w = -1 * sum(DW, 2);
    DW += diag(w);
    DW(n+1,:) = vd';
end

for k = [1:K]

% unroll degree data for spacing
d0 = edges((k - 1)* (n + 2 * n0 + 1) + 1: (k - 1) * (n + 2 * n0 + 1) + n0, :)'(:);

% need these later: vertices
vert = edges((k - 1) * (n + 2*n0 + 1) + n0 + 1: ...
    (k - 1) * (n + 2 * n0 + 1) + 2*n0 + 1, :);

% lamination
L = edges((k - 1) * (n + 2 * n0 + 1) + 2 * n0 + 2: k * (n + 2 * n0 + 1), :);

% new degree data for welding
d = ones(2 * n, 1);

a = [1:sub_div]'/sub_div; % the base vector to be copied 2 * n0 times around the circle

% beta function has growth asymptotic to x^(d/2) at 0 and 1, so that
% after the weld, vertices will be linearly spaced. 

b = betainc(a, d0(2 * n0)/2, d0(1)/2) * (pi * i / n0);

z = b;

for j = [1: 2 * n0 - 1]

    b = betainc(a, d0(j)/2, d0(j+1)/2) * (pi * i / n0);
             
    z = [z; b + pi * i * j / n0]; 
end

figure(1)
z = exp(z);

% uncomment to see pre-weld vertex spacing.

%plot(z)
%hold on
%plot(z, 'rx')
%hold off

% move to the half plane, add a point for infinity

% this is awkward because -1 is always a starting vertex.
% Doesn't seem to be source of failures, however.

z = [z; M; -M; i * M; -i * M];
z = (z - 1) ./ (z + 1);
z = [z; 1];

% the welding loop
for j = [1:n-1]
    % i1,i2,i3, i4 are the indices of the endpoints of the next interval to be
    % welded
    i1 = L(j, 1);
    i2 = mod(L(j, 1) - 2, 2 * n) + 1;
    i3 = L(j ,2);
    i4 = mod(L(j, 2) - 2, 2 * n) + 1;
    % the values on the circle
	x1 = z(m * i1);
	x2 = .5 * z(m * i2) + .5 * z(m * i3);
	x3 = z(m * i4);
	% the angle to evenly space the new vertex
	alpha = d(i4) / (d(i4) + d(i1));
	% perform the weld
	z = alpha_moeb(z, x1, x2, x3, alpha);
	z = slit_map(z, alpha);
	% degree of new vertex
	d(i1) += d(i4);
	d(i4) = d(i1);
end

% back to the unit circle

z = (z + 1) ./ (z - 1);

i1 = L(n, 1);
i2 = mod(L(n, 1) - 2, 2 * n) + 1;
i3 = L(n ,2);
i4 = mod(L(n, 2) - 2, 2 * n) + 1;

x1 = .5 * z(m * i1) + .5 * z(m * i4);
x2 = .5 * z(m * i2) + .5 * z(m * i3);

% two remaining endpoints map to -1, 1

z = (x1 + x2 + i * (z - x1) * abs(x1 + x2) - x1 * z * conj(x1 + x2)) ./ ...
    (x1 + x2 - i * (z - x1) * abs(x1 + x2) - x1 * z * conj(x1 + x2));

% weld the circle

z = (z + ones(size(z)) ./ z);

x = z(2 * n * m + 5);
z = z(1: 2 * n * m + 4);

% normalize by resetting infinity

z = ones(size(z)) ./ (z - x);

% then normalize derivative and shift at infinity

% normalization with M, -M, etc.
%b = sum(z) / (2 * n * m);
%a = sum(z ./ (exp([1:2 * n * m]' /(n * m) * pi * i))) / (2 * n * m);

b = 1/4 * (z(2 * n * m + 1) + z(2 * n * m + 2) + z(2 * n * m + 3) + z(2 * n * m + 4));
a = 1/(4 * M) * (z(2 * n * m + 1) - z(2 * n * m + 2) - i * z(2 * n * m + 3) + i * z(2 * n * m + 4));
z = z(1 : 2 * n * m);
z = (z - b) / a;
Z = [Z; z];

% calculate the second moment of the new data set
c = [c; sum(z.**2) / (4 * n * m)];

if abs(c) > 1
    fprintf('error: bad weld\n')
%else
%    v = [z(vert(:,1)*sub_div)+z(vert(:,2)*sub_div)]/2;
%    vd = [d0(vert(:,1))];
%    iterations = 8;
%    for i = [1:iterations]
%        % Newton's method
%        W = error_fun(v, vd);
%        DW = error_fun_derivative(v, vd);
%        dv = DW \ W;
%        v -= dv;
%    end
    % to calculate Shabat polynomials, we have the vector of updated vertex
    % data with the degree of each vertex. It would be better to have a list of 
    % the 2n images after improvement, but this would require some additional work
endif

end
c(abs(c) > 1) = 0;

if K > 1
    % plot the c data. Each data point has argument twice the argument of the "long axis"
    % of the tree, with magnitude the square of 
    figure(1);
    plot(c, '.')
    hold on;
    plot([exp(pi * i * [1:100] / 50)])
    hold off;
    figure(2);
    hist(sqrt(abs(c)), 20)
end
figure(3);
hold off
rand_ind = floor(rand() * K);
z = Z(rand_ind * 2 * n * m + 2: (rand_ind + 1) * 2 * n * m + 1, :);
plot(z)
hold on
% plot(v, 'or')
% uncomment to see vertex spacing
% hold on
% plot(z, 'rx')
% hold off
axis([-2.5 2.5 -2.5 2.5])
axis('equal')