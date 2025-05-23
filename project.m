% We build a random topology by starting with a base tree, and
% extend it by one terminal at a time.
% We do this by taking a terminal at random, converting it to a Steiner
% point, and adding two new terminals with edges to the point we picked
% before. This has the net effect of adding one terminal and one
% Steiner point to the graph.

% n is the number of terminals
n = 3;

% Start with the most simple base tree. Two terminals joined by an edge.
% A is our adjacency matrix
A = [
    0 1
    1 0
];

% For each node, we store its x and y positions, and whether it is a
% Steiner point
NodeTable = table({0;0}, {0;1}, {0;0}, VariableNames={'x', 'y', 'Steiner'});
T = graph(A, NodeTable);

% We extend the graph by adding n-2 terminals (base graph has 2 terminals)
for i = 1:n-2
    T = extend(T);
end

% This function displays the graph with the points shown in their correct
% positions
function display(T)
    x = [];
    y = [];

    for i = 1:height(T.Edges.EndNodes)
        p1 = T.Edges.EndNodes(i,1);
        p2 = T.Edges.EndNodes(i,2);
        x1 = T.Nodes.x{p1,1};
        x2 = T.Nodes.x{p2,1};
        y1 = T.Nodes.y{p1,1};
        y2 = T.Nodes.y{p2,1};
    
        x = [x [x1; x2]];
        y = [y [y1; y2]];
    end

    plot(x, y);
end

% Rotates a vector by a given angle in degrees
% This function is used to help create a ckean graph
function v = rotate(v, theta)
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    v = R*v;
end

% This function extends our tree by adding a terminal
function T = extend(T)
   % First, randomly pick the terminal and change it to a Steiner point
   h = height(T.Nodes);
   index = randi(h);
   while T.Nodes.Steiner{index,1} == 1
       index = randi(height(T.Nodes));
   end
   T.Nodes.Steiner(index) = {1};

   % Add the two new terminals
   % We first decide the positions of the new points
   % One way to do this is to put the two new points one unit away from the
   % Steiner, at about 30 degrees offset from the direction to this point's
   % existing neighbour
   neighbor = neighbors(T, index);
   s_1 = [T.Nodes.x{index,1}; T.Nodes.y{index,1}];
   s_2 = [T.Nodes.x{neighbor,1}; T.Nodes.y{neighbor,1}];
   v = s_1 - s_2;
   v_1 = rotate(v, 30);
   v_2 = rotate(v, -30);
   x_1 = v_1(1) / norm(v_1) + s_1(1) + (rand() - 0.5) * 0.3;
   y_1 = v_1(2) / norm(v_1) + s_1(2) + (rand() - 0.5) * 0.3;
   x_2 = v_2(1) / norm(v_2) + s_1(1) + (rand() - 0.5) * 0.3;
   y_2 = v_2(2) / norm(v_2) + s_1(2) + (rand() - 0.5) * 0.3;

   % Must properly calculate locations
   NodeProperties_1 = table({x_1}, {y_1}, {0}, VariableNames={'x', 'y', 'Steiner'});
   T = addnode(T, NodeProperties_1);
   NodeProperties_2 = table({x_2}, {y_2}, {0}, VariableNames={'x', 'y', 'Steiner'});
   T = addnode(T, NodeProperties_2);

   % Finally, we add the edges between our Steiner and new terminals
   T = addedge(T, [index index], [h+1 h+2], [0 0]);
end

% The edge weights in our graph are equal to the length of the edge
% This is useful for later analysis
function weights = calculate_weights(T)
weights = [];
    for i = 1:height(T.Edges)
        a = T.Edges.EndNodes(i,1);
        b = T.Edges.EndNodes(i,2);
        p1 = [T.Nodes.x{a, 1}; T.Nodes.y{a, 1}];
        p2 = [T.Nodes.x{b, 1}; T.Nodes.y{b, 1}];
        v = p2 - p1;
        dist = norm(v);
        weights(i) = dist;
    end
end

function weights = calculate_weights_at_x(T, x)
    Tcopy = T;
    for i=1:height(Tcopy.Nodes)
        Tcopy.Nodes.x{i} = x(2*i-1);
        Tcopy.Nodes.y{i} = x(2*i);
    end
    weights = [];
    for i = 1:height(Tcopy.Edges)
        a = Tcopy.Edges.EndNodes(i,1);
        b = Tcopy.Edges.EndNodes(i,2);
        p1 = [Tcopy.Nodes.x{a, 1}; Tcopy.Nodes.y{a, 1}];
        p2 = [Tcopy.Nodes.x{b, 1}; Tcopy.Nodes.y{b, 1}];
        v = p2 - p1;
        dist = norm(v);
        weights(i) = dist;
    end
end

function gradf = gradf(T, x)
    update_graph(T, x);
    gradf = zeros([height(x) 1]);
    for i = 1:height(T.Edges)
        a = T.Edges.EndNodes(i,1);
        b = T.Edges.EndNodes(i,2);
        grad = zeros([height(x) 1]);
        grad(a*2-1) = 2*x(a*2-1) - 2*x(b*2-1);
        grad(b*2-1) = 2*x(b*2-1) - 2*x(a*2-1);
        grad(a*2) = 2*x(a*2) - 2*x(b*2);
        grad(b*2) = 2*x(b*2) - 2*x(a*2);
        gradf = gradf + grad;
    end
end

%weights = calculate_weights(T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now we have a tree T. Let's solve the objective function

b = 2;

% We'll construct our problem so all points are considered decision
% variables, but all fixed points will have their coordinates fixed by
% constraints.

% x coords are at x(2n-1)
% y coords are at x(2n)
% x = [x_1 y_2 x_2 y_2 ... x_n y_n]

function [c, ceq] = cons(x, T, b)
    % Inequality constraints. Edge lengths no greater than b
    c = [];
    for e = 1:height(T.Edges)
      i = T.Edges.EndNodes(e,1);
      j = T.Edges.EndNodes(e,2);
      c(e) = (x(2*i-1)-x(2*j-1))^2 + (x(2*i)-x(2*j))^2 - b^2;
    end

    % Equality constraints. Fix the locations of the fixed points
    ceq = [];
    for i=1:height(T.Nodes)
      if T.Nodes.Steiner{i}==0
        ceq(i*2-1) = x(2*i-1) - T.Nodes.x{i};
        ceq(i*2)   = x(2*i) - T.Nodes.y{i};
      end
    end
end

% Create a conditions function with the treet T and edge length b built in
specific_cons = @(x)cons(x, T, b);

% Define the objective function
function obj = f(x, T)
    obj = 0;
    for k=1:height(T.Edges)
        i = T.Edges.EndNodes(k,1);
        j = T.Edges.EndNodes(k,2);
        obj = obj + (x(2*i-1)-x(2*j-1))^2 + (x(2*i)-x(2*j))^2;
    end
end

% Create an objective function with T built in
objfun = @(x)f(x, T);

% We set our starting points to what the current configuration is
x0 = [];
for i=1:height(T.Nodes)
    x0(2*i-1) = T.Nodes.x{i};
    x0(2*i)   = T.Nodes.y{i};
end

% Now solve the objective function
x = fmincon(objfun, x0, [], [], [], [], [], [], specific_cons);

Tmatlab = T;

% Using the optimal coordinates we just found, update the graph.
for i=1:height(Tmatlab.Nodes)
    Tmatlab.Nodes.x{i} = x(2*i-1);
    Tmatlab.Nodes.y{i} = x(2*i);
end

%display(Tmatlab);
matlab_weights = calculate_weights(Tmatlab)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Tcopy = update_graph(T, x)
    Tcopy = T;
    for i=1:height(Tcopy.Nodes)
        Tcopy.Nodes.x{i} = x(2*i-1);
        Tcopy.Nodes.y{i} = x(2*i);
    end
end

% Inequality constraints
function g = g(x, T, b)
    T = update_graph(T, x);
    g = [];
    for e = 1:height(T.Edges)
      i = T.Edges.EndNodes(e,1);
      j = T.Edges.EndNodes(e,2);
      g(e) = (x(2*i-1)-x(2*j-1))^2 + (x(2*i)-x(2*j))^2 - b^2;
    end
end

% Gradient of inequality constraints
function gradg = gradg(x, T)
    T = update_graph(T, x);
    gradg = [];
    for e = 1:height(T.Edges)
      i = T.Edges.EndNodes(e,1);
      j = T.Edges.EndNodes(e,2);
      gradgi = zeros([height(x) 1]);
      gradgi(2*i-1) = 2*x(2*i-1) - 2*x(2*j-1);
      gradgi(2*j-1) = 2*x(2*j-1) - 2*x(2*i-1);
      gradgi(2*i)   = 2*x(2*i)   - 2*x(2*j);
      gradgi(2*j)   = 2*x(2*j)   - 2*x(2*i);
      gradg = [gradg gradgi];
    end
end

% Equality constraints
function h = h(x, T)
    %T = update_graph(T, x);
    h = [];
    for i=1:height(T.Nodes)
      if T.Nodes.Steiner{i}==0
        h = [h (x(2*i-1) - T.Nodes.x{i})];
        h = [h (x(2*i)   - T.Nodes.y{i})];
      end
    end
end

% Gradient of equality constraints
function gradh = gradh(x, T)
    %T = update_graph(T, x);
    gradh = [];
    for i=1:height(T.Nodes)
        if T.Nodes.Steiner{i}==0
            x = zeros([height(x) 1]);
            x(2*i-1) = 1;
            y = zeros([height(x) 1]);
            y(2*i) = 1;
            gradh = [gradh x y];
        end
    end
end


% Now that we've solved the problem using MATLAB's implementation
% We will use the steepest descent method after applying a log barrier
% We will also use BFGS to solve the single variable step size t

Tpenalty = T;

% To generate the penalty function, we need the constraint functions
% We can use the function from above

function [p, gradp] = P(x, T, b, k)
    c = g(x, T, b);
    ceq = h(x, T);
    gradc = gradg(x, T);
    gradceq = gradh(x, T);

    weights = calculate_weights_at_x(T, x);

    p = sum(weights.^2) - (1/k)*sum(log(-c)) + (k/2)*sum(ceq.^2);
    gradp = gradf(T, x) - (gradc * (1./(k.*c)).') + (gradceq * k*ceq.');
end

k = 10;
p = @(x) P_p(x, Tpenalty, b, k);
gradp = @(x) P_grad(x, Tpenalty, b, k);

function out = P_p(x, Tpenalty, b, k)
    [out,~] = P(x,Tpenalty,b,k);
end
function out = P_grad(x, Tpenalty, b, k)
    [~,out] = P(x,Tpenalty,b,k);
end


tolerance1 = 0.1;
tolerance2 = 0.1;

x0 = x0.';

display(Tpenalty);

for i=0:4
    k = 2^i;
    p = @(x) P_p(x, Tpenalty, b, k);
    gradp = @(x) P_grad(x, Tpenalty, b, k);
    [xminEstimate, fminEstimate,k] = steepestDescentMethod(p, gradp, x0, tolerance1, tolerance2, 1)
    %weights = calculate_weights_at_x(T, xminEstimate)
    %c = g(xminEstimate, T, b)
    %ceq = h(xminEstimate, T)
    %sum(weights.^2)
    %(1/k)*sum(log(-c))
    %(k/2)*sum(ceq.^2)
end

xminEstimate
display(update_graph(Tpenalty, xminEstimate));
