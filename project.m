k = 100;
tolerance1 = 0.1;
tolerance2 = 0.01;
t = 0.1;

function v = rotate(v, theta)
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    v = R*v;
end

function T = extend(T)
   h = height(T.Nodes);
   index = randi(h);
   while T.Nodes.Steiner{index,1} == 1
       index = randi(height(T.Nodes));
   end
   T.Nodes.Steiner(index) = {1};

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

function Tcopy = update_graph(T, x)
    Tcopy = T;
    for i=1:height(Tcopy.Nodes)
        Tcopy.Nodes.x{i} = x(2*i-1);
        Tcopy.Nodes.y{i} = x(2*i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [c, ceq] = cons(x, T, b)
    % Inequality constraints. Edge lengths no greater than b
    c = [];
    Tnew = update_graph(T, x);
    for e = 1:height(Tnew.Edges)
      i = Tnew.Edges.EndNodes(e,1);
      j = Tnew.Edges.EndNodes(e,2);
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

function obj = f(x, T)
    obj = 0;
    for k=1:height(T.Edges)
        i = T.Edges.EndNodes(k,1);
        j = T.Edges.EndNodes(k,2);
        obj = obj + (x(2*i-1)-x(2*j-1))^2 + (x(2*i)-x(2*j))^2;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

function [p, gradp] = P(x, T, b, k)
    c = g(x, T, b);
    ceq = h(x, T);
    gradc = gradg(x, T);
    gradceq = gradh(x, T);

    weights = calculate_weights_at_x(T, x);

    p = sum(weights.^2) - (1/k)*sum(log(-c)) + (k/2)*sum(ceq.^2);
    gradp = gradf(T, x) - (gradc * (1./(k.*c)).') + (gradceq * k*ceq.');
end

function out = P_p(x, Tpenalty, b, k)
    [out,~] = P(x,Tpenalty,b,k);
end
function out = P_grad(x, Tpenalty, b, k)
    [~,out] = P(x,Tpenalty,b,k);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%


function T = generate_graph(n)
    A = [0 1; 1 0];
    NodeTable = table({0;0}, {0;1}, {0;0}, VariableNames={'x', 'y', 'Steiner'});
    T = graph(A, NodeTable);
    for i = 1:n-2
        T = extend(T);
    end
end

function xminEstimate = matlab_solver(T, b)

    % Create a conditions function with the treet T and edge length b built in
    specific_cons = @(x)cons(x, T, b);

    % Define the objective function
    
    objfun = @(x)f(x, T);

    x0 = [];
    for i=1:height(T.Nodes)
        x0(2*i-1) = T.Nodes.x{i};
        x0(2*i)   = T.Nodes.y{i};
    end

    xminEstimate = fmincon(objfun, x0, [], [], [], [], [], [], specific_cons);
end

function xminEstimate = penalty_method(T, b, k, tolerance1, tolerance2, t)
    x0 = [];
    for i=1:height(T.Nodes)
        x0(2*i-1) = T.Nodes.x{i};
        x0(2*i)   = T.Nodes.y{i};
    end
    x0 = x0.';
    p = @(x) P_p(x, T, b, k);
    gradp = @(x) P_grad(x, T, b, k);
    [xminEstimate, fminEstimate,k] = steepestDescentMethod(p, gradp, x0, tolerance1, tolerance2, t);
end

m = [];
p = [];
ns = [3 4 5 7 10 15];
%ns = [3 4];

for n=ns
    mtotal = 0;
    ptotal = 0;
    for i=1:3
        [n i]
        b = rand + 1.5;
        T = generate_graph(n);
        tic;
        mmins = matlab_solver(T, b);
        tmatlab = toc
        mtotal = mtotal + tmatlab;
        tic;
        pmins = penalty_method(T, b, k, tolerance1, tolerance2, t);
        tpenalty = toc
        ptotal = ptotal + tpenalty;
    end
    m = [m mtotal];
    p = [p ptotal];
end

plot(ns,m);
hold on
title('Total time taken to run 3 tests for each n value')
xlabel('n (Number of terminals)')
ylabel('Total time taken to run 3 tests (Seconds)')
plot(ns,p);
hold off
legend("MATLAB's built in NLP solver",'Our log-barrier penalty method implementation')
