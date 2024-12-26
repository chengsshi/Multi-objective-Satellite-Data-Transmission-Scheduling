
N = 100;  
M = 2;  
omega = 0.5;  
c1 = 1.5; 
c2 = 1.5; 
lambda = eye(M); 

P = initializePopulation(N, M);


A = [];

maxGenerations = 100;
generation = 0;
while generation < maxGenerations
   
    A = archiveUpdate(A, P, N);
    
  
    C = cloneOperation(A, N);
    
   
    P = PSO_based_search(P, A, lambda, omega, c1, c2);
    
    generation = generation + 1;
end


display(A);

function P = initializePopulation(N, M)

    P = struct('decs', [], 'velocity', [], 'objs', []);
    for i = 1:N
        P(i).decs = rand(1, M); 
        P(i).velocity = zeros(1, M);  
        P(i).objs = calculateObjective(P(i).decs); 
    end
end

function objs = calculateObjective(decs)

    objs = [sum(decs.^2), sum((decs-2).^2)];
end

function P = PSO_based_search(P, A, lambda, omega, c1, c2)
    for i = 1:length(P)
        xi = P(i);

        pbesti = selectPBest(A, xi, lambda);
     
        gbesti = selectGBest(A, xi);
        
  
        theta = vectorAngle(xi.objs, lambda(:,1));  
        Fir2 = (pi/2 - theta) / (pi/2 + fitness(pbesti) / fitness(gbesti));
        xi.velocity = omega * xi.velocity + c1 * (pbesti.decs - xi.decs) + c2 * (fitness(gbesti) - fitness(xi)) * Fir2 * (gbesti.decs - xi.decs);
        
 
        xi.decs = xi.decs + xi.velocity;
        
    
        xi.decs = max(min(xi.decs, 1), 0);
    
        xi.objs = calculateObjective(xi.decs);
        
        P(i) = xi;
    end
end

function pbesti = selectPBest(A, xi, lambda)
 
    distances = zeros(1, length(A));
    for i = 1:length(A)
        distances(i) = pbi(A(i), xi, lambda(:,1));  
    end
    [~, idx] = min(distances);
    pbesti = A(idx);
end

function distance = pbi(sol, xi, lambda)
  
    theta = 1;  
    d1 = dot(sol.objs, lambda) / norm(lambda);
    d2 = norm(sol.objs - d1 * lambda / norm(lambda));
    distance = d1 + theta * d2;
end

function gbesti = selectGBest(A, xi)
 
    angles = zeros(1, length(A));
    for i = 1:length(A)
        angles(i) = vectorAngle(A(i).objs, xi.objs);
    end
    [~, idx] = min(angles);
    gbesti = A(idx);
end

function f = fitness(solution)
  
    f = sqrt(sum(solution.objs .^ 2));
end

function A = maxMinAngleSelection(Cnd, N)

    Cnd = removeDRS(Cnd);
    

    angles = calculateAngles(Cnd);
 
    A = addExtremeSolutions(Cnd);
   
    while length(A) < N
        [~, idx] = min(angles);
        A = [A, Cnd(idx)];
        Cnd(idx) = [];
        angles(idx) = [];
    end
end

function Cnd = removeDRS(Cnd)
 
    DRS_indices = [];
    for i = 1:length(Cnd)
        if isDRS(Cnd(i))
            DRS_indices = [DRS_indices, i];
        end
    end
    Cnd(DRS_indices) = [];
end

function is_DRS = isDRS(solution)
  
    is_DRS = any(solution.objs == 1);
end

function angles = calculateAngles(Cnd)
 
    angles = zeros(1, length(Cnd));
    for i = 1:length(Cnd)
        min_angle = inf;
        for j = 1:length(Cnd)
            if i ~= j
                angle = vectorAngle(Cnd(i).objs, Cnd(j).objs);
                if angle < min_angle
                    min_angle = angle;
                end
            end
        end
        angles(i) = min_angle;
    end
end

function angle = vectorAngle(vec1, vec2)
  
    cos_theta = dot(vec1, vec2) / (norm(vec1) * norm(vec2));
    angle = acos(cos_theta);
end

function A = addExtremeSolutions(Cnd)

    extreme_solutions = findExtremeSolutions(Cnd);
    A = extreme_solutions;
    Cnd = setdiff(Cnd, extreme_solutions);
end

function extreme_solutions = findExtremeSolutions(Cnd)
 
    extreme_solutions = [];
    ref_vectors = eye(size(Cnd(1).objs, 1));
    for i = 1:size(ref_vectors, 2)
        min_angle = inf;
        extreme_solution = [];
        for j = 1:length(Cnd)
            angle = vectorAngle(Cnd(j).objs, ref_vectors(:, i)');
            if angle < min_angle
                min_angle = angle;
                extreme_solution = Cnd(j);
            end
        end
        extreme_solutions = [extreme_solutions, extreme_solution];
    end
end

function A = archiveUpdate(A, P, N)
    C = [A, P];  
    C = normalizePopulation(C); 
    
    fitness_values = calculateFitness(C);  
    
  
    Cnd = selectNondominated(C);
    
    if length(Cnd) > N
        A = maxMinAngleSelection(Cnd, N);
    else
        A = Cnd;
    end
end

function C = normalizePopulation(C)
  
    objs = reshape([C.objs], [], length(C))';
    minObjs = min(objs, [], 1);
    maxObjs = max(objs, [], 1);
    for i = 1:length(C)
        C(i).objs = (C(i).objs - minObjs) ./ (maxObjs - minObjs);
    end
end

function fitness_values = calculateFitness(C)
    
    fitness_values = zeros(1, length(C));
    for i = 1:length(C)
        fitness_values(i) = sqrt(sum(C(i).objs .^ 2));
    end
end

function Cnd = selectNondominated(C)
   
    Cnd = [];
    for i = 1:length(C)
        dominated = false;
        for j = 1:length(C)
            if i ~= j && dominates(C(j), C(i))
                dominated = true;
                break;
            end
        end
        if ~dominated
            Cnd = [Cnd, C(i)];
        end
    end
end

function is_dominated = dominates(sol1, sol2)
   
    is_dominated = all(sol1.objs <= sol2.objs) && any(sol1.objs < sol2.objs);
end

function C = cloneOperation(A, N)
    B = A;  
    C = []; 
    
    SDE_distances = calculateSDEDistances(B); 
    
    SDE_sum = sum(SDE_distances);
    NC = length(B);
    ni = ceil(N * SDE_distances / SDE_sum);
    
    
    for i = 1:NC
        for j = 1:ni(i)
            C = [C, B(i)];
        end
    end
    
   
    C = evolvePopulation(C);
end

function SDE_distances = calculateSDEDistances(B)
   
    NC = length(B);
    SDE_distances = zeros(1, NC);
    for i = 1:NC
  
        SDE_distances(i) = calculateSDE(B, B(i));
    end
end

function distance = calculateSDE(B, solution)
  
    distance = 0;
    for i = 1:length(B)
        if ~isequal(B(i), solution)
            distance = distance + norm(solution.decs - B(i).decs);
        end
    end
end

function C = evolvePopulation(C)
  
    for i = 1:length(C)
        if rand < 0.8
            % 交叉操作
            parent1 = C(i).decs;
            parent2 = C(randi(length(C))).decs;
            crossoverPoint = randi(length(parent1));
            C(i).decs = [parent1(1:crossoverPoint), parent2(crossoverPoint+1:end)];
        end
        if rand < 0.4
            % 变异操作
            mutationPoint = randi(length(C(i).decs));
            C(i).decs(mutationPoint) = rand;
        end
    end
    C = CalObj(C);
end

function C = CalObj(C)
    
    for i = 1:length(C)
        C(i).objs = calculateObjective(C(i).decs);
    end
end

        
