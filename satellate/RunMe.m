
clc;
clear;
close all;
tic; 


global Global;
Global.num_satellite = 4;  
Global.num_object = 70;   
%% 读取数据
a = readmatrix('data/G.csv','Range','B2:B5');
Global.rank_satellite = a';
a = readmatrix('data/P.csv','Range','B2:B97');
Global.rank_object = a'; 
a = readmatrix('data/need.csv','Range','B2:B97');
Global.sat_need_time = a'; 

Global.visible_window = cell(Global.num_object, Global.num_satellite);
Global.num_visible_window = zeros(Global.num_object, Global.num_satellite);


for i = 1:Global.num_object
    datfile = ['data/sat' num2str(i) '.csv'];
    a = readmatrix(datfile, 'Range', 'B1:U4');
    for j = 1:Global.num_satellite
        index = a(j,:) ~= 0; % 创建一个索引向量 index，该向量包含了 a(j,:) 中所有非零元素的下标
        Global.visible_window{i, j} = a(j, index);
        Global.num_visible_window(i, j) = numel(Global.visible_window{i, j}) / 2;
    end
end

%% 算法参数
maxgen = 300;
popsize = 100;
population = Init(popsize);
tempt=0;
trace_obj = zeros(1, maxgen);
trace_con = zeros(1, maxgen);
priority_product = zeros(1, maxgen); 

% 参数初始化
N = 100;  
M = 2;  
omega = 0.5;  
c1 = 1.5; 
c2 = 1.5;  
lambda = eye(M); 


P = initializePopulation(N, M);
A = [];
HV_values = zeros(1, maxgen);

%% 进化开始
for i = 1:maxgen
    % 档案集更新
    A = archiveUpdate(A, P, N);
    % 克隆操作
    C = cloneOperation(A, N);
    % 最小度粒子群搜索
    P = PSO_based_search(P, A, lambda, omega, c1, c2);
    offspring = Mutate(population, i / maxgen);
    population = Select(population, offspring, popsize);
    bestobj = population(1).objs;
    trace_obj(i) = bestobj;
    trace_con(i) = population(1).cons;
    best_object_list = population(1).decs(1:Global.num_object);
    best_satellite_list = population(1).decs(Global.num_object + 1:end);

    total_rank = 0;
    p = 0;
    for j = 1:length(best_object_list)
        cur_object = best_object_list(j);
        cur_satellite = best_satellite_list(j); 
        p = p + Global.rank_satellite(cur_satellite) * Global.rank_object(cur_object);
        if p > total_rank
            total_rank = p;
        end
    end
    priority_product(i) = total_rank; 

    ref_point = [max([population.objs]) + 1, max([population.cons]) + 1]; 
    HV_values(i) = CalculateHV(population, ref_point);
    
    if ~mod(i, 50)
        cons = [population.cons];
        num = sum(cons == 0);
        avgcons = mean(cons);
        disp(['第' num2str(i) '代，满足约束个体数量：' num2str(num), '，最佳个体：' num2str(bestobj)])
    end
end

figure
plot(trace_obj)
title('最优目标值进化示意图')

figure
plot(HV_values)
title('Hypervolume (HV) 进化示意图')
xlabel('代数')
ylabel('HV值')



for i=1:maxgen
tempt=trace_obj(i);
trace_obj(i)=trace_obj(maxgen-i+1);
trace_obj(maxgen-i+1)=tempt;
end

bestsol = population(1);

% 保存每代优先级乘积数据
save('priority_product_data.mat', 'priority_product');
save('HV_values_data.mat','HV_values');

pareto_plot(trace_obj, priority_product,maxgen);
drawresult(bestsol);
toc

disp(['最终的HV值: ', num2str(HV_values(end))]);

disp(['运行时间: ', num2str(toc)]);



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
    % 使用自适应PBI方法选择pbesti
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
    % 选择与xi向量角最小的解作为gbesti
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
           
            parent1 = C(i).decs;
            parent2 = C(randi(length(C))).decs;
            crossoverPoint = randi(length(parent1));
            C(i).decs = [parent1(1:crossoverPoint), parent2(crossoverPoint+1:end)];
        end
        if rand < 0.4
          
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

        

