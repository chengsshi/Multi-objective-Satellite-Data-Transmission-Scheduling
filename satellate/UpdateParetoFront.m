function pareto_front = UpdateParetoFront(pareto_front, population)

    new_solutions = [pareto_front, population];

    [pareto_front, ~] = NondominatedSort(new_solutions);
end

function [front, ranks] = NondominatedSort(population)

    N = length(population);
    ranks = zeros(N, 1);
    front = [];
    
    for i = 1:N
        dominated = false;
        for j = 1:N
            if Dominates(population(j), population(i))
                dominated = true;
                break;
            end
        end
        if ~dominated
            front = [front, population(i)];
        end
    end
end

function is_dominated = Dominates(ind1, ind2)

    is_dominated = all(ind1.objs <= ind2.objs) && any(ind1.objs < ind2.objs);
end

function hv = CalculateHV(pareto_front, ref_point)
 
    hv = 0;
    for i = 1:length(pareto_front)
        hv = hv + prod(ref_point - pareto_front(i).objs);
    end
end
