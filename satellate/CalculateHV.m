function HV = CalculateHV(population, ref_point)
    objs = [population.objs];
    
    feasible_objs = objs([population.cons] == 0);
    
    if isempty(feasible_objs)
        HV = 0;
    else
      
        feasible_objs = reshape(feasible_objs, [], 1);
        feasible_objs = [feasible_objs, zeros(size(feasible_objs))];
        
        HV = Hypervolume(feasible_objs, ref_point);
    end
end
