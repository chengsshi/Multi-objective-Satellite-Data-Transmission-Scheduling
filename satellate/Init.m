function population = Init(N)

global Global
empty.decs = [];
empty.objs = [];
empty.cons = [];

population = repmat(empty,1,N);
for i=1:N
    population(i).decs = [randperm(Global.num_object,Global.num_object) ...
        randi([1,Global.num_satellite],1,Global.num_object)];
end
population = CalObj(population);