function population=Select(population,offspring,N)


joint = [population,offspring];
objs = [joint.objs]';
cons = [joint.cons]';

[~,index] = sortrows([cons,objs]);

joint = joint(index);


del = [];
for i=1:length(joint)-1
    if find(i==del)
        continue;
    end
   for j=i+1:length(joint)
       if isequal(joint(i).decs,joint(j).decs)
          del = [del j]; 
       end
   end
end
joint(del) = [];

population = joint(1:N);