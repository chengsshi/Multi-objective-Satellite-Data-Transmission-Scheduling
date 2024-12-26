function population = CalObj(population)
% population = Init(1);

global Global

N = length(population);

for i=1:N
    ind = population(i);
    object_list = ind.decs(1:Global.num_object);
    satellite_list = ind.decs(Global.num_object+1:end);
    
    satellite_next_release_time = zeros(1,Global.num_satellite);
    
    time_start_guance = zeros(1,Global.num_object);
    time_end_guance = zeros(1,Global.num_object);
    index_window_guance = zeros(1,Global.num_object); 
    cons = 0;
    
    for j = 1:Global.num_object
        cur_object = object_list(j); 
        cur_satellite = satellite_list(j); 
        
  
        flag = 0;
        for m = 1:Global.num_visible_window(cur_object,cur_satellite)
            time_start = Global.visible_window{cur_object,cur_satellite}(2*m-1); 
            time_end = Global.visible_window{cur_object,cur_satellite}(2*m); 
            if satellite_next_release_time(cur_satellite) > time_end
                continue;
            end
            time_begin = max(satellite_next_release_time(cur_satellite),time_start);
            
            if time_begin < time_end - Global.sat_need_time(cur_object) %当前观测窗口可以用来观测
                time_start_guance(cur_object) = time_begin;
                time_end_guance(cur_object) = time_start_guance(cur_object) + Global.sat_need_time(cur_object);
                satellite_next_release_time(cur_satellite) = time_end_guance(cur_object) + 60; %卫星的反应时间
                index_window_guance(cur_object) = m;
                flag=1;
                break;
            end
        end
        if flag == 0 %说明这个观测目标没有被安排到地面站
            cons = cons + Global.sat_need_time(cur_object);
        end
    end
    
    %% 计算目标函数
    T = max(time_end_guance);
    total_rank = 0;
    for j = 1:Global.num_object
        cur_object = object_list(j); %当前观测目标
        cur_satellite = satellite_list(j); %当前观测目标选择的地面站
        total_rank = total_rank + Global.rank_satellite(cur_satellite) * Global.rank_object(cur_object);
    end
    
    %% 传出值
    population(i).objs = T - 10*total_rank;
    population(i).cons = cons;
    population(i).time_start_guance = time_start_guance;
    population(i).time_end_guance = time_end_guance;
    population(i).satellite_list = satellite_list;
    population(i).index_window_guance = index_window_guance;
    population(i).satellite_next_release_time = satellite_next_release_time;
end