function HV = Hypervolume(front, ref_point)

    
    if isempty(front)
        HV = 0;
        return;
    end

 
    if size(front, 2) ~= 2
        error('front 必须是2维矩阵');
    end
    
    % 排序
    front = sortrows(front, 1);
    
    % 初始化HV
    HV = 0;
    previous = ref_point(1);
    
    for i = 1:size(front, 1)
        HV = HV + (previous - front(i, 1)) * (ref_point(2) - front(i, 2));
        previous = front(i, 1);
    end
end
