function[idx_time1,idx_time2]=find_idx_time_4ts(time_1,time_2,v_time)


idx_time1=find(v_time==time_1);
if isempty(idx_time1)
    v_idx_time1=find(v_time>=time_1);
    idx_time1=v_idx_time1(1);
end

idx_time2=find(v_time==time_2);
if isempty(idx_time2)
    v_idx_time2=find(v_time>=time_2);
   
    if isempty(v_idx_time2)
        v_idx_time2=find(v_time<=time_2);
        idx_time2=v_idx_time2(end);
    else
         idx_time2=v_idx_time2(1);
    end
end