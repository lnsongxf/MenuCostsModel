function [xi hazard survival] = hazard_discrete(durations)

    % compute hazard only at points in durations
    
    % first compute survival function
    points   = 0:9;
    survival = zeros(length(points),1);
    for i=1:length(points)
        point       = points(i);
        temp        = durations(durations > point);
        survival(i) = length(temp)/length(durations);
    end

    % hazard function
    hazard = 1 - survival(2:end)./survival(1:end-1);
    xi     = points(2:end);
    
end