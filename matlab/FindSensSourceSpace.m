function [sensorspace,sourcespace] = FindSensSourceSpace(fwd)
c = struct2cell(fwd);
sourcespace = c{9};
channelinfo = c{11};
sensorspace = zeros(306,3);
for i = (1:306)
    temp = channelinfo(i).loc;
    chaloc = temp(1:3)';
    sensorspace(i,:) = chaloc;
end


end