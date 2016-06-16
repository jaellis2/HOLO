function worked = addSource(newfilename, newsource, Np_update, setNp, data)

%%%%%%%%
%   CURRENTLY MULTIPLIES NP BY FACTOR EACH NONLINEAR ITERATION
%

NP_FACTOR = data.np_factor;



worked = 0;

fileID = fopen(newfilename, 'r');

tline = fgetl(fileID);
HO_data = cell(0,1);
while ischar(tline)
    HO_data{end+1,1} = tline;
    tline = fgetl(fileID);
end

fclose(fileID);


%%%%%%%%%
%
%   Fixing XML Lines
%

%   <Parameter name="general_source" type="TwoDArray(double)" value="3x2:{1.0, 2.0, 3.0, 4.0, 5.0, 6.0}"/>


front = '';

loc = strfind(HO_data, 'general_source');
row = find(~cellfun(@isempty, loc) == 1);

line = HO_data(row);
last = line;

for j = 1:6
    
    [first, last] = strtok(last, '"');
    
    if j < 6
        front = strcat(front,first,'"');
    end
end


%%%%%%%%%%%%%
%
%  Add Source
%

sourceStr = first;

[sfirst, oldsource] = strtok(sourceStr, '{}');


newsourceStr = sprintf('%f,' , newsource);
newsourceStr = newsourceStr(1:end-1);

newlineSrc = strcat(front, sfirst,'{',newsourceStr,'}',last);

HO_data(row) = newlineSrc;


















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(Np_update)
    
    
    loc = strfind(HO_data, 'Np');
    row = find(~cellfun(@isempty, loc) == 1);
    
    line = HO_data(row);
    last = line;
    
    front = '';
    
    for j = 1:6
        
        [first, last] = strtok(last, '"');
        
        
        if j < 6
            front = strcat(front,first,'"');
        end
        
    end
    
    Np_num = eval([ first{1} ]);

    
    %%%%%%%%%%%%%%%%%%
    %   Np FACTOR
    %
    
    
    
    newNp = round(Np_num * NP_FACTOR, 5, 'Significant');
    
    if(setNp(1) == 1)
        newNp = setNp(2);
    end
    
    
    
    newNpStr = sprintf('%i' , newNp);
   
    
    newlineNp = strcat(front, newNpStr,last);
    
    HO_data(row) = newlineNp;
    %%%%%%%%%%%%%
    %
    %  Add Np
    %
    
    
end


%%%%%%
%
%   OVERWRITE OLD FILE
%
%

fileID = fopen(newfilename, 'w');

formatSpec = '%s \n';

for i = 1:length(HO_data)
   fprintf(fileID, formatSpec, HO_data{i});
end

fclose(fileID);






worked = 1;

end


