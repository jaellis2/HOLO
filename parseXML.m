function theStruct = parseXML(filename, newfilename, holo_elem, mat_elem)



theStruct = struct;

%%%%%%%%%%%%%%%%%%%%%%%%
%   READ HOLO FILE
%

fileID = fopen(filename, 'r');

tline = fgetl(fileID);
HO_data = cell(0,1);
while ischar(tline)
    HO_data{end+1,1} = tline;
    tline = fgetl(fileID);
end

fclose(fileID);

%%%%%%%%%%%%%%%%
%   STORE HOLO DATA
%

 values = cell(length(holo_elem),1);

for i = 1:length(holo_elem)
    
    loc = strfind(HO_data, char(holo_elem(i)));
    row = find(~cellfun(@isempty, loc) == 1);

    line = HO_data(row(1));
    last = line;
    for j = 1:6
        
    [first, last] = strtok(last, '"');
    
    end
    
    values(i) = first;
end

v1_cell = regexprep( values(1), '[{;}]', '');
v1_str = v1_cell{1};
v1_int = eval([ '[', v1_str, ']' ]);

theStruct.tau = v1_int(end) - v1_int(1);
theStruct.xgrid = v1_int';
theStruct.Np = eval(values{3});
theStruct.prob_name = values{4};

% v2_cell = regexprep( values(2), '[{;}]', '');
% v2_str = v2_cell{1};
% v2_int = eval([ '[', v2_str, ']' ]);
% 
% 
% theStruct.source = v1_int(2) - v1_int(1);



library_filename = values{2};

%%%%%%%%%%%%%%%%%%%%%%%%
%   READ HOLO FILE
%

fileID = fopen(library_filename, 'r');

tline = fgetl(fileID);
mat_data = cell(0,1);
while ischar(tline)
    mat_data{end+1,1} = tline;
    tline = fgetl(fileID);
end

fclose(fileID);

 values = cell(length(holo_elem),1);

%%%%%%%%%%%%%%%%
%   STORE HOLO DATA
%
 
mat_values = cell(length(mat_elem),1);
 
 
for i = 1:length(mat_elem)
    
    loc = strfind(mat_data, char(mat_elem(i)));
    row = find(~cellfun(@isempty, loc) == 1);

    line = mat_data(row);
    last = line;
    for j = 1:6
        
    [first, last] = strtok(last, '"');
    
    end
    
    mat_values(i) = first;
end


mv1_cell = regexprep( mat_values(1), '[{;}]', '');
mv1_str = mv1_cell{1};
mv1_int = eval([ '[', mv1_str, ']' ]);


theStruct.sig_t = mv1_int;

[first2, mv_2] = strtok(mat_values(2), ' ');
[mv_2, last2] = strtok(mv_2, '}');

theStruct.sig_s = eval([ '[',mv_2{1},']' ]);

%%%%%%%%%%%%%%%
%   Create running xml file
%

fileID = fopen(newfilename, 'w');

formatSpec = '%s \n';

for i = 1:length(HO_data)
   fprintf(fileID, formatSpec, HO_data{i});
end

fclose(fileID);





end