function M = cell2mat_amd(C)

M = NaN(size(C));
for i = 1:size(C,1)
  for j = 1:size(C,2)
    if ~isempty(C{i,j}), if ~isnumeric(C{i,j}), M(i,j) = str2num(C{i,j}); else M(i,j) = C{i,j}; end; end
  end
end

return

% AMD Depricated version

M = NaN(size(C));
ivec = find(find_cell(C));
try
  M(ivec) = cell2mat(C(ivec)); % AMD: This is incorrect! Does not convert from strings to numbers.
catch
  for i = 1:size(C,1)
    for j = 1:size(C,2)
      if ~isempty(C{i,j}), if ~isnumeric(C{i,j}), M(i,j) = str2num(C{i,j}); else M(i,j) = C{i,j}; end; end
    end
  end
end
