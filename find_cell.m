function resvec = find_cell(cellmat)

tmpmat = zeros(size(cellmat));
for i = 1:size(tmpmat,1)
  for j = 1:size(tmpmat,2)
    tmpmat(i,j) = ~isempty(cellmat{i,j});
  end
end
%resvec = find(tmpmat);
resvec = tmpmat;
