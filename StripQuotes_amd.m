function  D = StripQuotes_amd(Din,bc_flag)
%function  D = StripQuotes_amd(Din,bc_flag)
%
% copied from ~dale/matlab/ADNI_DB on 11/16/09
%

if ~exist('bc_flag','var'), bc_flag = true; end

D = Din;
for i = 1:size(D,1)
  for j = 1:size(D,2)
    if ischar(D{i,j})
      s = strtrim(D{i,j});
      if length(s)>=2
        if s(1) == '"' & s(end) == '"'
          D{i,j} = strtrim(s(2:end-1));
        end
      end
      if isempty(D{i,j})
        D{i,j} = [];
      else
        if bc_flag % Convert elements to numeric, if possible
          tmp = str2double(D{i,j});
          if isfinite(tmp)
            D{i,j} = tmp;
          end
        end
      end
    end
  end
end
