function      [data, result]= readtext(fname, delimiter, comment, quotes, options)

%  Usage:     [data, result]= readtext(fname, delimiter, comment, quotes, options)
% 
% Whatever text file you give it, readtext returns an array of the contents (or send me a 
%   bug report). Matlab can't read variable length lines or variable type values with the standard 
%   library. readtext can read any text file. Any string (or even regexp) can be delimiting, 
%   default is a comma. Everything after (and including) a comment character, until the line end, 
%   is ignored. Quote characters may also be given, everything between them is treated as one item. 
%   There are options to control what will be converted to numbers and how empty items are saved. 
% 
% If you find any errors, please let me know: peder at axensten dot se
% 
% fname:      the file to be read.
% 
% delimiter:  (default: ',') any string. May be a regexp, but this is a bit slow on large files. 
% 
% comment:    (default: '') zero or one character. Anything after (and including) this character, 
%   until the end of the line, will be ignored. 
% 
% quotes:     (default: '') zero, one (opening quote equals closing), or two characters (opening 
%   and closing quote) to be treated as paired braces. Everything between the quotes will be 
%   treated as one item. The quotes will remain. Quotes may be nested.
% 
% options:    (default: '') may contain (concatenate combined options): 
% - 'textual': no numeric conversion ('data' is a cell array of strings only), 
% - 'numeric': everything is converted to a number or NaN ('data' is a numeric array, empty items 
%   are converted to NaNs unless 'empty2zero' is given), 
% - 'empty2zero': an empty field is saved as zero, and 
% - 'empty2NaN': an empty field is saved as NaN. 
% - 'usewaitbar': call waitbar to report progress. If you find the wait bar annoying, get 'waitbar 
%   alternative' at http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=11398
% 
% data:       A cell array containing the read text, divided into cells by delimiter and line 
%   endings. 'data' will be empty if the file is not found, could not be opened, or is empty. 
%   With the option 'numeric', 'data' will be a numeric array, with 'textual', 'data' will be a 
%   cell array of strings only, and otherwise it will be a mixed cell array. For Matlab < version 7, 
%   returned strings may contain leading white-space.
% 
% result:     a structure:
% .min: minimum number of columns found in a line.
% .max: number of columns in 'data', before removing empty columns.
% .rows: number of rows in 'data', before removing empty rows. 
% .numberMask: true, if numeric conversion ('NaN' converted to NaN counts).
% .number: number of numeric conversions ('NaN' converted to NaN counts).
% .emptyMask: true, if empty item in file.
% .empty: number of empty items in file.
% .stringMask: true, if non-number and non-empty.
% .string: number of non-number, non-empty items.
% .quote: number of quotes. 
% 
% INSPIRATION: loadcell.m (id 1965). The goal of readtext is to be at least as flexible (you be 
%   the judge) and quicker. On my test file (see below), readtext is about 3--4 times 
%   as quick, maybe even more on large files. In readtext you may use a regexp as 
%   delimiter and it can ignore comments in the text file. 
% 
% SPEED:      Reading a 1MB file (150000 items!) with 'numeric' takes about 100 seconds on a 
%   fairly slow system. Time scales approximately linearly with input file size. 
% - Conversion from string to numeric is slow (I can't do anything about this), but using the 
%   option 'textual' is a lot quicker (above case takes 12 seconds).
% - Using a regexp delimiter is slower (during initializing), it adds 250 seconds! 
% 
% EXAMPLE:    [a,b]= readtext('txtfile', '[,\t]', '#', '"', 'numeric-empty2zero')
% This will load the file 'txtfile' into variable a, treating any of tab or comma as
%   delimiters. Everything from and including # to the next newline will be ignored. 
%   Everything between two double quotes will be treated as a string. Everything will 
%   be converted to numbers and a numeric array returned. Non-numeric items will become 
%   NaNs and empty items are converted to zero. 
% 
% Copyright (C) Peder Axensten (peder at axensten dot se), 2006.

% HISTORY:
% Version 1.0, 2006-05-03.
% Version 1.1, 2006-05-07:
% - Made 'options' case independent. 
% - Removed result.errmess -- now use error(errmess) instead. 
% - Removed result.nan -- it was equivalent to result.string, so check this instead.
% - Added some rows', 'result' fields: 'numberMask', 'emptyMask', and 'stringMask' 
%   (see 'result:', above).
% - A few small bug fixes.
% Version 1.2, 2006-06-06:
% - Now works in Matlab 6.5.1 (R13SP1) (maybe 6.5 too), versions <6.5 will NOT work.
% Version 1.3, 2006-06-20:
% - Better error report when file open fails. 
% - Somewhat quicker. 
% - Recommends 'waitbar alternative'. Ok with Matlab orig. waitbar too, of course. 
% Version 1.4, 2006-07-14:
% - Close waitbar instead of deleting it, and some other minor waitbar compatibility fixes. 
% Version 1.5, 2006-08-13:
% - No more (La)TeX formatting of file names. 
% - Prefixed waitbar messages with '(readtext)'. 
% Version 1.6, 2006-10-02:
% - Better removal of comments. Could leave an empty first row before. 
% - Added a 'usewaitbar' option. 
% - Now removes empty last columns and rows. 
% 
% TO DO:
% - Write a better text macro. 
% - Add result.quoteMask.
% - Add 'removeemptycolumns' and 'removeemptyrows' options. 

% KEYWORDS:     import, read, load, text, delimited, cell, numeric, array, flexible
% 
% REQUIREMENTS: Works in Matlab 6.5.1 (R13SP1) (probably 6.5 too), versions <6.5 will NOT work.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% Read (or set to default) the input arguments.
	if((nargin < 1) || ~ischar(fname) || isempty(fname))		% Is there a file name?
		error('First argument must be a file name!'); 
	end
	if(nargin < 2), delimiter=	',';				end			% Default delimiter value.
	if(nargin < 3), comment=	'';					end			% Default comment value.
	if(nargin < 4), quotes=		'"';					end			% Default quotes value.
	if(nargin < 5), options=	[];					end			% Default options value.
	
	options=		lower(options);
	op_waitbar=		~isempty(strfind(options, 'usewaitbar'));	% Do waitbar calls. 
	op_numeric=		~isempty(strfind(options, 'numeric'));		% Read as numerical. 
	op_textual=		~isempty(strfind(options, 'textual')) && ~op_numeric;	% Read as textual. 
	op_empty=		[];											% Ignore empties, ...
	if(~isempty(strfind(options, 'empty2zero')))
		op_empty=		0;										% ... or replace by zero ...
	elseif(op_numeric || ~isempty(strfind(options, 'empty2nan')))
		op_empty=		NaN;									% ... or replace by NaN.
	end
	if(op_textual), op_empty= num2str(op_empty);	end			% Textual 'empty'.
	if(~ischar(comment) || (length(comment) > 1))
		error('Argument ''comment'' must be a string of maximum one character.');
	end
	if(~ischar(quotes) || (length(quotes) > 2))
		error('Argument ''quotes'' must be a string of maximum two characters.');
	end
	
	% Set the default return values.
	result.min=		Inf;
	result.max=		0;
	result.quote=	0;
	
	% Read the file.
	[fid, errmess]=	fopen(fname, 'r');							% Open the file.
	if(fid < 0), error(['Trying to open ' fname ': ' errmess]); end
	text=			fread(fid, 'uchar=>char')';					% Read the file.
	fclose(fid);												% Close the file.

        if 0
          for posi = 1:200
            tmp_char = text(posi);
            tmp_int8 = int8(tmp_char);
            fprintf(1,'%3d: %3d',posi,tmp_int8);
            if tmp_int8>13 & tmp_int8<128
              fprintf(1,' %3s',tmp_char);
            end
            fprintf(1,'\n');
          end
        end

	eol = char(10);
        eol_dos = [char(13) char(10)];
        if ~isempty(find(findstr(text,eol_dos))) % Replace DOS eol string
          text = strrep(text,eol_dos,eol);
        end
        eol_win = char(13);
        if ~isempty(find(findstr(text,eol_win))) % Replace Windows(?) eol string
          text = strrep(text,eol_win,eol);
        end

        % First, determine dimensions of array
        eolvec = (text==eol);
        if ~eolvec(end), eolvec = [eolvec 1]; end
        eolindvec = find(eolvec);
       
        nrows = length(eolindvec);
        commentvec = zeros(1,nrows);
        ncolvec = zeros(1,nrows);
        lp1 = 1;
        for n = 1:nrows
          lp2 = eolindvec(n)-1;
          linestr = text(lp1:lp2);
          if strcmp(' ',delimiter)
            linestr = strtrim(strrep(linestr,char(9),' ')); % Treat TAB as space
            linestr = linestr(min(find(linestr~=' ')):end);
            while ~isempty(strfind(linestr,'  '))
              linestr = strrep(linestr,'    ',' ');
              linestr = strrep(linestr,'  ',' ');
            end
          end
          if (length(linestr)>=length(comment) & strcmp(comment,linestr(1:length(comment)))) 
            commentvec(n) = true;
          else
            quotevec = ismember(linestr,quotes);  
            quotedvec = (mod(cumsum(quotevec),2)==1) | (quotevec>0);
            delimvec = ismember(linestr,delimiter) & ~quotedvec;
            delimvec = [delimvec 1];
            delimindvec = find(delimvec);
            ncols = length(delimindvec);
            ncolvec(n) = ncols;
%            fprintf(1,'%d: %s\n',n,linestr);
            if n==1
              fprintf('ncols=%d\n',ncols);
              fprintf(1,'%s\n',linestr);
            elseif ncolvec(n)~=ncolvec(1)
              fprintf('ncolvec(%d)=%d\n',n,ncols);
              fprintf(1,'linestr:%s\n',linestr);
              for posi = 1:200
                tmp_char = text(posi);
                tmp_int8 = int8(tmp_char);
                fprintf(1,'%3d: %3d',posi,tmp_int8);
                if tmp_int8>13 & tmp_int8<128
                  fprintf(1,' %3s',tmp_char);
                end
                fprintf(1,'\n');
              end
              keyboard
            end
          end
          lp1 = eolindvec(n)+1;
          if mod(n,100000)==0, fprintf(1,'n=%6d of %6d (%6.2f%%)\n',n,nrows,100*n/nrows); end
        end
        ncols = max(ncolvec);
        ncomments = length(find(commentvec));
        data = cell(nrows-ncomments,ncols);

        % Now, read data into cell array
        lp1 = 1;
        rownum = 0;
        for n = 1:nrows
          lp2 = eolindvec(n)-1;
          linestr = text(lp1:lp2);
          if strcmp(' ',delimiter)
            linestr = strrep(linestr,char(9),' '); % Treat TAB as space
            linestr = linestr(min(find(linestr~=' ')):end);
            while ~isempty(strfind(linestr,'  '))
              linestr = strrep(linestr,'    ',' ');
              linestr = strrep(linestr,'  ',' ');
            end
          end
          if ~commentvec(n)
            rownum = rownum+1;
            quotevec = ismember(linestr,quotes);
            quotedvec = (mod(cumsum(quotevec),2)==1) | (quotevec>0);
            delimvec = ismember(linestr,delimiter) & ~quotedvec;
            delimvec = [delimvec 1];
            delimindvec = find(delimvec);
            cp1 = 1;
            for m = 1:ncols
              cp2 = delimindvec(m)-1;
              tmp = linestr(cp1:cp2);
              data{rownum,m} = tmp;
              cp1 = delimindvec(m)+1;
            end
          end
          lp1 = eolindvec(n)+1;
          if mod(n,100000)==0, fprintf(1,'n=%6d of %6d (%6.2f%%)\n',n,nrows,100*n/nrows); end
        end
end
