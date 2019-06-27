function uTest_Shuffle(doSpeed, doBias)
% Automatic test: Shuffle
% This is a routine for automatic testing. It is not needed for processing and
% can be deleted or moved to a folder, where it does not bother.
%
% uTest_Shuffle(doSpeed, doBias)
% INPUT:
%   doSpeed: Optional logical flag to trigger time consuming speed tests.
%            Default: TRUE. If no speed test is defined, this is ignored.
%   doBias:  Time consuming search for a bias caused by a weak shuffle algorithm
%            or a bad random number generator. SHUFFLE.MEX calls two different
%            methods for operating on the 1st or any other dimension. Therefore
%            two figures are created to test each method for each type of input
%            (although it is unlikely the algorithm works well for DOUBLEs but
%            fails for SINGLEs). The displayed lines show the relative deviation
%            from the expected mean value. For an unbiased shuffling, the lines
%            should be horizontal with additional small (<< 1.0) noise.
%            This test should be performed after changing the algorithm. It is
%            not expected that the compiler influences this test!
%            Optional, default: FALSE.
%
% OUTPUT:
%   On failure the test stops with an error.
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP
% Author: Jan Simon, Heidelberg, (C) 2010-2011 j@n-simon.de

% $JRev: R-D V:029 Sum:5AJSaZbXRC9Q Date:07-Mar-2011 00:49:06 $
% $License: BSD (see Docs\BSD_License.txt) $
% $File: Tools\UnitTests_\uTest_Shuffle.m $

% Initialize: ==================================================================
if nargin == 0
   doSpeed = true;
   doBias  = false;
elseif nargin == 1
   % Time consuming creation of graphs to detect an eventual bias:
   doBias  = false;
end

if doSpeed
   TestTime = 1;
else  % "No speed tests" mean actually fast and inprecise tests:
   TestTime = 0.1;
end

% Hello:
ErrID = ['JSim:', mfilename];
whichShuffle = which('Shuffle');

disp(['==== Test Shuffle:  ', datestr(now, 0), char(10), ...
   'Version: ', whichShuffle, char(10)]);

if doBias
   Fig1H = figure('Name', 'uTest_Shuffle: Look for bias, dim 1', ...
      'NumberTitle', 'off', ...
      'NextPlot', 'add');
   tmpH = uicontrol(Fig1H, 'Style', 'text', ...
      'String', 'Wanted: noisy but horizontal lines, values << 1.0', ...
      'FontSize', 12, ...
      'BackgroundColor', [1,0.8, 0.8], ...
      'Units', 'normalized', ...
      'Position', [0, 0, 1, 0.02]);
   tmpExt = get(tmpH, 'Extent');
   Pos    = [0, 0, 1, tmpExt(4)];
   set(tmpH, 'Position', Pos);
   
   Fig2H = figure('Name', 'uTest_Shuffle: Look for bias, dim 2', ...
      'NumberTitle', 'off', ...
      'NextPlot', 'add');
   uicontrol(Fig2H, 'Style', 'text', ...
      'String', 'Wanted: noisy but horizontal lines, values << 1.0', ...
      'FontSize', 12, ...
      'BackgroundColor', [1,0.8, 0.8], ...
      'Units', 'normalized', ...
      'Position', Pos);
end

% Compare results after seeding: -----------------------------------------------
% If the random number generator was modified, this is not necessarily an error!
Shuffle(1234567890, 'seed');
for i = 1:1e5
   Shuffle(1:3);
end
m = Shuffle(1:32);

% For exact random integer mode:
Expected1 = [18, 28, 31, 29, 9, 27, 22, 10, 6, 11, 8, 7, 16, 1, 15, 32, 14, ...
   26, 4, 23, 25, 19, 12, 24, 5, 20, 30, 13, 21, 17, 2, 3];
% For fast random integer mode:
Expected2 = [7, 2, 5, 24, 20, 22, 6, 19, 21, 25, 27, 31, 23, 10, 14, 13, 9, ...
   28, 16, 3, 12, 4, 30, 32, 18, 15, 1, 26, 17, 11, 29, 8];

if isequal(m, Expected1)
   Expected = Expected1;
elseif isequal(m, Expected2)
   Expected = Expected2;
else
   warning([ErrID, ':UnexpectedValue'], ...
      ['Unexpected values after seeding. \n', ...
      'This is not necessarily an error.\n', ...
      'Perhaps you have modified random number generator ?!']);
   Expected = m;  % Trust the result...
end

% Known answer test: -----------------------------------------------------------
ClassList = {'double', 'single', 'int32', 'int16', 'int8', 'char'};
for iClass = 1:length(ClassList)
   aClass = ClassList{iClass};
   disp(['== ', aClass]);
   R = Shuffle(mycast([], aClass));
   if ~isequal(R, mycast([], aClass))
      error(ErrID, 'SHUFFLE failed for: []');
   end
   
   R = Shuffle(mycast(1, aClass));
   if ~isequal(R, mycast(1, aClass))
      error(ErrID, 'SHUFFLE failed for: [1]');
   end
   
   R = Shuffle(mycast([1, 2], aClass));
   if ~isequal(sort(R), mycast([1, 2], aClass))
      error(ErrID, 'SHUFFLE failed for: [1, 2]');
   end
   
   R = Shuffle(mycast([1; 2], aClass));
   if ~isequal(sort(R), mycast([1; 2], aClass))
      error(ErrID, 'SHUFFLE failed for: [1; 2]');
   end
   
   for i = 1:12
      R = Shuffle(mycast([1, 2, 3], aClass));
      if ~isequal(sort(R), mycast(1:3, aClass))
         error(ErrID, 'SHUFFLE failed for: [1, 2, 3]');
      end
   end
   
   % Test is SHUFFLE creates all possible permutations for a tiny vector of
   % length 5. 120 possible permutations are possible, so I hope that creating
   % 4000 trials is enough to get all, but it is a stochastic process!!!
   % This must fail, if FACTORIAL(N) is greater than 2^32!
   nTrials = 4000;
   m = mycast(zeros(nTrials, 5), aClass);
   x = mycast(3:7, aClass);
   for i = 1:nTrials
      m(i, :) = Shuffle(x);
   end
   m = unique(m, 'rows');
   if ~all(ismember(perms(x), m, 'rows'))
      error(ErrID, 'SHUFFLE does not create all permutations.');
   end
   
   x = mycast(1:127, aClass);
   R = Shuffle(x);
   if ~isequal(x, sort(R))
      error(ErrID, 'SHUFFLE failed for: [1:127].');
   end
   
   % Test values after seeding and 3e5 times calling the KISS:
   % If SHUFFLE was compiled with LCC v2.4 (shipped with Matlab), the 64 bit
   % integer arithmetics are wrong. Or the user has implemented another RNG?
   Shuffle(1234567890, 'seed');
   for i = 1:1e5
      Shuffle(1:3);
   end
   
   m = Shuffle(mycast(1:32, aClass));
   if ~isequal(m, mycast(Expected, aClass))
      warning([ErrID, ':UnexpectedValue'], ...
         'Unexpected values after seeding. Compiled with LCC v2.4 ?!');
   end
   
   Shuffle([now, rand*4294967295, abs(cputime), rand*4294967295], 'seed');
   
   x = mycast(1:127, aClass);
   R = Shuffle(x, 1);
   if ~isequal(x, R)
      error(ErrID, 'SHUFFLE failed for: ([1:127], 1).');
   end
   
   % No multi-dimensional CHAR array:
   x = reshape(1:(3*4*5*6*7*8), 3:8);
   for iDim = 1:ndims(x)
      r = Shuffle(x, iDim);
      if ~isequal(x, sort(r, iDim))
         error(ErrID, 'SHUFFLE failed multi-dim call.');
      end
   end
   
   disp('  ok: known answer tests');
   
   % Search for a bias - this is actually not needed after the algorithm was
   % tested by the author exhaustively. But do this again with n=[2,3,4,8,32] if
   % the method is modifed!
   if doBias
      disp('== Draw diagrams to look for bias');
      n  = 8;  % Number of points to shuffle
      r  = 3;  % Number of repetitions, so we have r * n * 1e5 tests
      P1 = zeros(r, n, n);
      P2 = zeros(r, n, n);
      x  = repmat(1:n, n, 1);
      y  = transpose(x);
      for k = 1:r
         m1 = zeros(n, n, 1e5);
         m2 = zeros(n, n, 1e5);
         for i = 1:1e5
            m1(:, :, i) = Shuffle(y, 1);
            m2(:, :, i) = Shuffle(x, 2);
         end
         aP1         = transpose(sum(m1, 3));
         aP2         = transpose(sum(m2, 3));
         P1(k, :, :) = aP1 / 1e5 - mean(1:n);
         P2(k, :, :) = aP2 / 1e5 - mean(1:n);
         drawnow;
      end
      
      figure(Fig1H);
      subplot(2, 3, iClass);
      plot(reshape(P1, n, r*n));
      xlim([1, n]);
      title(aClass);
      
      figure(Fig2H);
      subplot(2, 3, iClass);
      plot(reshape(P2, n, r*n));
      xlim([1, n]);
      title(aClass);
      drawnow;
   end
end

% Check index method: ----------------------------------------------------------
fprintf('\n== Shuffle(n, ''index'')\n');

r = Shuffle(0, 'index');
if ~isempty(r)
   error(ErrID, 'SHUFFLE(0, ''index'') failed!');
end

r = Shuffle(1, 'index');
if ~isequal(r, uint8(1))
   error(ErrID, 'SHUFFLE(1, ''index'') failed!');
end

r = Shuffle(2, 'index');
if ~isequal(r, uint8(1:2)) && ~isequal(r, uint8([2, 1]))
   error(ErrID, 'SHUFFLE(2, ''index'') failed!');
end

r = Shuffle(5, 'index');
if ~isequal(sort(r), uint8(1:5))
   error(ErrID, 'SHUFFLE(5, ''index'') failed!');
end

r = Shuffle(255, 'index');
if ~isequal(sort(r), uint8(1:255))
   error(ErrID, 'SHUFFLE(255, ''index'') failed!');
end

r = Shuffle(256, 'index');
if ~isequal(sort(r), uint16(1:256))
   error(ErrID, 'SHUFFLE(256, ''index'') failed!');
end

r = Shuffle(65535, 'index');
if ~isequal(sort(r), uint16(1:65535))
   error(ErrID, 'SHUFFLE(65535, ''index'') failed!');
end

r = Shuffle(70000, 'index');
if ~isequal(sort(r), uint32(1:70000))
   error(ErrID, 'SHUFFLE(70000, ''index'') failed!');
end

% Needs 35 GB memory:
% r = Shuffle(4294967295, 'index');
% if ~isequal(sort(r), uint32(1:4294967295))
%    error(ErrID, 'SHUFFLE(4294967295, ''index'') failed!');
% end

r = Shuffle(0, 'index', 0);
if ~isempty(r)
   error(ErrID, 'SHUFFLE(0, ''index'', 1) failed!');
end

r = Shuffle(1, 'index', 1);
if ~isequal(r, uint8(1))
   error(ErrID, 'SHUFFLE(1, ''index'', 1) failed!');
end

r = Shuffle(2, 'index', 1);
if ~(isequal(r, uint8(1)) || isequal(r, uint8(2)))
   error(ErrID, 'SHUFFLE(2, ''index'', 1) failed!');
end

r = Shuffle(5, 'index', 4);
if length(unique(r)) ~= 4 || ~all(ismember(r, uint8(1:5)))
   error(ErrID, 'SHUFFLE(5, ''index'', 4) failed!');
end

r = Shuffle(256, 'index', 100);
if length(unique(r)) ~= 100 || ~all(ismember(r, uint16(1:256)))
   error(ErrID, 'SHUFFLE(256, ''index'', 100) failed!');
end

r = Shuffle(70000, 'index', 100);
if length(unique(r)) ~= 100 || ~all(ismember(r, uint32(1:70000)))
   error(ErrID, 'SHUFFLE(70000, ''index'', 100) failed!');
end

disp('  ok: Shuffle(n, ''index'')');

% Check derange method: --------------------------------------------------------
fprintf('\n== Shuffle(n, ''derange'')\n');

r = Shuffle(2, 'derange');
if ~isequal(r, uint8([2, 1]))
   error(ErrID, 'SHUFFLE(2, ''derange'') failed!');
end

r = Shuffle(5, 'derange');
if ~isequal(sort(r), uint8(1:5)) || any(r == 1:5)
   error(ErrID, 'SHUFFLE(5, ''derange'') failed!');
end

r = Shuffle(255, 'derange');
if ~isequal(sort(r), uint8(1:255)) || any(r == 1:255)
   error(ErrID, 'SHUFFLE(255, ''derange'') failed!');
end

r = Shuffle(256, 'derange');
if ~isequal(sort(r), uint16(1:256)) || any(r == 1:256)
   error(ErrID, 'SHUFFLE(256, ''derange'') failed!');
end

r = Shuffle(65535, 'derange');
if ~isequal(sort(r), uint16(1:65535)) || any(r == 1:65535)
   error(ErrID, 'SHUFFLE(65535, ''derange'') failed!');
end

r = Shuffle(70000, 'derange');
if ~isequal(sort(r), uint32(1:70000)) || any(r == 1:70000)
   error(ErrID, 'SHUFFLE(70000, ''derange'') failed!');
end

% Needs 35 GB memory and some patients:
% r = Shuffle(4294967296, 'derange');
% if ~isequal(sort(r), uint32(1:4294967296))
%    error(ErrID, 'SHUFFLE(4294967296, ''derange'') failed!');
% end
% for i = 1:4294967296
%   if r(i) == i
%     error(ErrID, 'SHUFFLE(4294967296, ''derange'') failed!');
%   end
% end

r = Shuffle(2, 'derange', 0);
if ~isempty(r)
   error(ErrID, 'SHUFFLE(2, ''derange'', 0) failed!');
end

r = Shuffle(2, 'derange', 1);
if ~isequal(r, uint8(2))
   error(ErrID, 'SHUFFLE(2, ''derange'', 1) failed!');
end

r = Shuffle(2, 'derange', 2);
if ~isequal(r, uint8([2, 1]))
   error(ErrID, 'SHUFFLE(2, ''derange'', 2) failed!');
end

r = Shuffle(5, 'derange', 4);
if length(unique(r)) ~= 4 || ~all(ismember(r, uint8(1:5))) ...
      || any(r == 1:4)
   error(ErrID, 'SHUFFLE(5, ''derange'', 4) failed!');
end

r = Shuffle(256, 'derange', 100);
if length(unique(r)) ~= 100 || ~all(ismember(r, uint16(1:256))) ...
      
   error(ErrID, 'SHUFFLE(256, ''derange'', 100) failed!');
end

r = Shuffle(70000, 'derange', 100);
if length(unique(r)) ~= 100 || ~all(ismember(r, uint32(1:70000)))
   error(ErrID, 'SHUFFLE(70000, ''derange'', 100) failed!');
end

for Len = [4, 255, 256, 65535, 65536]
   iTime = cputime;
   while cputime - iTime < 1.0
      r1 = Shuffle(Len, 'derange');
      r2 = Shuffle(Len, 'derange', Len - 1);
      if any(r1 == 1:Len) || any(r2 == (1:Len - 1))
         error(ErrID, 'SHUFFLE(N, derange) is not an derangement!');
      end
   end
end

disp('  ok: Shuffle(n, ''derange'')');

% Provoke errors: --------------------------------------------------------------
fprintf('\n== Catching errors:\n');

tooLazy = 0;
try
   r = Shuffle(1:10, 1.5);
   tooLazy = 1;
catch
   disp('  ok: Shuffle(1:10, 1.5) failed as wanted.');
end
if tooLazy
   error(ErrID, 'SHUFFLE(1:10, 1.5) did not fail.');
end

try
   r = Shuffle(2, 'index', 3, '?');
   tooLazy = 1;
catch
   disp('  ok: Shuffle(2, index, 3, 4thInput) failed as wanted.');
end
if tooLazy
   error(ErrID, 'SHUFFLE(2, index, 3, 4thInput) did not fail.');
end

try
   r = Shuffle(-1, 'index');
   tooLazy = 1;
catch
   disp('  ok: Shuffle(-1, index) failed as wanted.');
end
if tooLazy
   error(ErrID, 'SHUFFLE(-1, index) did not fail.');
end

try
   r = Shuffle(2, 'index', 3);
   tooLazy = 1;
catch
   disp('  ok: Shuffle(2, index, 3) failed as wanted.');
end
if tooLazy
   error(ErrID, 'SHUFFLE(2, index, 3) did not fail.');
end

try
   r = Shuffle(1:2, 'index');
   tooLazy = 1;
catch
   disp('  ok: Shuffle(1:2, index) failed as wanted.');
end
if tooLazy
   error(ErrID, 'SHUFFLE(1:2, index) did not fail.');
end

try
   r = Shuffle(1:5, 'derange');
   tooLazy = 1;
catch
   disp('  ok: Shuffle(1:5, derange) failed as wanted.');
end
if tooLazy
   error(ErrID, 'SHUFFLE(1:5, derange) did not fail.');
end

try
   r = Shuffle(1, 'derange');
   tooLazy = 1;
catch
   disp('  ok: Shuffle(1, derange) failed as wanted.');
end
if tooLazy
   error(ErrID, 'SHUFFLE(1, derange) did not fail.');
end

try
   r = Shuffle(1, 'fitzliputzli');
   tooLazy = 1;
catch
   disp('  ok: Shuffle(1, "unknownString") failed as wanted.');
end
if tooLazy
   error(ErrID, 'SHUFFLE(1, "unknownString") did not fail.');
end

% Speed: -----------------------------------------------------------------------
fprintf('\n== Compare speed:\n');

for aDataLen = [10, 100, 1000, 10000, 100000, 1000000]
   Data = 1:aDataLen;
   
   % Determine number of loops:
   iTime = cputime;
   iLoop = 0;
   while cputime - iTime < TestTime
      v = Data(randperm(aDataLen));  %#ok<*NASGU>
      clear('v');   % Suppress JIT acceleration for realistic times
      iLoop = iLoop + 1;
   end
   nDigit = max(1, floor(log10(max(1, iLoop))) - 1);
   nLoop  = max(4, round(iLoop / 10 ^ nDigit) * 10 ^ nDigit);
   
   tic;
   for i = 1:nLoop
      v = Data(randperm(aDataLen));
      clear('v');
   end
   RandPermTime = toc;
   
   tic;
   for i = 1:nLoop
      v = Shuffle(Data);
      clear('v');
   end
   MexTime = toc;
   
   tic;
   for i = 1:nLoop
      v = Shuffle_M(Data);
      clear('v');
   end
   MTime = toc;
   
   tic;
   for i = 1:nLoop
      v = Data(Shuffle(aDataLen, 'index'));
      clear('v');
   end
   MexIndexTime = toc;
   
   tic;
   for i = 1:nLoop
      v = Data(Shuffle(aDataLen, 'derange'));
      clear('v');
   end
   MexDerangeTime = toc;
   
   fprintf('[1 x %d],  %d loops:\n', aDataLen, nLoop);
   fprintf('  X(RANDPERM):         %6.2f sec\n', RandPermTime);
   fprintf('  SHUFFLE(X) Matlab:   %6.2f sec\n', MTime);
   fprintf('  X(SHUFFLE(index)):   %6.2f sec\n', MexIndexTime);
   fprintf('  X(SHUFFLE(derange)): %6.2f sec\n', MexDerangeTime);
   fprintf('  SHUFFLE(X):          %6.2f sec  ==>  %.1f%% of RANDPERM\n', ...
      MexTime, 100.0 * MexTime / RandPermTime);
end

% Test partial index vector:
fprintf('\n');
DataLen = [10, 1000, 100000, 1000000];
PartLen = [4,  100,  1000,   1000];
for iData = 1:length(DataLen)
   aDataLen = DataLen(iData);
   nOut     = PartLen(iData);
   
   % Determine number of loops:
   iTime = cputime;
   iLoop = 0;
   while cputime - iTime < TestTime
      v = randperm(aDataLen);  %#ok<*NASGU>
      v = v(1:nOut);
      clear('v');   % Suppress JIT acceleration for realistic times
      iLoop = iLoop + 1;
   end
   nDigit = max(1, floor(log10(max(1, iLoop))) - 1);
   nLoop  = max(4, round(iLoop / 10 ^ nDigit) * 10 ^ nDigit);
   
   tic;
   for i = 1:nLoop
      v = randperm(aDataLen);
      v = v(1:nOut);
      clear('v');
   end
   MTime = toc;
   
   tic;
   for i = 1:nLoop
      v = Shuffle(aDataLen, 'index', nOut);
      clear('v');
   end
   MexTime = toc;
   
   fprintf('Choose %d from [1 x %d],  %d loops:\n', nOut, aDataLen, nLoop);
   fprintf('  RANDPERM(N), [1:nOut]:        %6.3f sec\n', MTime);
   fprintf('  SHUFFLE(N, index, nOut)) Mex: %6.3f sec', MexTime);
   fprintf('  ==>  %.1f%% of RANDPERM\n', 100.0 * MexTime / MTime);
end

fprintf('\nSHUFFLE seems to work well.\n');

return;

% ******************************************************************************
function A = mycast(A, ClassName)
% Simulate CAST for Matlab 6
A = feval(ClassName, A);

return;

% ******************************************************************************
function X = Shuffle_M(X)
% Fast X(RANDPERM(LENGTH(X))), Knuth's shuffle
% This Matlab version is faster than RANDPERM, but the MEX is faster.
% Author: Jan Simon, Heidelberg, (C) 2010-2011 j@n-simon.de

for i = 1:numel(X)
   w    = ceil(rand * i);
   t    = X(w);
   X(w) = X(i);
   X(i) = t;
end

return;
