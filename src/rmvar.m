function rmvar(ProjectPath,filename, varargin)
%RMVAR  Remove variable from MAT-File
%   RMVAR FILENAME VAR1 VAR2... removes the variables VAR1, VAR2... from
%   the Mat-File FILENAME. If FILENAME has no extension RMVAR looks for
%   FILENAME.mat
%
%   RMVAR('Filename','VAR1','VAR2' ...) does the same.
%
%   Example:
%      % Creates a file 'myfile.mat' containing 2 variables:
%      a='hello';
%      b=magic(3);
%      save myfile.mat a b
%      % Removes the variable 'a' and opens the result
%      clear
%      rmvar myfile a
%      load myfile
%      whos
%
%   F. Moisy
%   Revision: 1.00,  Date: 2008/03/31.
%
%   See also LOAD, SAVE.

% History:
% 2008/03/31: v1.00, first version.

% cd (ProjectPath)
% error(nargchk(2,inf,nargin));
% vars = rmfield(load(filename),varargin(:));
% save(filename,'-struct','vars');
