function y = blkdiagFast(rs, cs, varargin)
%BLKDIAG  Block diagonal concatenation of matrix input arguments.
%
%                                   |A 0 .. 0|
%   Y = BLKDIAG(A,B,...)  produces  |0 B .. 0|
%                                   |0 0 ..  |
%
%   Class support for inputs: 
%      float: double, single
%
%   See also DIAG, HORZCAT, VERTCAT

% Copyright 1984-2009 The MathWorks, Inc. 
% $Revision: 1.7.4.7 $Date: 2010/09/02 13:35:29 $

% if any(~cellfun('isclass',varargin,'double'))
%     y = [];
%     for k=1:nargin
%         x = varargin{k};
%         [p1,m1] = size(y); [p2,m2] = size(x);
%         y = [y zeros(p1,m2); zeros(p2,m1) x]; %#ok
%     end
%     return
% end

% isYsparse = false;

% p2 = zeros(nargin+1,1);
% m2 = zeros(1, nargin+1);
% for k=1:nargin
%     x = varargin{k};
%     [p2(k+1),m2(k+1)] = size(x); %Precompute matrix sizes
%     if ~isYsparse&&issparse(x)
%         isYsparse = true;
%     end
% end

%Precompute cumulative matrix sizes
p1 = cumsum([0; rs]);
m1 = cumsum([0, cs]);

% if isYsparse
%     y = blkdiagmex(varargin{:});
% else
    y = zeros(p1(end),m1(end)); %Preallocate for full doubles only
    for k=1:(nargin-2)
        y(p1(k)+1:p1(k+1),m1(k)+1:m1(k+1)) = varargin{k};
    end
% end
