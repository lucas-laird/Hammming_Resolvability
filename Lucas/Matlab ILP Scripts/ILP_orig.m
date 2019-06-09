%{
 April 2019, inspired by seeing Lucas Laird's thesis,
"Metric Dimension of Hamming Graphs and Applicationsn to Computational Bio"
(w. Manuel Lladser and Richard "Carter" Tilquist)
Stephen Becker, April 25 2019

To run this code, 
 (1) install Matlab
 (2) install CVX and the Gurobi package -- instructions for doing these
both together (and getting the two free academic licenses you'll need)
are at this website:
http://web.cvxr.com/cvx/doc/gurobi.html


For small sets, we have an exact deterministic strategy, c   = 2.^(1:d)
For larger sets, we have a probability 1 strategy, c = randn(1,d)
  where d = a*k

There are also two variants:

(1) solve max abs( c'*x ) s.t. x is feasible
  (and we break this into two integer linear programs (ILP)s,
   one for c'*x and one for -c'*x, since max abs( ... ) is not convex )

(2) solve the feasibility problem (e.g., max 0, s.t. x is feasible)
  with the constraints that abs( c'*x ) >= 1e-8 (or some small epsilon)
  (and again, break this into two ILPs)

  The advantage of this 2nd formulation is that if we are trying to
  certify that the set is *not* resolving, we just need a nonzero
  feasible solution, and finding a single feasible point in (2)
  is faster than solving (1) which unnecessarily tries to find
  a maximum.

  With the deterministic c (for small "d"), version (2) is still correct.
  With a random c, version (2) can still certify when a set is not
  resolving, but if it says that a set *is* resolving, this is only true
  with a given probability (say, 99%, depending on the exact value
  of epsilon). So in this case, we'd need to re-run the algorithm a few
  times to get the uncertainty to a "safe" level (like 0.0000001% chance
  of being wrong).


Actually, I realized that our feasible region is symmetric about 0,
 in the sense that if x is feasible, so is -x. Therefore we 
 don't really need to do the max abs(c'*x), we can just remove
 the absolute value!!

%}

%% Load the known length 82 resolving set for H_{8,20}
fid     = fopen('octamers.txt','r');
str     = textscan(fid,'%s');
str     = str{1};
fclose(fid);

nNodes = length(str);
k   = 8;
a   = 20;
A   = zeros( nNodes, k*a );
alphabet = unique( strcat( str{:} ) );
if length( alphabet ) ~= a, disp('Error, too many alphabet entries'); end
% Make a lookup for the one-hot encoding
alphabetNumeric = double( alphabet );
encoding        = zeros( max(alphabetNumeric), 1 );
encoding( alphabetNumeric ) = 1:a;

for i= 1:nNodes
    s   = str{i};
    z   = zeros(a,k);
    for j = 1:k
        c = s(j)+0;
        z( encoding(c), j ) = 1;
    end
    A(i,:) = z(:)';
end
A_82 = A;



%% Choose your problem

% k   = 3; % length k sequences (k-mers)
% a   = 2; % means binary alphabet if a=2
% 
% % This is the example H_{3,2} from page 5 of Lucas' thesis
% A   = [0,1,1,0,1,0;
%        0,1,1,0,0,1;
%        1,0,1,0,0,1];a=2; k=3;
%    
% % Or an example that is not resolving...
% A   = [0,1,1,0,1,0;
%        0,1,1,0,0,1;
%        1,0,0,1,0,1];a=2; k=3; % change 3rd vertex to "010" instead of "001"

% == Or, try scaling it up to something larger, like H_{8,20}
% %   with, say, 80 nodes (or 10 nodes)
%   21 nodes, rng(3), is NOT resolving.
nNodes = 16;
% nNodes  = 30;
% nNodes  = 40;
k   = 8;
a   = 20;
% rng(0);
rng(3);
A   = randi( [0,1], nNodes, k*a );

% == Try the known H_{8,20} one
% A   = A_82;
% % % A   = A_82(1:81,:); % NOT a resolving set!
% % A   = A_82(6:82,:); % IS a resolving set!
% A   = A_82(7:82,:); % NOT a resolving set!
% nNodes = size(A,1);
% k   = 8;
% a   = 20;




%% Solve
fprintf('\n=== RUNNING ILP TO CHECK IF SET IS RESOLVING ===\n');
d   = a*k;
mat = @(x) reshape( x, a, k ); % turns (a*dk x 1 into a x k )

% Pick a strategy
if d <= 16
    c   = 2.^(1:d);  % for d > 16, we have numerical issues with this...
else
    % Probability 1
    c = randn(1,d);
end

% if x has any non-zero entries, in range {-1,0,1},
%   then dot( c, x ) must be non-zero because there
%   are no cancellations.
% Thus, if we have a non-zero feasible point x,
%   then either dot(c,x) is > 0 or it is < 0
% So, we have two cases
%   (We have to do it this way since maximize sum(abs(x))
%                           or even  maximize dot(c,x)
%   isn't convex).

tic
cvx_begin quiet
    cvx_solver gurobi % handles integer constraints
    variable x(d) integer
    %minimize c*x
    minimize 0
    subject to
        abs(x) <= 1 % so {-1,0,1} range
        sum( mat(x) ) == 0 % Want blocks to sum to 0
        sum( mat(abs(x) ) ) <= 2
        A*x == 0
        c*x <= -1e-3 % quickly find certificate
cvx_end
toc
if strcmpi( cvx_status, 'infeasible' )
    resolving = true;
else
    if norm( x ) < 1e-8
        disp('Our certificate seems to have numerical problems');
    end
    resolving = false;
end

% If we claimed to be resolving, do a final check
if resolving
    disp('Feasiblity problem suggests this is a resolving set; running slower ILP to double-check');
    if d > 16
        % Draw a new random variable
        c = randn(1,d);
    end
    tic
    cvx_begin quiet
        cvx_solver gurobi % handles integer constraints
        variable x(d) integer
        minimize c*x
        subject to
            abs(x) <= 1 % so {-1,0,1} range
            sum( mat(x) ) == 0 % Want blocks to sum to 0
            sum( mat(abs(x) ) ) <= 2
            A*x == 0
    cvx_end
    toc
    if ~all( x== 0 )
        disp('Slower ILP found a non-zero solution, so this is actually not resolving');
        resolving = false;
    end
end


if ~resolving
    fprintf(2,'Found a certificate that this is NOT a resolving set\n');
    X   = mat(x)
else
    fprintf(2,'Only solution is the zero solution, so this IS a resolving set\n');
    X   = mat(x);
end


% Realized that we don't need to loop over sgn
% resolving = true;
% for sgn = [-1,1]
%   tic
%   cvx_begin quiet
%     cvx_solver gurobi % handles integer constraints
%     variable x(d) integer
%     %minimize sgn*(c*x)
%     minimize 0
%     subject to
%         abs(x) <= 1 % so {-1,0,1} range
%         sum( mat(x) ) == 0 % Want blocks to sum to 0
%         sum( mat(abs(x) ) ) <= 2
%         A*x == 0
%         sgn*(c*x) <= -sgn*1e-3 % quickly find certificate
%   cvx_end
%   if ~all( x== 0 )
%       resolving = false;
%       break;
%   end
%   toc
%   fprintf('Finished sgn = %d case\n', sgn );
% end
% if ~resolving
%     disp('Found a certificate that this is NOT a resolving set');
%     X   = mat(x)
% else
%     disp('Only solution is the zero solution, so this IS a resolving set');
%     X   = mat(x);
% end
% 

