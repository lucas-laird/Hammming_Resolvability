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

%% Generate a smaller resolving subset
n = 500;
sizes = zeros(1,n);
sets = cell(1,n);
for j = 1:n
    fprintf('Starting set %i \n',j)
    tic
    is_resolving = false;
    shuffled_order = randperm(82);
    integer_list = [];
    curr_size = 0;
    counter = 1;
    A_rank = 0;
    temp = [];
    X = false;
    while ~is_resolving
        temp = [integer_list,shuffled_order(counter)];
        counter = counter+1;
        A = A_82(temp,:);
        is_valid = check_valid(A,A_rank,X);
        if is_valid
            A_rank = rank(A);
            integer_list = temp;
            curr_size = curr_size+1;
            %fprintf('\nSet size = %i \n',curr_size);
            [is_resolving,X] = ILP_resolve(k,a,A);
        end
    end

    new_set = cell(curr_size,1);
    for i = 1:curr_size
        s = str{integer_list(i)};
        new_set{i} = s;
    end
    sets{j} = new_set;
    sizes(j) = curr_size;
    toc
end