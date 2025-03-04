function parset = genparset ( chain );
% Generates a 2D matrix parset from 3D array chain

% Determine how many elements in chain
[ T , d , N ] = size ( chain ); 

% Initalize parset
parset = [];

% If save in memory -> No -- parset is empty
if (T == 0),
    % Do nothing
else
    % parset derived from all chain
    for n = 1 : N,
        parset = [ parset; chain(:,:,n) (1:T)'];
    end
    parset = sortrows(parset,[d+1]); parset = parset(:,1:d);
end