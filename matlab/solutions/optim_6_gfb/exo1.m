randu = randn(N/2,N/2,2); % a random vector field
b = 2;
for i=1:4
    LLs_u = Li{i}( LiS{i}( randu ) );
    % relative error should be very small
    norm( abs( LLs_u(:) - b*randu(:) ) )/norm( randu(:) );
end
