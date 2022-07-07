N = length(W(1,:));   % Number of column vectors
a = 1/N;

B = zeros(3,3);
for i=1:N,
    B = B + a*(W(:,i)*V(:,i)');
end

S   = B+B';
rho = trace(B);

Z = zeros(3,1);
for i=1:N,
    Z = Z + a*cross(W(:,i),V(:,i));
end

K = [ S-rho*eye(3) Z
      Z'           rho ];

[E,D] = eig(K);          % Find the eigenvector for largest eigenvalue of K
[a,b] = max(diag(D));

qq      = E(:,b);         % Compute the unit quaternion
qq      = qq/norm(qq);
q       = [qq(4) 
          -qq(1:3)];

