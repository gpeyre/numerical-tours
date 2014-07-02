function [x,f,it,B] = perform_linprog(A,b,c,k,maxit,tol)

%LINPROG  This code uses the revised simplex method to solve the linear
%       programming problem: Minimize the cost c'x subject to 
%       equations Ax=b and nonnegativity x >=0:
%
%                f =  Min { c'x; Ax = b, x >= 0 },
%
%       You must define an m by n matrix A , a column vector b with 
%       m components, and a column vector c with n components.
%       You may define k to change the first k equations of Ax=b to
%       inequalities: thus [A(1:k,:)]x <= b(1:k).   
%
%       The output vector x gives the minimum cost, which is the output f.
%
%                [x,f,itn] = linprog(A,b,c,k,maxit,tol)
%
%       At most "maxit" iterations (default 10*length(b)) are applied
%       and the actual number of iterations is returned in "itn".
%
%       If the optimal solution is unbounded or the constraints are
%       inconsistent then a diagnostic is displayed.  
%       Bland's rule is used to resolve degeneracies.  In exact
%       arithmetic cycling is not possible.  But in real life!!!
%       If x has more than 20 components it is returned in the
%       sparse format. Note that if A has many zeros it is worth
%       passing it to linprog in sparse format. 
%
%       Although written for teaching purposes this routine has
%       successfully solved some problems with size(A) = [50,100000]!
%
%       Please report any difficulties to: idc@math.canterbury.ac.nz

%       New version                  (c) I.D.Coope, 1988, 1993


[m,n]=size(A); b=b(:); c=c(:); it=0;
if (length(c)~=n | length(b)~=m),error('wrong dimensions'); end
if (nargin<6), tol=1e-10; end
if (nargin<5), maxit=10*m; end
if (nargin<4), k=0; elseif isempty(k), k=0; end
D=sign(sign(b)+.5); if k, D(1:k)=ones(k,1); end
D = diag(D);                    % initial (inverse) basis matrix
A = [A D];                      % incorporate slack/artificial variables
B = n+1:n+m;                    % initial basis
N = 1:n;                        % non-basis
[bmin,j]=min(b(1:k));
if bmin<0,
   phase=1; xb=ones(m,1); s=[zeros(n+k,1);ones(m-k+1,1)];   % supercost
   N=[N,B(j)]; J=B; J(j)=[]; B(j)=n+m+1;
   a=b-sum(A(:,J)')'; A=[A a];
   D(:,j)= -a/a(j); D(j,j)=1/a(j);
elseif k==m,
   phase=2; xb=b; s=[c;zeros(m,1)];        % cost function
else                            % k==[] or bmin>=0
   phase=1; xb=abs(b); s=[zeros(n+k,1);ones(m-k,1)];   % supercost
end
while phase<3,
   df=-1; t=inf;
   yb= D'*s(B);                    % multipliers for Ax=b
   while (it < maxit)
      if isempty(N), break, end     % no freedom for minimization
      r = s(N) - [A(:,N)]'*yb;      % reduced costs
      [rmin,q] = min(r);            % determine new basic variable
      if rmin>=-tol*(norm(s(N),inf)+1), break, end % optimal!
      it=it+1;
      if df>=0                      % apply Bland's rule to avoid cycling
         if maxit==inf,
            disp(['LINPROG(',int2str(it),'): warning! degenerate vertex']);
         end
         J=find(r<0); Nq=min(N(J)); q=find(N==Nq); 
      end 
      d = D*A(:,N(q));
      I=find(d>tol); %    I=find(d>0);
      if isempty(I), disp('Solution is unbounded'); it=-it; break; end
      xbd=xb(I)./d(I); [r,p]=min(xbd); p=I(p);
      if df>=0,                     % apply Bland's rule to avoid cycling
         J=find(xbd==r); Bp=min(B(I(J))); p=find(B==Bp); 
      end 
      xb= xb - r*d; xb(p)=r;        % update x
      df=r*rmin;                    % change in f 
      v = D(p,:)/d(p);              % row vector
      yb= yb + v'*( s(N(q)) - d'*s(B) );
      d(p)=d(p)-1;
      D = D - d*v;                  % update inverse basis matrix
      t=B(p); B(p)=N(q);
      if t>n+k, N(q)=[]; else N(q)=t; end
   end                             % end of phase
   xb=xb+D*(b-A(:,B)*xb);          % iterative refinement
   I=find(xb<0);                   % must be due to rounding error
   if I, xb(I)=xb(I)-xb(I); end    % so correct
   if phase==2 | it<0, break; end; % B, xb,n,m,res=A(:,B)*xb-b
   if xb'*s(B)>tol,it=-it; disp('no feasible solution'); break;  end
   phase=phase+1;                  % re-initialise for Phase 2
   s=1e6*norm(c,'inf')*s; s(1:n)=c;% tol=tol*norm(s,inf);
end
x=sparse(n,1); x(B)=xb; x=x(1:n); if n<21, x=full(x); end
f=c'*x;
if it>=maxit, disp('too many iterations'); it=-it; end




