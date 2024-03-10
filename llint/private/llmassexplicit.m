function [odeFcn,odeFArgs,odeFxcn,odeFxArgs] = llmassexplicit( FcnHandlesUsed,FcnHandlesUsedFx,...
    massType,odeFcn,odeFxcn,odeArgs,massFcn,massM)  
%ODEMASSEXPLICIT  Helper function for handling the mass matrix
%   For explicit ll solvers -- incorporate the mass matrix into the ll
%   function.   
%
massL=[];
massU=[];
massp=[];
if FcnHandlesUsed
    switch massType
      case 1  % use LU factors of constant M 
        if issparse(massM) 
            %[massL,massU,massP,massQ,massR] = lu(massM);
            [massL,massU] = ilu(massM,struct('type','nofill'));
            %odeArgs = [{odeFcn,massL,massU,massP,massQ,massR},odeArgs]; 
            odeFArgs = [{odeFcn,massL,massU},odeArgs];
            odeFcn = @ExplicitSolverHandleMass1sparse;       
        else % M full
            [massL,massU,massp] = lu(massM,'vector');
            odeFArgs = [{odeFcn,massL,massU,massp},odeArgs];    
            odeFcn = @ExplicitSolverHandleMass1;
        end  
      case 2
        odeFArgs = [{odeFcn,massFcn},odeArgs];    
        odeFcn = @ExplicitSolverHandleMass2;
      otherwise % case {3,4}
        odeFArgs = [{odeFcn,massFcn},odeArgs];    
        odeFcn = @ExplicitSolverHandleMass34;
    end
else % ode-file:  F(t,y,'mass',p1,p2...)    
    if massType == 1   % use LU factors of constant M 
        if issparse(massM) 
            %[massL,massU,massP,massQ,massR] = lu(massM);
            [massL,massU] = ilu(Mt,struct('type','nofill'));
            %odeArgs = [{odeFcn,massL,massU,massP,massQ,massR},odeArgs];
            odeFArgs = [{odeFcn,massL,massU},odeArgs];
            odeFcn = @ExplicitSolverHandleMass1sparse;       
        else % M full
            [massL,massU,massp] = lu(massM,'vector');
            odeFArgs = [{odeFcn,massL,massU,massp},odeArgs];    
            odeFcn = @ExplicitSolverHandleMass1;
        end  
    else  
        odeFArgs = [{odeFcn},odeArgs];  
        odeFcn = @ExplicitSolverHandleMassOld;   
    end
end

if FcnHandlesUsedFx
    switch massType
      case 1  % use LU factors of constant M 
        if issparse(massM) 
            %[massL,massU,massP,massQ,massR] = lu(massM);
%             [massL,massU] = ilu(Mt,struct('type','nofill'));
            %odeArgs = [{odeFcn,massL,massU,massP,massQ,massR},odeArgs]; 
            odeFxArgs = [{odeFxcn,massL,massU},odeArgs];
            odeFxcn = @ExplicitSolverHandleMass1sparse;       
        else % M full
%             [massL,massU,massp] = lu(massM,'vector');
            odeFxArgs = [{odeFxcn,massL,massU,massp},odeArgs];    
            odeFxcn = @ExplicitSolverHandleMass11;
        end  
      case 2
        odeFxArgs = [{odeFxcn,massFcn},odeArgs];    
        odeFxcn = @ExplicitSolverHandleMass2;
      otherwise % case {3,4}
        odeFxArgs = [{odeFxcn,massFcn},odeArgs];    
        odeFxcn = @ExplicitSolverHandleMass34;
    end
else % ode-file:  F(t,y,'mass',p1,p2...)    
    if massType == 1   % use LU factors of constant M 
        if issparse(massM) 
            %[massL,massU,massP,massQ,massR] = lu(massM);
%             [massL,massU] = ilu(Mt,struct('type','nofill'));
            %odeArgs = [{odeFcn,massL,massU,massP,massQ,massR},odeArgs];
            odeFxArgs = [{odeFxcn,massL,massU},odeArgs];
            odeFxcn = @ExplicitSolverHandleMass1sparse;       
        else % M full
%             [massL,massU,massp] = lu(massM,'vector');
            odeFxArgs = [{odeFxcn,massL,massU,massp},odeArgs];    
            odeFxcn = @ExplicitSolverHandleMass1;
        end  
    else  
        odeFxArgs = [{odeFxcn},odeArgs];  
        odeFxcn = @ExplicitSolverHandleMassOld;   
    end
end

% --------------------------------------------------------------------------

function yp = ExplicitSolverHandleMass1(t,y,odeFcn,L,U,p,varargin)
  ode = feval(odeFcn,t,y,varargin{:});
  yp = U \ (L \ ode(p));
  
function yp = ExplicitSolverHandleMass11(t,y,odeFcn,L,U,p,varargin)
  ode = feval(odeFcn,t,y,varargin{:});
  yp = U \ (L \ ode(p,:));

% --------------------------------------------------------------------------

function yp = ExplicitSolverHandleMass1sparse(t,y,odeFcn,L,U,varargin)
  yp =  U \ (L  \ feval(odeFcn,t,y,varargin{:}));
 
% --------------------------------------------------------------------------
  
function yp = ExplicitSolverHandleMass2(t,y,odeFcn,massFcn,varargin)
  yp = feval(massFcn,t,varargin{:}) \ feval(odeFcn,t,y,varargin{:});
  
% --------------------------------------------------------------------------  

function yp = ExplicitSolverHandleMass34(t,y,odeFcn,massFcn,varargin)
  yp = feval(massFcn,t,y,varargin{:}) \ feval(odeFcn,t,y,varargin{:});

% --------------------------------------------------------------------------  
  
function yp = ExplicitSolverHandleMassOld(t,y,odeFcn,varargin)
  yp = feval(odeFcn,t,y,'mass',varargin{2:end}) \ ...
       feval(odeFcn,t,y,varargin{:});
  
% --------------------------------------------------------------------------  
