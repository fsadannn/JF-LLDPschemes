function solver_output = llfinalize(solver, sol,...
                                     outfun, outargs,...
                                     printstats, statvect,...
                                     nout, tout, yout,...
                                     haveeventfun, teout, yeout, ieout,...
                                     interp_data)
%ODEFINALIZE Helper function called by LLDP solvers at the end of integration.
%

if ~isempty(outfun)
  feval(outfun,[],[],'done',outargs{:});
end

fullstats = (length(statvect) > 3);  % faster than 'switch' or 'ismember'

stats = struct('nsteps',statvect(1),'nfailed',statvect(2),'nfevals',statvect(3)); 
if fullstats
  stats.npds     = statvect(4);
  stats.ndecomps = statvect(5);
  stats.nsolves  = statvect(6);  
else 
  statvect(4:6) = 0;   % Backwards compatibility
end  

if printstats
  fprintf(getString(message('llint:llfinalize:LogSuccessfulSteps', sprintf('%g',stats.nsteps))));
  fprintf(getString(message('llint:llfinalize:LogFailedAttempts', sprintf('%g',stats.nfailed))));
  fprintf(getString(message('llint:llfinalize:LogFunctionEvaluations', sprintf('%g',stats.nfevals))));
  if fullstats
    fprintf(getString(message('llint:llfinalize:LogPartialDerivatives', sprintf('%g',stats.npds))));
    fprintf(getString(message('llint:llfinalize:LogLUDecompositions', sprintf('%g',stats.ndecomps))));
    fprintf(getString(message('llint:llfinalize:LogSolutionsOfLinearSystems', sprintf('%g',stats.nsolves))));
  end
end

solver_output = {};

if (nout > 0) % produce output
  if isempty(sol) % output [t,y,...]
    solver_output{1} = tout(1:nout).';
    solver_output{2} = yout(:,1:nout).';
    if haveeventfun
      solver_output{3} = teout.';
      solver_output{4} = yeout.';
      solver_output{5} = ieout.';
    end
    solver_output{end+1} = statvect(:);  % Column vector
  else % output sol  
    % Add remaining fields
    sol.x = tout(1:nout);
    sol.y = yout(:,1:nout);
    if haveeventfun
      sol.xe = teout;
      sol.ye = yeout;
      sol.ie = ieout;
    end
    sol.stats = stats;
    switch solver
     case 'LLDP_Kphi1'
      [f3d,idxNonNegative] = deal(interp_data{:});
%       sol.idata.f3d = f3d(:,:,1:nout);      
      sol.idata.idxNonNegative = idxNonNegative;
     case 'LLDP1'
      [f3d,idxNonNegative] = deal(interp_data{:});
%       sol.idata.f3d = f3d(:,:,1:nout);      
      sol.idata.idxNonNegative = idxNonNegative;
      case 'LLDP2'
      [f3d,idxNonNegative] = deal(interp_data{:});
%       sol.idata.f3d = f3d(:,:,1:nout);      
      sol.idata.idxNonNegative = idxNonNegative;
    case 'LLDP_exact'
      [f3d,idxNonNegative] = deal(interp_data{:});
%       sol.idata.f3d = f3d(:,:,1:nout);      
      sol.idata.idxNonNegative = idxNonNegative;
    case 'ode15sk'      
      [kvec,dif3d,idxNonNegative] = deal(interp_data{:});
      sol.idata.kvec = kvec(1:nout);
      maxkvec = max(sol.idata.kvec);
      %sol.idata.dif3d = dif3d(:,1:maxkvec+2,1:nout);
      sol.idata.idxNonNegative = idxNonNegative;
    case 'ode15sk-fj'      
      [kvec,dif3d,idxNonNegative] = deal(interp_data{:});
      sol.idata.kvec = kvec(1:nout);
      maxkvec = max(sol.idata.kvec);
      %sol.idata.dif3d = dif3d(:,1:maxkvec+2,1:nout);
      sol.idata.idxNonNegative = idxNonNegative;
     otherwise
      error(message('llint:llfinalize:UnrecognizedSolver', solver));
    end  
    solver_output{1} = sol; 
  end
end    
