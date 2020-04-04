function res = solve_quasi_random(this)
% solve_quasi_random works with quasi-random sampling
%
global switching_semantics;
seed = this.solver_options.quasi_rand_seed;
num_samples =  this.solver_options.num_quasi_rand_samples;

if ~strcmp(this.display,'off')
    fprintf('\n\n++++++++++++++++++++++++++++++++++++++++++++\nRunning %g quasi-random samples with seed %g\n', num_samples, seed);
    this.display_status_header();
end

BrQ = this.BrSet.copy();
BrQ.ResetParamSet();
BrQ.SetParamRanges(this.params, [this.lb this.ub])
BrQ.QuasiRandomSample(num_samples, seed);
X0 = BrQ.GetParam(this.params);

switching_semantics = 'semantic1';
if strcmp(switching_semantics, 'semantic1')
    res1 = this.FevalInit(X0);
    [fval1,j1] = sort(res1.fval(1,:));
    X1 = X0(:,j1);
end

switching_semantics = 'semantic2';
if strcmp(switching_semantics, 'semantic2')
    res2 = this.FevalInit(X0);
    [fval2,j2] = sort(res2.fval(1,:));
    X2 = X0(:,j2);
end

difX1 = (X1(:,end)-X1(:,1));
normX1 = sqrt(sum(difX1(:,1).^2));
difX2 = (X2(:,end)-X2(:,1));
normX2 = sqrt(sum(difX2(:,1).^2));

difF1 = (fval1(1,end)-fval1(1,1));
difF2 = (fval2(1,end)-fval2(1,1));

dis1 = difF1/normX1;
dis2 = difF2/normX2;
done = 0;
if dis1 > dis2
    switching_semantics = 'semantic1';
    res = res1;
    done = 1;
    
elseif dis1 < dis2
    switching_semantics = 'semantic2';
    res = res2;
    done = 1;
else
    for k = 1:length(fval1(1,:))
        
        if fval1(1,k) < fval2(1,k)
            
            switching_semantics = 'semantic1';
            res = res1;
            done = 1;
            break;
            
        elseif fval2(1,k) < fval1(1,k)
            
            switching_semantics = 'semantic2';
            res = res2;
            done = 1;
            break;
            
        end
        
    end
end

if done == 1
    
else
    switching_semantics = 'semantic1';
    res = res1;
    
end

this.add_res(res);
end


