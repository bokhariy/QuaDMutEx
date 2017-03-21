%----------------------------------------------------------------------------------
%	A function uses Gurobi optimization softwar to solve the binary quadratic program.
%----------------------------------------------------------------------------------
function [x, objective_function] = Gurobi_optimization(G,names,c,k)
	Q=transpose(G) * G;
	H=.5*(Q+transpose(Q));
	H=.5*(k+1)*H;
	%fullMat=full(G);
	%[m,n] = size(fullMat);
	[m,n] = size(G);
	f=(-((k+3)/2) *ones(1,m))*G;
	model.Q=H;
	model.varnames=names;
	model.Obj=f+c.*ones(size(f));
	model.vtype='B';
	model.moddsense='min';
	model.sense='>';
	I=eye(n);
	Ix=sparse(I);
	model.A=Ix;
	zrs=zeros(1,n);
	model.rhs=zrs;

	%gurobi_write(model, 'qp.lp');
	%params.resultfile = 'mip1.lp';

	params.outputflag = 0;
    params.threads=1;
	results = gurobi(model,params);

	x = results.x;
	objective_function = results.objval;
	
