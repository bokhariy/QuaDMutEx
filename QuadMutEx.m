function [objective_function,x, solutionSetGeneNames] = QuadMutEx(G,GenesNames,t,n,c,k,GenesSelProb)
%
%Yahya Bokhari  3/1/2017
%
%
% matlab dirctory should be where Gurobi sofrware saved (ex: /opt/gurobi600/linux64/matlab/)
%
%
% Describtion of QuadMutEx_FromFiles arguments:
% 1- G: 
% 	mutation matrix in a sparse format. Patients have a serial number from 1 to m. Genes have serial number from 1 to n.
% 2- GenesNames:
%	vector of genes orderd as in the G sparse matrix(from 1 to n).
% 3- t:
%	Number of iterations desired.
% 4- n:
% 	Number of genes to be included in the submatrix G' each iteration.
% 5- c:
% 	Parameter c (As c gets bigger as the number of genes in solution shrinks).
% 6- k:
%	Parameter k (High value of k panalize for overlap).
% 7- GenesSelProb: (optional)
%	vector of probabilities of each gene to be chosen randomly.The probabilities should be in one column orderd as 
%	created serial number in the G_file(from 1 to n). The sum of the 
%	colmn should be equal 1.
%	
%
if (nargin<6) 
 		disp('Usage: QuadMutEx(G,GenesNames,t,n,c,k,GenesSelProb)')
 		disp(' G: Mutatation sparse matrix.') 
 		disp(' GenesNames: a vector contains list of genes orderd as in the G sparse matrix..')
 		disp(' t: Number of iterations.')
 		disp(' n: Number of genes to be included in the submatrix G''')
 		disp(' c: Parameter c.')
 		disp(' k: Parameter k.')
 		disp(' GenesSelProb: a vector of Genes selection probabilities.')
 		disp(' (Type: help QuadMutEx OR see README.txt for more details).')
 	return;
 end 



mySparseMatrix=G;
[N,M]= size(mySparseMatrix);

genes = GenesNames;

if (nargin<7)
	selection_prob=[1:length(genes)]./length(genes)
else
	selection_prob=GenesSelProb;
end

best_obj_f=10000000;  % we want to minimize
X = zeros(1,n);  % n is the number of genes to be tried each iteration.
loc_cand=[];
Objective_function=0;
prev_objective_function=Objective_function ;

best_x=zeros(1,n);
prev_x=X;
best_names='';
prev_names='';
randy=[];
best_randy=[];
prev_randy=[];
prev_BEST_G=[];

 
for j = 1 :t

	[randy , nms] =selectSubsetVer(X, randy, selection_prob);

	names=genes(randy);
	G=mySparseMatrix(:,randy);  %  new gene names and Gsubmatrix  

	[X,O]= Gurobi_optimization(G,names,c,k);  %Gurobi part

	Objective_function= N+O;
	Objective_function;
	best_obj_f;
	if (Objective_function<best_obj_f)
        %store to be returned
		best_obj_f=Objective_function;
		best_x=X;
		best_names=names;
		BEST_G=G;
		best_randy=randy;

        %always accept
		prev_objective_function=Objective_function;
		prev_x=X;
		prev_BEST_G=G;
		prev_names=names;

	end

	power =   Objective_function  - prev_objective_function;  
	T=2;
	(power/T)*-1;
	luck = min(1,exp(-(power/T)));
	coin= rand(1);

	if (coin <= luck)	%accept
		prev_objective_function=Objective_function;
		prev_x=X;
		prev_names=names;
		prev_BEST_G=G;
		prev_randy=randy;
	else % discard; to save last best sol and obj function
		X = prev_x; 
		names= prev_names;
		G=prev_BEST_G;
		randy=prev_randy;
		Objective_function=prev_objective_function;

	end  %END of MH LOOP


	loc_cand = find(best_x>0) ; %locations in THIS x solution = the real locations of the genes in the matrix if we started the genes that we started with were from the sorted list in term of covrage 


	
end  % end of iter

%%############  RESULTS  ###############


objective_function=best_obj_f;

solutionSetGeneNames=best_names([find(best_x ==1)]);

Gene_solution_set_location=[];
for z = 1 :length(solutionSetGeneNames)
	Gene_solution_set_location=[Gene_solution_set_location; find(strcmp(genes,solutionSetGeneNames(z))==1)];
end 
x=Gene_solution_set_location;
