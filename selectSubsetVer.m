%-----------------------------------------------------------------------------------------------------------------------
%	A function returns randomly seleceted genes based on number of genes in the solution set and selection probabilitis 
%-----------------------------------------------------------------------------------------------------------------------
function [selectIdx,selectedIDs ] = selectSubsetVer(scores, randy, weights, IDs, scr_IDs, scr_weights );
n=length(scores); % QP solution from Gurobi

if (nargin<6)
	scr_weights=ones(1,n)./n;
end

if (nargin<5)
	scr_IDs=[1:n];
end

if (nargin<4) %n number/  ar arguments / input
	IDs=[1:length(weights)];
end

cnt=1; % how many genes to be replaced in case number of once in the sol. > number of zeros

pre1=find(scores>0) ; % QP genes >0
pre0=find(scores==0); % QP genes ==0

if not (isempty(randy))
	current_genes_number_sol=randy(pre1);
else
	current_genes_number_sol=[];
end

%if 1's are <50% of current solution - replace all zeros with new genes
if (length(pre1)<length(pre0))
	selectIdx = current_genes_number_sol;
	[index selectedIDs]= utlRoulette(weights,length(pre0),IDs,selectIdx);
	selectIdx = unique([selectIdx index]);


else %otherwise, change one of the 1's to zero, and replace all zeros with new genes, 
	sub_length = cnt;   %length of 1s that need to be replaced randomly
	randpre1 = randperm(length(pre1)) ; 
	subpre1 = randpre1(1:sub_length);

	toBeZero=pre1(subpre1);
	scores(toBeZero)=0;
	pre1(subpre1)=[];
	current_genes_number_sol=randy(pre1);  %gene ID numbers

	pre0=find(scores==0);
	selectIdx = current_genes_number_sol;
	[index selectedIDs]= utlRoulette(weights,length(pre0),IDs,selectIdx);

	selectIdx = unique([selectIdx index]);


	current_genes_number_sol;
	randy;
	selectIdx;
end
