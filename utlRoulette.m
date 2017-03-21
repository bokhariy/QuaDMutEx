%--------------------------------------------------------------------------------------
%	A function returns genes randomly seleceted based on thier probabilitis
%--------------------------------------------------------------------------------------
function [selectedIdx,selectedIDs] = utlRoulette(weights,cnt,IDs,selectedIdx);
	if (nargin<3) %n number/  ar arguments / input
 		IDs=[1:length(weights)];
    end
    if (nargin<4)
       	selectedIdx=[];
    end
    cweights=weights; %cumsum is already done outside, weights are actually cumsum(weights)
	composit=length(weights);
	k=0;
	while (k<cnt)
		r = rand(1);   % 31 random number between 0 and 1
		bingo=min(find(cweights>=r));
		if isempty(bingo)
			bingo=composit;
		end
		present = find(bingo==selectedIdx);
		if isempty(present)
			selectedIdx= [selectedIdx bingo];  % store the location of the coverage_weight which is the same as location of the gene set
			k=k+1;
		end	
	end
	selectedIDs= IDs(selectedIdx);





