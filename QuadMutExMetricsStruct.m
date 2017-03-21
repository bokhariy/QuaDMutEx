function [met]=QuadMutExMetricsStruct(sol,mx,c,k)
%returns a structure met with quantitative evaluation of a solution vector (sol) given the mutation matrix (mx), and parameters c & k (see manuscript)
%    met.qObj=qObj; % quadratic objective function
%    met.lObj=lObj; % Dendrix score
%    met.covR=covR; % coverage
%    met.excessR=excessR; % excess coverage
%    met.mCnt=mCnt; % mutations in the selected genes
%    met.sCnt=sCnt; % samples #
%    met.gCnt=gCnt; % genes in the solution
%    met.covTot=covTot; % # of covered patients (samples)
%    met.covOvr=covOvr; % # of mutations in excess of 1 in a patient (i.e. coverage overlap as defined by Dendrix)

    solMx=full(mx(:,sol));
    [sCnt,gCnt]=size(solMx);
    covMx=sum(solMx,2);
    mCnt=sum(covMx);
    covTot=length(find(covMx>0));
    covR=covTot/sCnt;
    covOvr=mCnt-covTot;
    excessR=length(find(covMx>1))./covTot;
    lObj=covTot-covOvr;
    
    G=solMx;
    Q=transpose(G) * G;
	H=.5*(Q+transpose(Q));
	H=.5*(k+1)*H;
	f=(-((k+3)/2) *ones(1,sCnt))*G;
    f=f+c.*ones(size(f));
    x=ones(1,gCnt);
    qObj=sCnt+x*H*(x')+f*(x');
    
    met.qObj=qObj; % quadratic objective function
    met.lObj=lObj; % Dendrix score
    met.covR=covR; % coverage
    met.excessR=excessR; % excess coverage
    met.mCnt=mCnt; % mutations in the selected genes
    met.sCnt=sCnt; % samples #
    met.gCnt=gCnt; % genes in the solution
    met.covTot=covTot; % # of covered patients (samples)
    met.covOvr=covOvr; % # of mutations in excess of 1 in a patient (i.e. coverage overlap as defined by Dendrix)
    