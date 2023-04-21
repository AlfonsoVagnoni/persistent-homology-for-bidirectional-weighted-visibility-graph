function [VG,weight_input] = fast_BWVG(signal_var,coords_var,weight,B_period)

%%  INPUT CHECK
    VG=[];
    if isvector(signal_var)==0 || isvector(coords_var)==0
        disp('Error size: series and coordinates must be vectors')
        return;
    end
    if length(signal_var)~=length(coords_var)
        disp('Error size: series and coordinates must have the same length')
        return;
    end
    if iscolumn(signal_var)==0
        signal_var=signal_var';
    end
    if iscolumn(coords_var)==0
        coords_var=coords_var';
    end
    if exist('weight','var')==0
        weight='u';
    elseif ischar(weight)==0 || (strcmp(weight,'u')+strcmp(weight,'w'))==0
        disp('Error weight: ''u''=unweighted; ''w''=weighted ')
        return;
    end    
    if exist('B_period','var')==0
        B_period=0;
    elseif isscalar(B_period)==0 || sum(ismember([0,1],B_period))==0
        disp('Error boundary periodicity: ''0''=no periodicity; ''1''=periodicity ')
        return;
    end  
    if B_period==1 && length(unique(diff(coords_var)))~=1
        disp('Error: the vector with coordinates must be uniformly-spaced. Use NaN in the signal to account for non-uniform sampling.')
        disp('E.g.: signal_var=[1,NaN,4.2,NaN,NaN,8.4]; coords_var=[1,2,3,4,5,6]')
        return;
    end
    
    
%% PRE-PROCESSING
    
    if B_period==0
        TS2map=signal_var;
        TT2map=coords_var;
    elseif B_period==1
        %re-arrange the signal to account for boundary periodicity
        Nt=size(signal_var,1);
        [~,peaks_signal_ind]=max(signal_var);
        [~,sorted_indx]=sort([(peaks_signal_ind:Nt)';(1:peaks_signal_ind-1)'],'ascend'); %keep the indexing order
        TS2map=[signal_var(peaks_signal_ind:end);signal_var(1:peaks_signal_ind)]; % the peaks_signal_ind node is put both at the beginning and at the end of the new signal
        TT2map=(1:length(TS2map))';
    end
    
    weight=find(ismember({'u','w'},weight)); %weight=1 or 2
    clear coords_var
%% RUNNING
    N=size(TS2map,1);
    B=cell(N,weight);
    B=NVG_alg(TS2map,1,N,B,TT2map,weight);
    
%% CONVERTING INTO SPARSE MATRIX

    Adj_tmp=cell(N,1);
    for ii=1:N
        Adj_tmp{ii}=ones(length(B{ii,1}),1)*ii;
    end
    target1=vertcat(Adj_tmp{:});
    clear A_tmp
    source1=vertcat(B{:,1});
    target=[target1;source1];
    source=[source1;target1];
    weight_input=zeros(size(target));
    for iii=1:length(target)
        if signal_var(target(iii))==signal_var(source(iii))
            weight_input(iii)=1e-25;
        elseif signal_var(target(iii))>signal_var(source(iii))
            weight_input(iii)=abs(atan((signal_var(target(iii))-signal_var(source(iii)))/(target(iii)-source(iii))));
        elseif signal_var(target(iii))<signal_var(source(iii))
            weight_input(iii)=abs(atan((target(iii)-source(iii))/(signal_var(target(iii))-signal_var(source(iii)))));
        end
    end

    if ~isa(source,'double')
        source=double(source); 
    end
    VG=sparse(source,target,weight_input,N,N);

    %% restore the original indexing and remove duplicate links
    if B_period==1
        % remove the link between the maximum with the auxiliary node (they are the same node)
        VG(1,end)=0;
        VG(end,1)=0;
        %merge the links from the maximum (1st node) and the auxiliary (last node)
        VG(:,1)=or(VG(:,1),VG(:,end));
        VG(1,:)=VG(:,1)'; %=or(A(1,:),A(end,:));
        % remore the auxiliary node
        VG(:,end)=[];
        VG(end,:)=[];
        % restore the initial indexing
        VG=VG(sorted_indx,sorted_indx);
    end  

end

