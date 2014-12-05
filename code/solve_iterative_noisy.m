function [output_freq,l2_freq] = solve_iterative_noisy(dat0,Config)

if strcmp(Config.read_type,'PE')
    error('Make sure PE is supported properly (filter????)')
end

nR = length(dat0.M);
nB_all = size(dat0.M{1},2);


% Filter columns (bacterias)
% We keep bacteria that are present in all the 0MM regions, any of 2MM regions
if Config.verbose
    disp('Filter out columns (bacteria) with no reads')
end
missing_0MM_region = false(nR,nB_all);
A_mat = dat0.A;
sumAPerBactPerRegion = zeros(nR,nB_all);

% This is the minimal probabilty mass that should be present in reads 
% (it could be larger than this for a cirtain bacteria if some other kmers are mapping to kmer of given bacteria)
% total_pr_thresh = (1-Config.pe)^Config.readLen;
total_pr_thresh = ((1-Config.pe).^(Config.readLen-Config.nMM_cut)).*((Config.pe/3).^Config.nMM_cut); 
perfect_match_pr_SE = (1-Config.pe)^Config.readLen/2;
max_error_pr_SE = total_pr_thresh/2;
for i = 1:nR
    % Calculate the total observed probabilty per bacteria in region
    sumAPerBactPerRegion(i,:) = sum(A_mat{i},1);
    
    if Config.do_filter == 1
        if 0
            % Mark bacteria that dont have both FWD and RVS in 0MM region as absent
            missing_0MM_region(i,:) = dat0.is_perfect_match(:,i)' & sumAPerBactPerRegion(i,:) < total_pr_thresh;
            
            % Mark bacteria that dont have both FWD and RVS as not present in 2MM region
            SE_ind = ~dat0.is_perfect_match(:,i)' & sumAPerBactPerRegion(i,:) < total_pr_thresh;
            A_mat{i}(:,SE_ind) = 0;

            % Count frequency in reads
            sumAPerBactPerRegion(i,:) = sum(A_mat{i},1);
        elseif 0
            missing_0MM_region(i,:) = dat0.is_amplified(:,i)' & sumAPerBactPerRegion(i,:) < total_pr_thresh;
        else
            num_perfect_matches = sum(A_mat{i} > perfect_match_pr_SE-0.1*max_error_pr_SE-eps,1);
            warning('TAKE PROPER CARE OF NOT AMPLIFIED REGIONS')
            missing_0MM_region(i,:) = dat0.is_amplified(:,i)' & num_perfect_matches<2; 
        end
    
    end
end

% Remove zero rows of A (even if there is a read there
F_vec = dat0.F;
DEBUG_F_vec = zeros(sum(cellfun(@(x) size(x,1),A_mat)),1); k= 0;
for i = 1:nR
    % If we have a row of zeros in A and a non zero read - we assume that this
    % is due to the fact that pair of this read was not observed due to low
    % coverage and we erased the single one from A matrix because if the region
    % was amplified we expect to get 2 reads (for SE) - so now we just erase
    % this row - assuming that the double read (PE) was not observed
    remove_ind = find(sum(A_mat{i},2)==0);
    F_vec{i}(remove_ind) = [];
    A_mat{i}(remove_ind,:) = [];    

    DEBUG_F_vec(k+1:k+size(A_mat{i},1)) = F_vec{i};
    k = k+ size(A_mat{i},1);
end


% Keep only bacterias present in all the regions 0MM regions and at least 1 region
keep_col = find(sum(missing_0MM_region,1)==0 & sum(sumAPerBactPerRegion,1) > 0);
nB = length(keep_col);

% Normalize frequency counts
if Config.verbose
    disp('Normalize frequency counts')
end
numReadsPerRegion = cellfun(@sum,F_vec);
freq_vec = {};
for i = 1:nR
    
    % Normalize the read counts
    if strcmp(Config.mixture_type,'Multiplex')
        freq_vec{i} = F_vec{i}/sum(numReadsPerRegion);
    elseif strcmp(Config.mixture_type,'RegionByRegion')
        freq_vec{i} = F_vec{i}/numReadsPerRegion(i)/nR;
    end
    
end



% Write one matrix for all the measurements
if Config.verbose
    disp('Build matrix A_L2')
end
A_L2 = zeros(sum(cellfun(@(x) size(x,1),A_mat)),nB);
y_vec = zeros(size(A_L2,1),1);
DEBUG_F_vec = zeros(size(A_L2,1),1);
Rg_vec = zeros(size(A_L2,1),1);
k = 0;
for i = 1:nR
    A_L2(k+1:k+size(A_mat{i},1),1:nB) = A_mat{i}(:,keep_col);
    y_vec(k+1:k+size(A_mat{i},1)) = freq_vec{i};
    DEBUG_F_vec(k+1:k+size(A_mat{i},1)) = F_vec{i};
    Rg_vec(k+1:k+size(A_mat{i},1)) = i;
    k = k+ size(A_mat{i},1);
end
A_L2 = sparse(A_L2);
A_L2 = A_L2/nR;

% Find unique columns
if Config.do_unique == 1
    if Config.verbose
        disp('Making columns of A unique...')
    end
    [~,I,J] = unique(A_L2','rows');
    %     my_A = bsxfun(@rdivide,A_L2,sqrt(sum(A_L2.^2,1)));
    %     my_I = zeros(1,nB);
    %     tic
    %     for ii = nB:-1:1
    %         if ii == 10000*fix(ii/10000)
    %             disp(['Checked uniqueness for ' num2str(ii) ' out of ' num2str(nB) ' in ' num2str(toc) ' sec'])
    %         end
    %         tmp = my_A(:,ii)'*my_A(:,ii:end);
    %         min_ind_uni = find(tmp>=(1-sqrt(eps)),1,'last');
    %         if min_ind_uni == 1
    %             my_I(ii) = 1;
    %         end
    %     end

    
    % It is important to mark all the bacteria as present at the end since due
    % to filter could be mismatch
    A_L2 = A_L2(:,I);
    keep_col = keep_col(I);
    nB = length(keep_col);
end
    
% Filter out bacterias "included" in other bacteria
if Config.filter_included_bacteria == 1
    if Config.verbose
        disp('Removing included bacterias...')
    end
    A_L2_TF = double((A_L2 > 0));
    num_non_zeros = full(sum(A_L2_TF,1));
    n_ii_jj_mat = full(A_L2_TF'*A_L2_TF);
    remove_included = false(1,nB);
    for ii = 1:nB
        jj_ii_diff = num_non_zeros - num_non_zeros(ii);
        this_row_min_n = num_non_zeros - (jj_ii_diff>0).*jj_ii_diff;
        if sum((n_ii_jj_mat(ii,:)-this_row_min_n) == 0 & jj_ii_diff > 0)
            remove_included(ii) = true;
        end
        %         for jj = ii+1:nB
        %             if n_ii_jj_mat(ii,jj) == this_row_min_n(jj)
        %                 if num_non_zeros(ii) < num_non_zeros(jj)
        %                     remove_included(ii) = true;
        %                 elseif num_non_zeros(ii) > num_non_zeros(jj)
        %                     remove_included(jj) = true;
        %                 end
        %             end
        %         end
    end
    A_L2(:,remove_included) = [];
    keep_col(remove_included) = [];
    if Config.verbose
        disp(['Removed ' num2str(sum(remove_included)) ' out of ' num2str(length(remove_included))])
    end
end

num_reads_per_bact = sum(A_L2>0,1);
non_even_count = sum(num_reads_per_bact~=2*fix(num_reads_per_bact/2));
if Config.verbose
    warning(['Found ' num2str(non_even_count) ' bacterias with non even number of reads mapped'])
    disp('Starting iterations...')
end

% Start the iterations
tic
tStart = tic;
% nB = length(keep_col);
% bact_freq = ones(nB,1)/nB;
bact_freq = A_L2'*y_vec;
bact_freq = bact_freq/sum(bact_freq);
for jj = 1:Config.numIter
    
    % work region by region
    %         % Init the factor
    %         bact_factor = zeros(1,nB);
    %
    %         for i = 1:nR
    %             % Estimate theta
    %             theta_i = A_mat{i}*bact_freq;
    %
    %             % Reweight the counts
    %             r_weighted = freq_vec{i}./theta_i;
    %
    %             % Update the the factor
    %             bact_factor = bact_factor + r_weighted'*A_mat{i};
    %         end
    
    % Estimate theta
    theta_i = A_L2*bact_freq;
    
    % Reweight the counts
    r_weighted = y_vec./theta_i;
    
    % Update the the factor
    bact_factor = r_weighted'*A_L2;

    
    % Calculate the error
    L1_error = abs(1-bact_factor)*bact_freq;
    L2_error = norm((1-bact_factor)'.*bact_freq);
    FACT_error = mean(abs(1 - bact_factor(bact_freq>1e-10)));

    % Update the frequency
    bact_freq = bact_freq.*bact_factor';
    
    % Check if rached tolerance
    if L1_error < Config.tol
        break
    end
    
    tTmp = toc;
    if tTmp > 60
        % Remove bacterias with zero frequency
        remove_ind = find(bact_freq < 1e-10);
        A_L2(:,remove_ind) = [];
        keep_col(remove_ind) = [];
        bact_freq(remove_ind) = [];
        
        % Print
        disp(['Iter:' num2str(jj) '. Error reduction of X (L1 norm): ' num2str(L1_error) '. Factor error: ' num2str(FACT_error)])
        tic
    end
    
    

end
if Config.verbose
    disp(['Total iterations time: ' num2str(toc(tStart))])
end

% Hard thresholding
bact_freq(bact_freq<1e-10) = 0;
bact_freq = bact_freq/sum(bact_freq);


% Sove LS for comparison
l2_freq = zeros(nB_all,1);
if Config.solve_L2 == 1
    %     F_mat = zeros(size(A_L2,1),nR);
    %     k = 0;
    %     for i = 1:nR
    %         F_mat(k+1:k+size(A_mat{i},1),nB+i) = freq_vec{i}*nR;
    %         k = k+ size(A_mat{i},1);
    %     end
    %
    %     tic
    %     x = solve_L2_fourth(A_L2, F_mat, 2);
    %     toc
    %
    %     l2_freq(keep_col) = x;
    l2_bact_freq = lsqnonneg(A_L2,y_vec);
    l2_freq(keep_col) = l2_bact_freq;
end


%
%
%
% % Dummy result
% bact_freq = At*y_vec;
% bact_freq = bact_freq/sum(bact_freq);
% 
% 

% Set the output vector
output_freq = zeros(nB_all,1);
output_freq(keep_col) = bact_freq;

