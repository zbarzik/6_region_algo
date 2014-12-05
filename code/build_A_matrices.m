function [dat0] = build_A_matrices(bactData, reads, AlgoConfig)


if ~AlgoConfig.barcoded_regions
    dat0 = [];
    return
end

pe = AlgoConfig.pe;
% p_cut = AlgoConfig.p_cut;
nMM_cut = AlgoConfig.nMM_cut;
% nMM_fix = AlgoConfig.nMM_fix;

rL = AlgoConfig.readLen;
nR = length(bactData.kmers);
nB = size(bactData.indInSeqs,1);


for rr = 1:nR
    
    % Filter low abundance reads
    if AlgoConfig.filter_reads == 1
        tot_count = sum(reads{rr}.uniqueReads_count);
        filter_thresh = max(AlgoConfig.min_read_freq*tot_count,AlgoConfig.min_read_count);
        remove_ind = reads{rr}.uniqueReads_count < filter_thresh;
        reads{rr}.uniqueReads(remove_ind,:) = [];
        reads{rr}.uniqueReads_count(remove_ind) = [];
    end
    
    
    % Build M matrix
    if AlgoConfig.verbose
        disp(['Building matrix M for region ' num2str(rr) ' out of ' num2str(nR)])
    end
    bamp_in_reg = find(bactData.indInSeqs(:,rr)>0);
    
    if strcmp(AlgoConfig.read_type,'PE')
        Kd = size(bactData.kmers{rr},1);
        dat0.M{rr} = sparse(bactData.indInSeqs(bamp_in_reg,rr),bamp_in_reg,ones(1,length(bamp_in_reg)),Kd,nB);
        dat0.kmers{rr} = bactData.kmers{rr};
    elseif strcmp(AlgoConfig.read_type,'SE')
        % Fwd
        [values_fwd, ~, indInUni_fwd] = unique(bactData.kmers{rr}(:,1:rL),'rows');
        Kd_fwd = length(values_fwd);
        indInValue_fwd = zeros(nB,1);
        indInValue_fwd(bamp_in_reg) = indInUni_fwd(bactData.indInSeqs(bamp_in_reg,rr));
        Ad_fwd = sparse(indInValue_fwd(bamp_in_reg),bamp_in_reg,ones(1,length(bamp_in_reg)),Kd_fwd,nB);
        
        % Rvs
        [values_rvs, ~, indInUni_rvs] = unique(bactData.kmers{rr}(:,rL+1:end),'rows');
        Kd = length(values_rvs);
        indInValue_rvs = zeros(nB,1);
        indInValue_rvs(bamp_in_reg) = indInUni_rvs(bactData.indInSeqs(bamp_in_reg,rr));
        Ad_rvs = sparse(indInValue_rvs(bamp_in_reg),bamp_in_reg,ones(1,length(bamp_in_reg)),Kd,nB);
        
        % Write the matrix
        dat0.M{rr} = [Ad_fwd;Ad_rvs];
        dat0.kmers{rr} = [values_fwd; values_rvs];
    end
    sumMPerBactPerRegion = sum(dat0.M{rr},1);
    dat0.M{rr} = bsxfun(@rdivide,dat0.M{rr},(sumMPerBactPerRegion+eps));
    
    
    % Count errors between reads and kmers
    if AlgoConfig.verbose
        disp(['Building matrix A for region ' num2str(rr) ' out of ' num2str(nR)])
    end
    nY = length(reads{rr}.uniqueReads_count);
    cell_A_i = cell(1,nY);
    cell_A_j = cell(1,nY);
    cell_A_s = cell(1,nY);
    non_zero_counts = zeros(1,nY);

    nK = size(dat0.M{rr},1);
    debug_dist_mat = zeros(nY,size(dat0.kmers{rr},1));
    for yy = 1:nY
        if AlgoConfig.verbose && yy == 1000*fix(yy/1000)
            disp(['Mapped ' num2str(yy) ' reads out of ' num2str(nY)])
        end
        
        distvec = sum(bsxfun(@ne,dat0.kmers{rr},reads{rr}.uniqueReads(yy,:)),2);
        debug_dist_mat(yy,:) = distvec;
        min_dist_2db = min(distvec);
        
        %         % Correct the database
        %         dbfix_thresh = AlgoConfig.min_dbfix_freq*tot_count;
        %         if reads{rr}.uniqueReads_count(yy) > dbfix_thresh && min_dist_2db>0 && min_dist_2db<=nMM_fix
        %             kmer_2fix_ind = find(distvec == min_dist_2db);
        %             yy
        %             disp('OK-fix the DB')
        %         end

    
        if min_dist_2db <= nMM_cut
            %             Pr_read_given_kmer = ((1-pe).^(rL-distvec)).*((pe/3).^distvec);
            %             Pr_read_given_kmer(distvec>nMM_cut) = 0;
            %             %         Pr_read_given_kmer(Pr_read_given_kmer<p_cut) = 0;
            %             Pr_read_given_j = (Pr_read_given_kmer'*dat0.M{rr})./(sumMPerBactPerRegion+eps);
            
            mapped_kmer = distvec<=nMM_cut;
            tmp_Pr_read_given_kmer = ((1-pe).^(rL-distvec(mapped_kmer))).*((pe/3).^distvec(mapped_kmer));
            Pr_read_given_kmer = sparse(find(mapped_kmer),ones(sum(mapped_kmer),1),tmp_Pr_read_given_kmer,nK,1);
%             Pr_read_given_j = (Pr_read_given_kmer'*dat0.M{rr})./(sumMPerBactPerRegion+eps);
            Pr_read_given_j = Pr_read_given_kmer'*dat0.M{rr};
            
            mapped_bact = Pr_read_given_j>0;
            cell_A_i{yy} = yy*ones(1,sum(mapped_bact));
            cell_A_j{yy} = find(mapped_bact);
            cell_A_s{yy} = full(Pr_read_given_j(mapped_bact));
            non_zero_counts(yy) = sum(mapped_bact);
        end
        
    end
    
    
    % Build one sparse matrix
    i_all = [cell_A_i{:}];
    j_all = [cell_A_j{:}];
    s_all = [cell_A_s{:}];
    dat0.A{rr} = sparse(i_all,j_all,s_all,nY,nB);
    dat0.reads{rr} = reads{rr}.uniqueReads;
    dat0.F{rr} = reads{rr}.uniqueReads_count;

    
    % DEBUG
    if 0
        max(dat0.F{rr}(sum(dat0.A{rr}>0,2)==0))
        
        tmp = find(dat0.M{rr}>0);
        kmersind = mod(tmp,size(dat0.M{rr},1));
        kmersind(kmersind==0) = size(dat0.M{rr},1);
        [min_dist,I] = min(debug_dist_mat,[],1);
        best_match_read_ind = I(kmersind);
        best_match_read_dist = min_dist(kmersind);
        best_match_read_count = dat0.F{rr}(best_match_read_ind);
        
        read_count_difference = cell(1,nY);
        read_index_difference = cell(1,nY);
        for yy = 1:nY
            yy
            tmp_ind = find(dat0.A{rr}(yy,:)>0);
            xx_ind_mat = zeros(0,2);
            for xx = 1:length(tmp_ind)
                xx_ind = find(dat0.A{rr}(:,tmp_ind(xx))>(1-pe).^rL/2-eps);
                if length(xx_ind)>1
                    xx_ind_mat = [xx_ind_mat;xx_ind'];
                end
            end
            read_index_difference{yy} = unique(xx_ind_mat,'rows');
            
            read_count_difference{yy} = [dat0.F{rr}(read_index_difference{yy}(:,1)) dat0.F{rr}(read_index_difference{yy}(:,2))];
        end
    end

end
dat0.is_perfect_match = bactData.is_perfect_match;
dat0.is_amplified = bactData.indInSeqs>0;

