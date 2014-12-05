clear
simulate_reads = 0;

% Load the regions data
dbPath = '../Database';
dbFileName = 'GreenGenes_201305_unique_up_to_3_ambiguous_16S_amp6Regions_2mm';

regions_files = dir([dbPath '/' dbFileName '*.mat']);
nR = length(regions_files);
load([dbPath '/' regions_files(1).name],'Dictionary','Header_amp');
nB = length(Header_amp);
bactData.Header_Dictionary = Dictionary;
bactData.Header_amp = Header_amp;
bactData.kmers = cell(1,nR);
bactData.indInSeqs = zeros(nB,nR);
bactData.is_perfect_match = zeros(nB,nR);
for rr = 1:nR
    disp(['Loading bacterial DB for region ' num2str(rr) ' out of ' num2str(nR)])
    load([dbPath '/' regions_files(rr).name],'values','is_perfect_match','indInValue');
    bactData.kmers{rr} = values;
    bactData.is_perfect_match(:,rr) = is_perfect_match;
    bactData.indInSeqs(:,rr) = indInValue;
end


% *************************** LOAD EXPERIMENTAL READS ********************
sample_num = 20;

readsPath = '../Reads';
resDir = '../Results';

files = dir([readsPath '/sample_' num2str(sample_num) '_region*.mat']);
nR = length(files);
experimental_reads{nR} = struct('Suni', 0, 'freq', 0);
for rr = 1:nR
    tmp = load([readsPath '/sample_' num2str(sample_num) '_region_' num2str(rr) '_unireads.mat']);
    
    % Make sure the reads are unique
    [S,xi] = sortrows(tmp.readsuni);
    freq1 = tmp.frequni(xi);
    [Suni, ia] = unique(S,'rows', 'first');
    [~, ib] = unique(S,'rows', 'last');
    cumcount = [0;cumsum(freq1)];
    freq = cumcount(ib+1)-cumcount(ia);
    
    
    experimental_reads{rr}.uniqueReads = Suni;
    experimental_reads{rr}.uniqueReads_count = freq;
end



% *************************** RUN THE ALGORITHM ********************

% Config the algorithm
AlgoConfig.min_read_freq = 1e-4;
AlgoConfig.min_read_count = 2;

AlgoConfig.nMM_fix = 1; 
AlgoConfig.min_dbfix_freq = 1e-3;

AlgoConfig.pe = 0.005;
AlgoConfig.nMM_cut = 2; 


AlgoConfig.verbose = 1;
AlgoConfig.filter_reads = 1;
AlgoConfig.do_filter = 1;
AlgoConfig.do_unique = 1;
AlgoConfig.solve_L2 = 0;
AlgoConfig.filter_included_bacteria = 0;

AlgoConfig.read_type = 'SE';
AlgoConfig.mixture_type = 'Multiplex'; %'RegionByRegion';
AlgoConfig.readLen = 76;
AlgoConfig.barcoded_regions = 1;

AlgoConfig.tol = 5e-7;
AlgoConfig.numIter = 10000;

if 1
    % Build matrices Ad
    if simulate_reads
        dat0 = build_A_matrices(bactData, noisy_reads, AlgoConfig);
    else
        dat0 = build_A_matrices(bactData, experimental_reads, AlgoConfig);
        %         save([resDir '/sample_' num2str(sample_num) '_mapped_reads.mat'],'dat0','-v7.3')
    end
end

% Solve the mixture
[est_frequency,l2_frequency] = solve_iterative_noisy(dat0, AlgoConfig);
%---------------------------------------------------------------------------------




